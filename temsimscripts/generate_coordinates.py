"""
This script will generate the coordinates for a specific tomogram.
"""

import mrcfile
import argparse
import tempfile
from typing import List, Protocol, Tuple, Callable
from pathlib import Path
from dataclasses import dataclass
import os
from skimage.transform import rescale, resize, downscale_local_mean
from biotite import structure as struc
import biotite.structure.io as strucio
import biotite.database.rcsb as rcsb
import biotite.structure.io.pdbx as pdbx
import time
import numpy as np
import json

import scipy.ndimage as snd
import scipy.spatial.transform as sstr


from numpy.typing import NDArray
from tqdm import tqdm



class WrongDTypeException(Exception):
    """Expection for wrong dtype"""
    ...


class IncompatibleShapesException(Exception):
    """Exception when shapes do not match"""
    ...

@dataclass
class OccupancyVolume:
    """
    Saves where placed particles already occupy the volume.
    """

    volume: NDArray

    def size(self) -> int:
        return self.volume.size

    def shape(self) -> tuple:
        return self.volume.shape

    def is_occupied(self, vol_slice: NDArray, mask: NDArray) -> bool:
        '''
        Check if position in already occupied.
        '''

        return np.any(self.volume[vol_slice][mask] >= 1)

    def occupy(self, vol_slice: NDArray, volume: NDArray, value: int) -> bool:
        '''
        Occupies a part of the volume
        '''
        mask = volume > 0
        if self.is_occupied(vol_slice, mask):
            return False

        self.volume[vol_slice][mask] = value
        return True


@dataclass
class ParticleLike:
    """
    Represents a one or many particle like objects that needs to placed inside a volume
    """

    name: str
    diameter: float  # in angstrom
    n_occurence: int  # Number of particles of this type
    positions: NDArray = None  # Positions of one or many particles in a coordinate system with lower left corner as origin
    positions_tsim: NDArray = None  # Positions in the TEM simulator format
    source: str = None  # Can be path to pdb or a constant string
    mrc_file: str = None  # Can be path to mrc or None
    mrc_vol: NDArray = None
    fill: bool = True

    def get_mrc_vol(self):
        if self.mrc_vol is None and self.mrc_file is not None:
            with mrcfile.open(self.mrc_file, permissive=True) as mrc:
                # / 10 to nm, swap x and z axis to be in the same coordinate system
                vsize = np.array([mrc.voxel_size.x,mrc.voxel_size.y,mrc.voxel_size.z])
                scale = np.array(vsize)/10
                mrc_data = resize(mrc.data, (mrc.data.shape[0] // (1/scale[0]), mrc.data.shape[1] // (1/scale[1]), mrc.data.shape[2] // (1/scale[2])),
                                       anti_aliasing=True, order=5).swapaxes(0, 2)

                #mrc_data = snd.zoom(mrc.data, 1 / 10, order=5).swapaxes(0, 2)
            box_size = int(1.5 * max(mrc_data.shape))
            pad_list = []
            for size in mrc_data.shape:
                pad = box_size - size
                min_pad = pad // 2
                max_pad = pad // 2 + (box_size != size) * abs(pad % 2)
                pad_list.append((min_pad, max_pad))
            self.mrc_vol = np.pad(mrc_data, pad_list, mode='minimum')
        return self.mrc_vol

    def rotate_volume(self, orient: Tuple[float, float, float]):

        """
        Rotates 3D image around image center
        INPUTS
          orient: list of Euler angles (phi,the,psi)
        OUTPUT
          arrayR: rotated 3D numpy array
        by E. Moebel, 2020
        """
        phi = orient[0]
        the = orient[1]
        psi = orient[2]

        # create meshgrid
        dim = self.mrc_vol.shape
        ax = np.arange(dim[0])
        ay = np.arange(dim[1])
        az = np.arange(dim[2])
        coords = np.meshgrid(ax, ay, az)

        # stack the meshgrid to position vectors, center them around 0 by substracting dim/2
        xyz = np.vstack([coords[0].reshape(-1) - float(dim[0]) / 2,  # x coordinate, centered
                         coords[1].reshape(-1) - float(dim[1]) / 2,  # y coordinate, centered
                         coords[2].reshape(-1) - float(dim[2]) / 2])  # z coordinate, centered

        # create transformation matrix
        r = sstr.Rotation.from_euler('ZYZ', [phi, the, psi], degrees=True)
        mat = r.as_matrix()

        # apply transformation
        transformed_xyz = np.dot(mat, xyz)

        # extract coordinates
        x = transformed_xyz[0, :] + float(dim[0]) / 2
        y = transformed_xyz[1, :] + float(dim[1]) / 2
        z = transformed_xyz[2, :] + float(dim[2]) / 2

        x = x.reshape((dim[1], dim[0], dim[2]))
        y = y.reshape((dim[1], dim[0], dim[2]))
        z = z.reshape((dim[1], dim[0], dim[2]))  # reason for strange ordering: see next line

        # the coordinate system seems to be strange, it has to be ordered like this
        new_xyz = [y, x, z]

        # sample
        return snd.map_coordinates(self.mrc_vol, new_xyz, order=1)


class Binarizer(Protocol):
    """Protocol for binarizr. These methods convert a particle into a binary mask."""

    def __call__(
        self, particle_type: ParticleLike, particle_position: Tuple
    ) -> NDArray:
        """
        Converts a particle into a binary volume.
        :param particle_type: Particle type for that a binary volume must be generated
        :param particle_position: 6 entry tuple (x,y,z,phi,theta,psi)
        :return: Binary volume with dtype bool.
        """
        ...

def get_max_diameter_volume(volume: NDArray) -> float:
    from skimage.measure import label, regionprops
    '''
    blobs_labels = label(volume, background=0)
    print("Labeling done")
    props = regionprops(blobs_labels)
    print("done")
    diameter = props[0].feret_diameter_max # WAY TO SLOW
    print("Diameter", diameter) 
    '''

    return np.max(volume.shape)#diameter


def get_max_diamater_by_pdb(pdb_id: str) -> float:
    """
    Method to calculate the diameter of a PDB/mmtf
    :param pdb_id:
    :return:
    """

    path = rcsb.fetch(pdb_id, "pdbx", tempfile.gettempdir())
    return get_max_diameter_pdb(Path(path))


def get_max_diameter_pdb(path: Path) -> float:
    """
    Method to calculate the diameter of a PDB/mmtf
    :param path: Path to pdb
    :return: diameter in Angstrom
    """

    atom_array = None
    if path.suffix == ".pdbx":
        # need to catch cases where I've assemblis
        pdbx_file = pdbx.PDBxFile.read(str(path))

        atom_array = pdbx.get_assembly(pdbx_file, assembly_id="1", model=1)
        assem = pdbx.list_assemblies(pdbx_file)
        print(f"For {path.name} assembly 1 is used: {assem['1']}")

    if atom_array is None:
        atom_array = strucio.load_structure(str(path))
    # Remove all non-amino acids
    try:
        atom_array = atom_array.get_array(0)
    except AttributeError:
        pass  # ignore
    atom_array = atom_array[struc.filter_amino_acids(atom_array)]
    # coord is a NumPy array
    coord = atom_array.coord
    min_xyz = np.amin(coord, axis=0)
    max_xyz = np.amax(coord, axis=0)
    return np.linalg.norm(max_xyz - min_xyz)


def create_bin_sphere(
    particle_type: ParticleLike, particle_position: Tuple[float, float, float, float, float, float], **_
) -> NDArray:
    """
    Creates binary sphere
    """
    diameter_in_nm = particle_type.diameter / 10  # to nm
    radius = diameter_in_nm / 2  # radius
    volume_size = int(radius * 2)+1
    coords = np.ogrid[: volume_size, : volume_size, : volume_size]
    distance = np.sqrt(
        (coords[0] - volume_size/2) ** 2
        + (coords[1] - volume_size/2) ** 2
        + (coords[2] - volume_size/2) ** 2
    )
    return distance <= radius

def calc_quantil(volume: NDArray, quantil = 0.1) -> float:
    '''

    :param volume: volume to find threshold for
    :param quantil: Relative amount of pixels (> min volume) below the threshold
    :return:
    '''
    flat = volume[volume > volume.min()].flatten()
    return np.quantile(flat, quantil)

def dilate(volume: NDArray, ndilations: int):
    volume = np.pad(volume, [[ndilations * 2 + 1] * 2] * 3)
    if ndilations > 0:
        volume = snd.binary_dilation(volume, iterations=ndilations)
    return volume

def fill(volume: NDArray):
    return snd.binary_fill_holes(volume)



def create_bin_mrc(
    particle_type: ParticleLike, particle_position: Tuple[float, float, float, float, float, float]
) -> NDArray:
    """
    Creates binary mask of a volume
    """

    volume = particle_type.get_mrc_vol()

    threshold = np.max(volume) - (np.max(volume)-np.min(volume))*0.9 #calc_quantil(volume, quantil=0.93)
    rotated_volume = particle_type.rotate_volume(particle_position[3:])

    return rotated_volume > threshold

def get_binary_func(particle: ParticleLike) -> Binarizer:
    if particle.mrc_file:
        return create_bin_mrc
    else:
        return create_bin_sphere

def get_ndilations(particle: ParticleLike, ndilations: int) -> int:
    if particle.mrc_file:
        return ndilations
    else:
        return 0

def generate_positions(
    particles: List[ParticleLike],
    occupancy: OccupancyVolume,
    ndilations: int,
    allow_clip: bool,
    tilt_range: float,
    max_trials_rot: int = 5,
    max_trials_pos: int = 100,
    value_offset: int = 0,
) -> Tuple[List[ParticleLike], dict]:
    """
    Generates positions for all particles within the volume. This method assumes a flat sample.
    :param particles: Particles to place
    :param occupancy: Occupancy volume
    :param ndilations: number of dilations to on each binarized particle
    :param max_trials_pos: Maximum number of not successfull placing trials before giving up.
    :param max_trials_rot: Maximum number of rotation to try before giving up.
    :param value_offset: offset for particle id
    :return: Tuple of two elements: List of ParticleLike with updated positions, dictonary that maps particle id to filename
    """

    def gen_orientation(tilt_range: float) -> Tuple[float, float, float]:

        phi = np.random.rand() * 360  # phi
        thetha = np.random.rand() * 2 * tilt_range  # thetha
        psi = np.random.rand() * 360  # psi
        thetha += 90 - tilt_range

        return (phi, thetha, psi)

    def minimal_volume(volshape: NDArray, ptcl: NDArray, pos: Tuple[float,float,float], allow_clip):
        volume_idx = np.where(ptcl > 0)
        current_mask = np.ones(len(volume_idx[0]), dtype=bool)
        for idx in range(len(ptcl.shape)):
            ptcl_coords = int(pos[idx]) - ptcl.shape[idx] // 2 + volume_idx[idx]
            mask = (0 <= ptcl_coords) & (ptcl_coords < volshape[idx])
            # Skip if clipped
            if not mask.all() and not allow_clip:
                return None, None
            current_mask = current_mask & mask

        if np.count_nonzero(current_mask) / len(volume_idx[0]) < 0.4:
            return None, None

        vol_coords = []
        ptc_coords = []
        for idx in range(len(ptcl.shape)):
            ptcl_coords = int(pos[idx]) - ptcl.shape[idx] // 2 + volume_idx[idx]
            vol_coords.append(ptcl_coords[current_mask])
            ptc_coords.append(volume_idx[idx][current_mask])

        try:
            return (
                np.s_[
                min(vol_coords[0]):max(vol_coords[0]) + 1,
                min(vol_coords[1]):max(vol_coords[1]) + 1,
                min(vol_coords[2]):max(vol_coords[2]) + 1
                ],
                np.s_[
                min(ptc_coords[0]):max(ptc_coords[0]) + 1,
                min(ptc_coords[1]):max(ptc_coords[1]) + 1,
                min(ptc_coords[2]):max(ptc_coords[2]) + 1
                ]
            )
        except:
            return None, None

    def gen_bin(particle, current_ndilations, phi, thetha, psi):
        binarized_particle_volume = binarize_func(
            particle_type=particle,
            particle_position=[0,0,0,phi,thetha,psi],
        )
        if particle.fill:
            binarized_particle_volume = dilate(binarized_particle_volume, current_ndilations)
            binarized_particle_volume = fill(binarized_particle_volume)
        return binarized_particle_volume


    #total_size = occupancy.size()
    occ_id_particle_map = {}
    for particle_id, particle in enumerate(particles, value_offset + 1):
        if particle_id not in occ_id_particle_map:
            occ_id_particle_map[particle_id] = particle.name
        #particle_diameter = particle.diameter / 10  # diameter in nm
        number_occ_particle = particle.n_occurence
        pcoords = np.zeros(shape=(number_occ_particle, 6))
        particles_placed = 0


        pbar = tqdm(total=number_occ_particle, desc=f"Simulating {particle.name}")
        phi, thetha, psi = gen_orientation(tilt_range)
        binarize_func = get_binary_func(particle)
        current_ndilations = get_ndilations(particle, ndilations)
        binarized_particle_volume = gen_bin(particle, current_ndilations, phi, thetha, psi)
        unsucessfull_tries = 0
        unsucessfull_tries_rot = 0
        while particles_placed < number_occ_particle:
            if unsucessfull_tries >= max_trials_pos:
                unsucessfull_tries = 0
                unsucessfull_tries_rot += 1
                if unsucessfull_tries_rot >= max_trials_rot:
                    print(f"Too many unsuccessfull position trials. Giving up placing {particle.name}. Continue with next protein.")
                    break
                phi, thetha, psi = gen_orientation(tilt_range)
                binarized_particle_volume = gen_bin(particle, current_ndilations, phi, thetha, psi)
            center_0 = np.random.rand() * occupancy.shape()[0]
            center_1 = np.random.rand() * occupancy.shape()[1]
            center_2 = np.random.rand() * occupancy.shape()[2]

            pos = (center_0, center_1, center_2, phi, thetha, psi)


            #binarized_particle_volume = dilate_and_fill(binarized_particle_volume, ndilations=current_ndilations)




            # Returns the valid slices for volume and particle
            try:
                vol_slice, ptcl_slice = minimal_volume(occupancy.shape(), binarized_particle_volume, pos[:3], allow_clip)
            except ZeroDivisionError:
                print("Zero devision error!!")
                print(particle)
                import sys
                sys.exit()
            if vol_slice is None:
                unsucessfull_tries = unsucessfull_tries + 1
                continue
            elif occupancy.occupy(vol_slice, binarized_particle_volume[ptcl_slice], particle_id):
                pcoords[particles_placed, :] = pos
                particles_placed = particles_placed + 1
                pbar.update(1)

                # Generate a new orientation
                phi, thetha, psi = gen_orientation(tilt_range)
                binarized_particle_volume = gen_bin(particle, current_ndilations, phi, thetha, psi)
                unsucessfull_tries = 0
                unsucessfull_tries_rot = 0
            else:
                unsucessfull_tries = unsucessfull_tries + 1
                continue

        # print(f"Filled: {np.round(np.count_nonzero(occupancy.volume) / total_size * 100, 2)}%")
        particle.positions = pcoords
    return particles, occ_id_particle_map


def create_parser() -> argparse.ArgumentParser:
    """
    Create the argument parser
    :return: Argumentparser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdbs", type=str, nargs="+", help="Path to PDB/mmtf files"
    )
    parser.add_argument(
        "--npdbs",
        type=int,
        nargs="+",
        help="Number of particles per PDB",
    )

    parser.add_argument(
        "--maps", type=str, nargs="+", help="Path to mrc maps"
    )
    parser.add_argument(
        "--nmaps",
        type=int,
        nargs="+",
        help="Number of particles per map",
    )
    parser.add_argument(
        "--ndilations", type=int, default=1, help="Number of dilations to perform on the input volume if --ptcls is provided. Will be ignored otherwise."
    )
    parser.add_argument(
        "--vdiameter", type=int, default=512, help="Volume diameter in nm"
    )
    parser.add_argument("--vheight", type=int, default=125, help="Volume height in nm")
    parser.add_argument("--imodheight", type=int, default=200, help="Imod output Volume height in nm")
    parser.add_argument("--random_seed", type=int, help="Random seed")
    parser.add_argument(
        "--ptcls", type=str, nargs="+", help="Path to mrc ptcl files"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output folder"
    )
    parser.add_argument(
        "--occupancy", type=str, help="Path to a previous occupancy map"
    )
    parser.add_argument(
        "--tilt_range", type=float, default=90, help="Maximum tilt 90+-tilt_range"
    )
    parser.add_argument(
        "--max_trials_pos", type=int, default=500, help="Maximum trials to find a position match"
    )
    parser.add_argument(
        "--max_trials_rot", type=int, default=1, help="Maximum trials to find a rotation match"
    )
    parser.add_argument(
        "--value_offset", type=int, default=0, help="Offset for the assigned values in the occupancy map."
    )
    parser.add_argument(
        "--write_occupancy", action='store_true', help="Write the occupancy map to disc"
    )
    parser.add_argument(
        "--write_raw_occupancy", action='store_true', help="Write the raw occupancy map to disc"
    )
    parser.add_argument(
        "--allow_clip", action='store_true', help="Allow clipping of volumes"
    )

    return parser


def create_particles(
    pdbs: List[str], npdbs: List[int], maps: List[str], nmaps: List[int], mrc_file: List[str]
) -> List[ParticleLike]:
    """
    Create the particles for that positions should be generated
    :param pdbs: List of PDB paths
    :param npdbs: Number of particles that should be generated per PDB path
    :return: List of ParticleLike objects
    """
    particles = []
    if pdbs:
        for i, pdb_path in enumerate(pdbs):
            print("Create particles for ", pdb_path)
            pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
            if os.path.exists(pdb_path):
                diameter = get_max_diameter_pdb(Path(pdb_path))
            else:
                diameter = get_max_diamater_by_pdb(pdb_path)
            particle = ParticleLike(
                name=pdb_name, diameter=diameter, n_occurence=npdbs[i], source=pdb_path, mrc_file=mrc_file[i]
            )
            particles.append(particle)

    if maps:
        for i, maps_path in enumerate(maps):
            print("Create particles for ", maps_path)
            maps_name = os.path.splitext(os.path.basename(maps_path))[0]
            particle = ParticleLike(
                name=maps_name, diameter=None, n_occurence=nmaps[i], source=maps_path, mrc_file=maps_path, fill=False
            )
            bin_func = get_binary_func(particle)
            bin_vol = bin_func(particle,(0,0,0,0,0,0))
            particle.diameter = get_max_diameter_volume(bin_vol)

            particles.append(particle)

    return particles


def convert_coords_to_tem_simulator_format(
    particles: List[ParticleLike], volume_diameter: int, volume_height: int
) -> List[ParticleLike]:
    """
    Method to convert the particle coord to tem simulator format
    :param particles: List of particles
    :param volume_diameter: Volume diameter in nm
    :param volume_height: Volume height in nm
    :return:
    """
    for particle in particles:
        pos = particle.positions
        tsim_pos = np.zeros(shape=(particle.n_occurence, 6))  # 6: x,y,z, phi, theta, psi

        volume_center = np.array(
            [volume_diameter / 2, volume_diameter / 2, volume_height / 2]
        )
        for i in range(len(pos)):
            shifted_pos = pos[i, :3] - volume_center
            tsim_pos[i, 0:3] = shifted_pos
            tsim_pos[i, 3:] = pos[i, 3:]  # phi, theta, psi
        particle.positions_tsim = tsim_pos
    return particles


def write_to_disk(particles: List[ParticleLike], output_dir: str, random_seed: int) -> None:
    """
    Write coordiantes in a TEM Simulator format to disk
    :param particles: Lsit of particles
    :param output_dir: Output dir
    :return: None
    """
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'random_seed'), 'w') as seed_file:
        seed_file.write(str(random_seed))

    for particle in particles:
        # Write target lines to disk
        with open(os.path.join(output_dir, particle.name + "_coords.txt"), "w") as coordinate_file:
            coordinate_file.write("#" + os.linesep)
            coordinate_file.write(f"  {particle.n_occurence} 6" + os.linesep)
            coordinate_file.write("#" + os.linesep)
            for i in range(particle.n_occurence):
                items = particle.positions_tsim[i, :].tolist()
                items_str = [str(v) for v in items]
                line = " ".join(items_str)
                coordinate_file.write(line + os.linesep)


def run(args) -> None:
    """
    Patch everything together and run it :-)
    """

    # Get all arguments
    pdbs = args.pdbs
    npdbs = args.npdbs
    maps = args.maps
    nmaps = args.nmaps
    ptcls = args.ptcls

    vdiameter = args.vdiameter
    vheight = args.vheight
    output = args.output
    ndilations = args.ndilations

    allow_clip = args.allow_clip
    tilt_range = args.tilt_range
    max_trials_rot = args.max_trials_rot
    max_trials_pos = args.max_trials_pos
    value_offset = args.value_offset

    if args.random_seed is None:
        random_seed = int(time.time())
    else:
        random_seed = args.random_seed
    np.random.seed(random_seed)

    if maps is None and pdbs is None:
        assert False, "Maps or pdbs needs to be provided!"

    if pdbs is not None:
        assert (
            len(pdbs) == len(npdbs) or len(npdbs) == 1
        ), "npdbs and pdbs must have the same number of elements."

    if npdbs and len(npdbs) == 1:
        npdbs = npdbs * len(pdbs)

    if pdbs and ptcls is not None:
        assert len(pdbs) == len(ptcls), "ptcls and pdbs must have the same number."
    else:
        ptcls = [None] * len(pdbs)

    if maps is not None:
        assert len(maps) == len(nmaps), "maps and nmaps have to have the same length"

    if args.occupancy is None:
        occupancy = OccupancyVolume(
            volume=np.zeros(
                (vdiameter, vdiameter, vheight), dtype=np.int8
            )
        )
    else:
        with mrcfile.open(args.occupancy, permissive=True) as mrc:
            occupancy = OccupancyVolume(volume=mrc.data)
            occupancy.volume.setflags(write=1)


    # Create ParticleLike
    print("create particles")
    particles = create_particles(
        pdbs=pdbs,
        npdbs=npdbs,
        maps=maps,
        nmaps=nmaps,
        mrc_file=ptcls
    )
    print("sort particles")
    particles.sort(key=lambda x:  np.sum(
        fill(dilate(get_binary_func(x)(x, particle_position=(0,0,0,0,0,0)), ndilations=1))
    ), reverse=True)

    # Generate positions
    print("generate positions")
    particles, occ_map = generate_positions(
        particles=particles,
        occupancy=occupancy,
        ndilations=ndilations,
        allow_clip=allow_clip,
        tilt_range=tilt_range,
        max_trials_rot=max_trials_rot,
        max_trials_pos=max_trials_pos,
        value_offset=value_offset,
    )

    # Convert positions
    particles = convert_coords_to_tem_simulator_format(
        particles, volume_diameter=vdiameter, volume_height=vheight
    )

    # Write to disk

    # TODO: Add calculation of pad factor
    write_to_disk(particles=particles, output_dir=output, random_seed=random_seed)


    f = open(os.path.join(output, 'occu_map.json'), "w")
    json.dump(occ_map, f, indent=4)
    f.close()

    if args.write_occupancy:
        with mrcfile.new(os.path.join(output, 'occupancy.mrc'), overwrite=True) as mrc:
            pad = args.imodheight - vheight
            pad_min = pad // 2
            pad_max = pad // 2 + (pad % 2)
            swap_volume = np.rot90(
                np.pad(
                    occupancy.volume.swapaxes(0, 1).swapaxes(0, 2).swapaxes(1, 2),
                    ((pad_min, pad_max), (0, 0), (0, 0))
                ),
                k=3,
                axes=(1, 2)
            )
            mrc.set_data(swap_volume)
            mrc.voxel_size = 10.20
            mrc.header.origin.z = -1024

    if args.write_raw_occupancy:
        with mrcfile.new(os.path.join(output, 'occupancy_raw.mrc'), overwrite=True) as mrc:
            mrc.set_data(occupancy.volume)
            mrc.voxel_size = 10.20
            mrc.header.origin.z = -1024


def _main_():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    _main_()
