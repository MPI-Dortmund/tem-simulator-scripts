import argparse
from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
from pathlib import Path
from enum import Enum, auto
import mrcfile
from dataclasses import dataclass
import os

'''
Extracts subvolumes and write them to disc
'''

class OutOfVolumeException(Exception):
    ''' out of volume exception'''
    pass

class Origin(Enum):
    CENTER = auto()
    LOWER_LEFT = auto()

@dataclass
class ParticleLike:
    '''
    Represents a particle like object that needs to placed inside a volume
    '''
    identifier: str
    x: float # center x
    y: float # center y
    z: float # center z
    box_size: int
    origin: Origin # Origin of the coorindate system

    def make_coord_LL(self, volume_shape: np.array):
        if self.origin == Origin.CENTER:
            self.x = volume_shape[2] / 2 + self.x
            self.y = (volume_shape[1] / 2 + self.y)
            self.z = volume_shape[0] / 2 + self.z
            self.origin = Origin.LOWER_LEFT

    def extract(self, volume: np.array) -> np.array:
        if self.origin == Origin.CENTER:
            self.make_coord_LL(volume_shape=volume.shape)
        nx1 = (int(self.x) - (self.box_size - 1) // 2)
        nx2 = (int(self.x) + (self.box_size - 1) // 2 + 1)
        ny1 = (int(self.y) - (self.box_size - 1) // 2)
        ny2 = (int(self.y) + (self.box_size - 1) // 2 + 1)
        nz1 = (int(self.z) - (self.box_size - 1) // 2)
        nz2 = (int(self.z) + (self.box_size - 1) // 2 + 1)
        subvol = volume[nz1: nz2, ny1: ny2, nx1: nx2]
        if subvol.shape != (self.box_size, self.box_size, self.box_size):
            raise OutOfVolumeException

        return subvol

def extract_and_write(particles: Dict[str,List[ParticleLike]], volume: np.array, output_dir: Path, prefix: str = "", apix=None) -> None:
    '''Writes coordinates to disk, so that our SHREC evaluation script can read'''
    #if not output_dir.exists():
    output_dir.mkdir(exist_ok=True, parents=True)
    output_dir_coords = output_dir.joinpath("coords/")
    output_dir_coords.mkdir(exist_ok=True, parents=True)


    pos_data = []

    for id in particles:
        coords = []
        for pnum, particle in enumerate(particles[id]):
            try:
                subvol = particle.extract(volume)
                subvol = subvol * -1
            except OutOfVolumeException:
                #print(f"Particle {particle} is out of volume.")
                continue
            num_str = str(pnum).zfill(3)
            coords.append([particle.x, particle.y, particle.z])
            pos_data.append([particle.identifier, particle.x, particle.y, particle.z, np.NaN, np.NaN, np.NaN])
            with mrcfile.new(output_dir.joinpath(prefix + particle.identifier + '_' + num_str + '.mrc')) as newmrc:
                newmrc.set_data(subvol)
                if apix:
                    newmrc.voxel_size = apix
        print("S:", id, np.array(coords).shape)
        np.savetxt(output_dir_coords.joinpath(id+".coords"), np.array(coords), fmt='%f')
    df = pd.DataFrame(pos_data)
    df.to_csv(output_dir.joinpath("particle_positions.txt"), header=False, index=False, na_rep="NaN")

def read_apix(path_vol: Path):

    if path_vol.exists() and path_vol.is_file():
        with mrcfile.mmap(path_vol) as mrc:
            return mrc.voxel_size
    else:
        raise ValueError("Your input volume does not exist.")
    return None

def read_volume(path_vol: Path) -> np.array:
    if path_vol.exists() and path_vol.is_file():
        with mrcfile.open(path_vol) as mrc:
            return mrc.data
    else:
        raise ValueError("Your input volume does not exist.")
    return None

def get_identifier(path: Path):
    import os
    name = os.path.splitext(path.name)[0]
    return name.split("_")[0]

def get_coordinates(path_coord: Path, pixel_size=1.02041) -> List[Tuple]:
    with open(path_coord, 'r') as file1:
        lines = file1.readlines()
    coords = []
    n_hashtags = 0
    for l in lines:
        if n_hashtags > 1:
            lsplitted = l.split(' ')
            x = float(lsplitted[0])/pixel_size
            y = float(lsplitted[1])/pixel_size
            z = float(lsplitted[2])/pixel_size

            # 90 degree anti-clockwise rotation
            help = x
            x = -y
            y = help

            coords.append((x,y,z))

        if l.startswith("#"):
            n_hashtags=n_hashtags+1
    return coords




def read_particles(path_coord: Path, box_size: int) -> Dict[str,List]:
    particles = {}
    if path_coord.exists() and path_coord.is_dir():
        pth_coords = path_coord.glob("*.txt")
        for pth in pth_coords:
            id = get_identifier(pth)
            coords = get_coordinates(pth)

            if id not in particles:
                particles[id] = []
            for coord in coords:
                p = ParticleLike(x=coord[0],y=coord[1],z=coord[2],box_size=box_size, identifier=id, origin=Origin.CENTER)
                particles[id].append(p)
    else:
        raise ValueError("Input path must be a directory")

    return particles

def offset_z(particles: Dict[str,List], zoffset) -> None:
    for id in particles:
        for particle in particles[id]:
            particle.z = particle.z + zoffset



def run(args) -> None:
    """Runs the script"""
    path_coord = Path(args.coord)
    path_vol = Path(args.vol)
    path_out = Path(args.out)
    box_size = args.boxsize
    zoff = args.zoff
    prefix = args.prefix

    volume = read_volume(path_vol)
    apix = read_apix(path_vol)

    particles = read_particles(path_coord, box_size)

    offset_z(particles, zoffset=zoff)

    for id in particles:
        for p in particles[id]:
            p.make_coord_LL(volume_shape=volume.shape)

    extract_and_write(particles=particles, volume=volume,output_dir=path_out, prefix=prefix, apix=apix)

def create_parser() -> argparse.ArgumentParser:
    """Creates parser"""

    parser = argparse.ArgumentParser()
    parser.add_argument("--coord", type=str, required=True, help="Path to folder with coordinaates")
    parser.add_argument("--vol", type=str, required=True, help="Path to tomogram")
    parser.add_argument("--out", type=str, required=True, help="Output path")
    parser.add_argument("--boxsize", type=int, default=37, help="Boxsize in px")
    parser.add_argument("--prefix", type=str, default="", help="Prefix of the filename")
    parser.add_argument("--zoff", type=float, default=0, help="Offset the coordinates in z-direction")


    return parser

def _main_() -> None:
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    _main_()