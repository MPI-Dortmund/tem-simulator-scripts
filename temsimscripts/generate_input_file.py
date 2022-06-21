"""
Script to generate input file for TEM-Simulator
"""
import os
import copy
from dataclasses import dataclass
from enum import auto, Enum
from typing import List
import argparse
import numpy as np


class UnknownSourceException(Exception):
    """Get raised when source is unknown"""

class InputMode(Enum):
    tomo = auto()
    particle = auto()

@dataclass
class ParticleSet:
    """
    Represents a set of particles
    """
    particle_type: str
    num_particles: int
    source: str = "pdb"
    pdb_file_in: str = None
    pdb_transf_file_in: str = None
    voxel_size: float = 0.1
    map_file_re_in: str = None
    map_file_format: str = None
    use_imag_pot: str = "yes"
    contrast_re: float = 20
    smoothness: float = 6
    make_positive: str = "yes"
    gen_nxyz: List = None
    particle_coords: str = "random"
    coord_file_in : str = None


@dataclass
class Detector:
    """
    Represents a detector
    """
    det_pix_x: int
    det_pix_y: int
    pixel_size: int  # microns
    gain: int
    use_quantization: str
    dqe: float
    mtf_a: float
    mtf_b: float
    mtf_c: float
    mtf_alpha: float
    mtf_beta: float
    image_file_out: str

    def gen_tem_sim_dict(self, no_noise: bool = False):

        thedect = self
        if no_noise:
            thedect = copy.deepcopy(self)
            thedect.use_quantization = "no"
            splitted = os.path.splitext(thedect.image_file_out)
            thedect.image_file_out = splitted[0] + "_nonoise" + splitted[1]

        return {
                "det_pix_x": thedect.det_pix_x,
                "det_pix_y": thedect.det_pix_y,
                "pixel_size": thedect.pixel_size, # microns
                "gain": thedect.gain,
                "use_quantization": thedect.use_quantization,
                "dqe": thedect.dqe,
                "mtf_a": thedect.mtf_a,
                "mtf_b": thedect.mtf_b,
                "mtf_c": thedect.mtf_c,
                "mtf_alpha": thedect.mtf_alpha,
                "mtf_beta": thedect.mtf_beta,
                "image_file_out": thedect.image_file_out
            }

@dataclass
class Optics:
    magnification: int
    cs: float #mm
    cc: float  #mm
    aperture: int # micrometer
    focal_length: float # mm
    cond_ap_angle: float # milliradian
    gen_defocus: str # yes/no
    defocus_nominal: float # micrometer

    def gen_tem_sim_dict(self) -> dict:
        return {
            "magnification": self.magnification,
            "cs": self.cs,  # mm
            "cc": self.cc,  # mm
            "aperture": self.aperture,  # micrometer
            "focal_length": self.focal_length,  # mm
            "cond_ap_angle": self.cond_ap_angle,  # milliradian
            "gen_defocus": self.gen_defocus,
            "defocus_nominal": self.defocus_nominal # micrometer
        }


default_input = {
    "simulation": {
        "generate_micrographs": "yes",
        "generate_particle_maps": "yes",
        "log_file": "simulator.log"
    },

    "sample": {
        "diameter": 2000, # nm
        "thickness_edge": 125, #nm 125
        "thickness_center": 125, #nm 125
    },

    "geometry": {
        "gen_tilt_data": "yes",
        "ntilts": 61,
        "theta_start": -60, # degree
        "theta_incr": 2, # degree
        "geom_errors": "none"
    },

    "electronbeam": {
        "acc_voltage": 300, #kv
        "energy_spread": 1.3, #eV
        "gen_dose": "yes",
        "total_dose": 15000, # electrons per square nanometer
    },
}


def write_input_file(output_file: str, input_dict: dict):
    output_list = []
    for key in input_dict:

        if key == "particles":
            particles = input_dict[key]
            for particle_id in particles:
                pname = input_dict["particlesets"][particle_id]["particle_type"]
                output_list.append(f"=== particle {pname} ===" + os.linesep)
                for setting in particles[particle_id]:
                    output_list.append(f"{setting} = " + str(particles[particle_id][setting]) + os.linesep)
                if len(particles)>1:
                    output_list.append("" + os.linesep)
        elif key == "particlesets":
            sets = input_dict[key]
            for setid in sets:
                output_list.append("=== particleset ===" + os.linesep)
                for setting in sets[setid]:
                    output_list.append(f"{setting} = " + str(sets[setid][setting]) + os.linesep)
                if len(sets)>1:
                    output_list.append("" + os.linesep)
        elif key == "detectors":
            detectors = input_dict[key]
            for detector_id in detectors:
                output_list.append("=== detector ===" + os.linesep)
                for setting in detectors[detector_id]:
                    output_list.append(f"{setting} = " + str(detectors[detector_id][setting]) + os.linesep)
                if len(detectors)>1:
                    output_list.append("" + os.linesep)
        else:
            output_list.append(f"=== {key} ===" + os.linesep)
            for key2 in input_dict[key]:
                output_list.append(f"{key2} = " + str(input_dict[key][key2]) + os.linesep)
        output_list.append(""+os.linesep)

    with open(output_file, 'w') as f:
        f.write(''.join(output_list))

def add_particle_set(sets: List[ParticleSet], input_dict: dict) -> None:

    def get_particle_set(pset: ParticleSet) -> dict:
        return {
            "particle_type": pset.particle_type,
            "num_particles": pset.num_particles,
            "particle_coords": pset.particle_coords,
            "where": "volume",
            "coord_file_out": pset.particle_type + "_coords.txt",
            "coord_file_in": pset.coord_file_in
        }
    def get_particle(pset: ParticleSet) -> dict:
        particle = {
            "source": pset.source,
            "voxel_size": pset.voxel_size, #nm
            "use_imag_pot": pset.use_imag_pot,
            "map_file_re_out": pset.particle_type + "_map_re.mrc",
            "map_file_im_out": pset.particle_type + "_map_im.mrc"
        }
        if pset.use_imag_pot == "no":
            particle["famp"] = 1.0

        if pset.source == "pdb":
            particle["pdb_file_in"] = pset.pdb_file_in
            if pset.pdb_transf_file_in is not None:
                particle["pdb_transf_file_in"]= pset.pdb_transf_file_in
        if pset.source == "map":
            particle["map_file_format"] = pset.map_file_format
            particle["map_file_re_in"] = pset.map_file_re_in
        if pset.source == "random":
            particle["contrast_re"] = pset.contrast_re
            particle["smoothness"] = pset.smoothness
            particle["make_positive"] = pset.make_positive
            particle["nx"] = pset.gen_nxyz[0]
            particle["ny"] = pset.gen_nxyz[1]
            particle["nz"] = pset.gen_nxyz[2]
        return particle

    if "particlesets" not in input_dict:
        input_dict["particlesets"] = {}
    if "particles" not in input_dict:
        input_dict["particles"] = {}

    for set_id, set in enumerate(sets):

        id = len(input_dict["particlesets"])+1+set_id
        input_dict["particlesets"][id] = get_particle_set(set)
        input_dict["particles"][id] = get_particle(set)


def add_detectors(detectors: List[Detector], input_dict: dict, add_noisefree=True):
    if "detectors" not in input_dict:
        input_dict["detectors"] = {}

    for d_id, detector in enumerate(detectors):

        id = len(input_dict["detectors"])+1+d_id
        input_dict["detectors"][id] = detector.gen_tem_sim_dict()
        if add_noisefree:
            input_dict["detectors"][id+1] = detector.gen_tem_sim_dict(no_noise=True)


def get_source(ext: str) -> str:
    if ext.endswith("pdb"):
        return "pdb"
    elif ext.endswith("mrc"):
        return "map"
    elif ext.endswith("random"):
        return "random"
    else:
        raise UnknownSourceException()

def get_num_part_coords(path: str) -> int:
    with open(path, 'r') as file1:
        lines = file1.readlines()

    num_coords = 0

    for l in lines:
        if l.startswith("#"):
            pass
        else:
            num_coords = num_coords +1
    return num_coords-1

def make_pdb_trans_map(pdb_files: List[str], trans_files: List[str]) -> dict:
    def get_match(pdbid: str, trans_files: List[str]):
        if trans_files is None:
            return None

        for tfile in trans_files:
            if pdbid.lower() in tfile.lower():
                return tfile
        return None
    pdb_trans_map = {}

    for pdb in pdb_files:
        pdb_id = os.path.splitext(os.path.basename(pdb))[0]
        pdb_trans_map[pdb] = get_match(pdb_id, trans_files)

    return pdb_trans_map

def create_parser() -> argparse.ArgumentParser:

    def create_tomo_parser(parser: argparse.ArgumentParser):
        parser.add_argument("--pdbs", type=str, required=True, nargs="+", help="Path to one or more pdbs/maps")
        parser.add_argument("--trans", type=str, default=None, nargs="+", help="Path to one or more transfer files")
        parser.add_argument("--coords", type=str, required=True, nargs="+", help="Path to one or more coordinate files.")
        parser.add_argument("--defocus_lower", type=float, default=5, help="Lower bound of defocus")
        parser.add_argument("--defocus_upper", type=float, default=5, help="Upper bound of defocus")
        parser.add_argument("--output_file", default="input.txt", help="Output file path")
        parser.add_argument("--thickness", default=125, help="Thickness in nm")

    def create_particle_parser(parser: argparse.ArgumentParser):
        parser.add_argument("--pdbs", type=str, required=True, nargs="+", help="Path to one or more pdbs/maps")
        parser.add_argument("--trans", type=str, default=None, nargs="+", help="Path to one or more transfer files")
        parser.add_argument("--output_file", default="input.txt", help="Output file path")


    parser_parent = argparse.ArgumentParser(
        description="Generate input file for tem-simulator"
    )
    subparsers = parser_parent.add_subparsers(help="sub-command help")

    tomo_parser = subparsers.add_parser(
        "tomogram",
        help="Simulate volume",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_tomo_parser(tomo_parser)

    particle_parser = subparsers.add_parser(
        "particle",
        help="Simulate particle maps",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_particle_parser(particle_parser)

    return parser_parent


def get_voxelsize_mrc(pth: str):
    import mrcfile
    with mrcfile.open(pth, permissive=True) as mrc:
        # / 10 to nm, swap x and z axis to be in the same coordinate system
        vsize = np.array([mrc.voxel_size.x, mrc.voxel_size.y, mrc.voxel_size.z])
    return vsize


def _main_():
    import sys
    parser = create_parser()
    args = parser.parse_args()

    pdbs = args.pdbs
    trans = args.trans
    output_file = args.output_file

    if "tomogram" in sys.argv[1]:
        coords = args.coords
        defocus_lower = args.defocus_lower
        defocus_upper = args.defocus_upper
        input_mode = InputMode.tomo
        thickness = args.thickness

    if "particle" in sys.argv[1]:
        input_mode = InputMode.particle
        defocus_lower = 5 # Set it to constants. It is not used anyway in this case
        defocus_upper = 5


    pdb_trans_map = make_pdb_trans_map(pdbs, trans)


    if input_mode == InputMode.particle:
        default_input["simulation"]["generate_micrographs"] = "no"
        tmp_file = os.path.join(os.path.dirname(output_file), '.tmp_coords')
        with open(tmp_file, 'w') as write:
            write.write(' '.join(map(str, [0]*6)))
        coords = [tmp_file] * len(pdbs)

    if input_mode == InputMode.tomo:
        default_input["sample"]["thickness_edge"] = thickness
        default_input["sample"]["thickness_center"] = thickness

    assert len(pdbs) == len(coords) , "Please provide path to coord file for each PDB"

    sets = []
    for i, pdb in enumerate(pdbs):
        pdb_type = os.path.splitext(os.path.basename(pdb))[0]
        ext = os.path.splitext(os.path.basename(pdb))[1]
        source = get_source(ext)
        if source == "pdb":
            particle_set = ParticleSet(
                particle_type=pdb_type,
                source=source,
                num_particles=get_num_part_coords(coords[i]),
                pdb_file_in=os.path.abspath(pdb),
                pdb_transf_file_in=pdb_trans_map[pdb],
                particle_coords="file",
                coord_file_in=coords[i],
                voxel_size=0.1)
            sets.append(particle_set)
        if source == "map":
            vsize_in_angst = get_voxelsize_mrc(os.path.abspath(pdb))
            particle_set = ParticleSet(
                particle_type=pdb_type,
                source=source,
                num_particles=get_num_part_coords(coords[i]),
                map_file_format="mrc",
                use_imag_pot="no",
                map_file_re_in=os.path.abspath(pdb),
                particle_coords="file",
                coord_file_in=coords[i],
                voxel_size=vsize_in_angst[0]/10)
            sets.append(particle_set)
        if source == "random":
            particle_set = ParticleSet(
                particle_type=pdb_type,
                source=source,
                num_particles=get_num_part_coords(coords[i]),
                use_imag_pot="no",
                contrast_re=0.3,
                smoothness=6,
                make_positive="no",
                voxel_size=0.1,
                particle_coords="file",
                coord_file_in=coords[i],
                gen_nxyz=[512,512,512])
            sets.append(particle_set)



    add_particle_set(sets=sets, input_dict=default_input)

    ###############
    # Optics
    ###############
    defocus = defocus_lower + np.random.rand()*(defocus_upper-defocus_lower)
    optics = Optics(
        magnification=4900*2,
        cs=2.7,
        cc=2,
        aperture=80,
        focal_length=3,
        cond_ap_angle=0.1,
        gen_defocus="yes",
        defocus_nominal=defocus
    )
    default_input["optics"] = optics.gen_tem_sim_dict()

    # Detectors
    ## Default
    '''
    default_detector = Detector(det_pix_x=512,
                                det_pix_y=512,
                                pixel_size=15,
                                gain=10,
                                use_quantization="yes",
                                dqe=0.4,
                                mtf_a=0.4,
                                mtf_b=0.7,
                                mtf_c=0.1,
                                mtf_alpha=10,
                                mtf_beta=40,
                                image_file_out="tiltseries.mrc")
    '''

    ## Gatan K3
    k3_detector = Detector(
        det_pix_x=1024,
        det_pix_y=1024,
        pixel_size=5,
        gain=10,
        use_quantization="yes",
        dqe=0.9,
        mtf_a=0.496,
        mtf_b=0.5,
        mtf_c=0,
        mtf_alpha=2.144,
        mtf_beta=2.144,
        image_file_out="tiltseries.mrc")

    add_detectors([k3_detector], default_input, add_noisefree = True)

    write_input_file(output_file, default_input)

if __name__ == "__main__":
    _main_()
