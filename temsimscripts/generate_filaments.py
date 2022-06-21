import time
import json
import os
import argparse
import dataclasses as dc

import scipy.spatial.transform as sstr
import pandas
import numpy
from numpy.typing import NDArray
import abc


@dc.dataclass
class FilamentSettings:
    """
    Settings related to a filamentous protein.
    :member rise: Helical rise of the helical symmetry parameters
    :member twist: Helical twist of the helical symmetry parameters
    """

    rise: float
    twist: float


class BaseBend(abc.ABC):
    """
    Base class serving as an interface to the world.
    """

    def __init__(self):
        self.tilt_shifts: list[float] = []
        self.tilt_angles: list[float] = []
        self.z_shifts: list[float] = []
        self.z_angles: list[float] = []
        self.subunits: int = 0

    def set_subunits(self, subunits):
        """
        Set the number of subunits
        :param subunits: Number of subunits for this filament.
        """
        self.subunits = subunits

    def set_tilt_shifts(self, tilt_shifts: list[float]):
        """
        Set the tilt_shifts value
        :param tilt_shifts: List containing the tilt_shifts
        """
        self.tilt_shifts = tilt_shifts

    def set_tilt_angles(self, tilt_angles: list[float]):
        """
        Set the tilt_angles value
        :param tilt_angles: List containing the tilt_angles
        """
        self.tilt_angles = tilt_angles

    def set_z_shifts(self, z_shifts: list[float]):
        """
        Set the z_shifts values
        :param z_shifts: List containing the z_shifts
        """
        self.z_shifts = z_shifts

    def set_z_angles(self, z_angles: list[float]):
        """
        Set the z_angles value
        :param z_angles: List containing the z_angles
        """
        self.z_angles = z_angles

    def get_subunits(self):
        """
        Return the number of subunits
        :return: Number of subunits
        """
        return self.subunits

    def get_params(self) -> tuple[list[float], list[float], list[float], list[float]]:
        """
        Return the stored parameters as a tuple
        :return: Parameters that describe the filament
        """
        return (self.tilt_shifts, self.tilt_angles, self.z_shifts, self.z_angles)


class CoshBend(BaseBend):
    """
    Bend the filament using a cosh function
    """

    def __init__(
        self,
        subunits: int,
        distance: float,
        filament_settings: FilamentSettings,
        twist_offset: float,
    ):
        """
        Initialize the bending machinery
        :param subunits: Number of subunits of the filament
        :param distance: Distance to sample on the cosh function
        :param filament_settings: Dataclass containing the filamentous settings
        :param twist_offset: Offset to the helical twist for a bit of variety
        """
        super().__init__()

        self.scale = 10
        self.set_subunits(subunits)

        # Oversample by 2 to approximate the tilt_angle of the edge subunits
        raw_index: NDArray = numpy.arange(-1, self.get_subunits() + 1)

        indices = self.generate_index(distance)

        z_shifts: NDArray = self.generate_z_shifts(raw_index, filament_settings.rise)
        # Remove the oversampling for the helical rise
        self.set_z_shifts(z_shifts[1:self.get_subunits() + 1])

        tilt_shifts: NDArray = self.generate_tilt_shifts(indices)
        self.set_tilt_shifts(tilt_shifts[1:self.get_subunits() + 1])

        self.generate_z_angle(raw_index[1:self.get_subunits() + 1], filament_settings.twist, twist_offset)
        self.generate_tilt_angle(z_shifts, tilt_shifts)


    def generate_index(self, distance: float):
        """
        Generate valid indices on the target cosh function.
        Each new index is chosen in a way that the arc length stays constant.
        """
        old_x = 0
        indices = [old_x]
        # Adjust the scaling for long filaments
        scale_value = self.scale * numpy.sqrt(self.get_subunits() / 50)
        for _ in range((self.get_subunits() + 2) // 2):
            # Calculation of the next index for a constant arclength of distance:
            # NOTE: This might not be mathematically not 100% correct.
            # integral sqrt(1 + f'(x)^2) dx
            # A scale factor makes this super ugly, so I omitted it

            # Formular:
            # integral cosh(x) dx from old_x to new_x = distance
            # => sinh(new_x) - s*sinh(old_x) = distance
            # => (sinh(new_x) - sinh(old_x)) = distance
            # => sinh(new_x) - sinh(old_x) = distance
            # => sinh(new_x) = distance + sinh(old_x)
            # => new_x = arcsinh( distance + sinh(old_x) )
            new_x = numpy.arcsinh(distance / scale_value + numpy.sinh(old_x))
            indices.append(new_x)
            indices.insert(0, -new_x)
            old_x = new_x
        return indices

    def generate_tilt_shifts(self, indices: NDArray):
        """
        Create tilt_shift values based on the cosh function.
        Remove the offset of the shift to avoid large shifts in the BDB file.
        :param indices: List of indices to use for the cosh function
        :return: Shift values
        """
        shifts = self.scale * numpy.cosh(indices)
        shifts -= numpy.min(shifts)
        return shifts

    def generate_tilt_angle(self, z_shifts: NDArray, tilt_shifts: NDArray):
        """
        Approximate tilt_angle values based on the provided x and y coordinates
        of an unknown function.
        :param z_shifts: x values of the approximate function
        :param tilt_shifts: y values of the approximate function
        """
        diff_y = tilt_shifts[2:2+self.get_subunits()] - tilt_shifts[:self.get_subunits()]
        diff_x = z_shifts[2:2+self.get_subunits()] - z_shifts[:self.get_subunits()]
        approx = diff_y[:diff_x.size] / diff_x
        self.set_tilt_angles(numpy.arctan(approx))

    def generate_z_shifts(self, indices: NDArray, rise: float) -> NDArray:
        """
        Calculate the z_shifts along the helical axis
        :param indices: indices of the individual subunits
        :param rise: Helical rise of the filament
        :return: Actual shift values of the filament for each subunit
        """
        return indices * rise

    def generate_z_angle(self, indices: NDArray, twist: float, twist_offset: float):
        """
        Calculate the helical twist along the helical axis.
        :param indices: indices of the individual subunits.
        :param twist: Helical twist of the filament.
        :param twist_offset: Offset of the helical twist.
        :return: Actual shift values of the filament for each subunit
        """
        self.set_z_angles(indices * twist + twist_offset)


def get_parser() -> argparse.ArgumentParser:
    """
    Return the argument parser for this script.
    :return: Argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdbs", type=str, required=True, nargs="+", help="Path to PDB files"
    )
    parser.add_argument(
        "--settings",
        type=str,
        required=True,
        nargs="+",
        help="Path to filamentous settings files. Need to contain two lines:shift SHIFT and twist TWIST",
    )
    parser.add_argument(
        "--nsubs",
        type=int,
        required=True,
        nargs="+",
        help="Number of subunits per PDB",
    )
    parser.add_argument(
        "--twist_offset",
        type=float,
        help="offset for the curvature",
    )
    parser.add_argument(
        "--mult",
        type=float,
        help="Multiplyer for the curvature function. The higher the value, the stronger the bend.",
    )
    parser.add_argument("--random_seed", type=int, help="Random seed")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output folder")
    return parser


def read_pdb(pdb_file: str) -> pandas.DataFrame:
    """
    Parse the provided pdb file.
    Only the entries ATOM, TER, and HETATM are of interest for the simluation.
    :param pdb_file: Path to the pdb file
    """

    with open(pdb_file) as read:
        atom_data_raw = [
            line.strip().split()
            for line in read.readlines()
            if line.startswith("ATOM")
            or line.startswith("TER")
            or line.startswith("HETATM")
        ]

    # Append a final TER, which separates chains from each other.
    atom_data_raw.append(atom_data_raw[-1][:])
    atom_data_raw[-1][0] = "TER"
    del atom_data_raw[-1][2]
    del atom_data_raw[-1][5:]

    return pandas.DataFrame(
        atom_data_raw,
        columns=[
            "label",
            "idx",
            "atom",
            "acid",
            "Dontknow",
            "residue",
            "x",
            "y",
            "z",
            "b",
            "c",
            "real_atom",
        ],
    )


def read_settings(settings_file: str) -> FilamentSettings:
    """
    Read the filament settings file.
    It needs to be in json format.
    :param settings_file: Path to the json settings file
    :return: FilamentSettings class containg the required information.
    """
    with open(settings_file) as read:
        settings_dict = json.load(read)

    filament_settings = FilamentSettings(**settings_dict)
    return filament_settings


def generate_pdb(atom_data: pandas.DataFrame, bender: BaseBend):
    """
    Generate a new pdb file with the proper length
    :param atom_data: Data frame containin the input pdb information.
    :param bender: Bender class providing the filamentous parameters of the function.
    :return: Data frame containing all pdb information.
    """
    coordinate_array = numpy.array(
        [
            atom_data["x"].astype(float),
            atom_data["y"].astype(float),
            atom_data["z"].astype(float),
        ]
    )
    output_data_frame = pandas.DataFrame(
        index=range(atom_data.shape[0] * bender.get_subunits()),
        columns=atom_data.columns,
    )

    for sub_idx, (tilt_shift, tilt_angle, z_shift, z_angle) in enumerate(
        zip(*bender.get_params())
    ):
        # Rotate the subunit around the global Z axis of the whole box
        rot_mat_phi = sstr.Rotation.from_euler("Z", [z_angle], degrees=True)
        new_array = rot_mat_phi.as_matrix()[0] @ coordinate_array

        # We need to center the subunit before applying the tilt to it
        # to rotate around the center of mass rather then the external
        # coordinate system.
        box_center = numpy.empty(3)
        for box_idx, entry in enumerate(new_array):
            box_min = numpy.min(entry[~numpy.isnan(entry)])
            box_max = numpy.max(entry[~numpy.isnan(entry)])
            center = (box_max + box_min) / 2
            box_center[box_idx] = center
        rot_mat_theta = sstr.Rotation.from_euler("Y", [tilt_angle])
        new_array = rot_mat_theta.as_matrix()[0] @ (new_array - box_center[..., None])
        new_array += box_center[..., None]

        new_array[0] += tilt_shift
        new_array[2] += z_shift

        tmp_data = atom_data.copy()
        tmp_data["x"] = new_array[0].round(3).astype(str)
        tmp_data["y"] = new_array[1].round(3).astype(str)
        tmp_data["z"] = new_array[2].round(3).astype(str)

        output_data_frame.iloc[
            sub_idx * atom_data.shape[0] : (sub_idx + 1) * atom_data.shape[0]
        ] = tmp_data
    return output_data_frame


def write_output(output_data: pandas.DataFrame, output_file: str):
    """
    Write the pdb output file.
    :param output_data: Data frame containing the rows to write
    :param output_file: File to write the output to
    """
    # PDB formatter string... It is what it is
    line_format = "{:<6s}{:>5s}  {:<4s}{:<4s}{:1s}{:>4s}{:>12s}{:>8s}{:>8s}{:>6s}{:>6s}{:>12s}\n"

    out_array = []
    for _, line in output_data.iterrows():
        valid_line = [entry if entry not in (None, "nan") else "" for entry in line]
        if valid_line[0] == "TER":
            valid_line.insert(2, "")
        out_array.append(line_format.format(*valid_line))

    # Remove the final TER and add an END
    del out_array[-1]
    out_array.append("END\n")

    with open(os.path.join(output_file), "w") as out:
        out.writelines(out_array)


def _main_():
    """
    Main function executing the logic :)
    """
    args = get_parser().parse_args()

    if args.random_seed is None:
        random_seed = int(time.time())
    else:
        random_seed = args.random_seed
    numpy.random.seed(random_seed)

    if args.mult is None:
        args.mult = 2 + numpy.random.random() * 3

    if args.twist_offset is None:
        args.twist_offset = numpy.random.random() * 360

    if len(args.pdbs) != len(args.nsubs):
        if len(args.nsubs) == 1:
            args.nsubs = args.nsubs * len(args.pdbs)
        else:
            assert False, (args.pdbs, args.nsubs, "do not have the same length!")
    assert len(args.pdbs) == len(args.nsubs), (
        args.pdbs,
        args.nsubs,
        "do not have the same length!",
    )
    assert len(args.pdbs) == len(args.settings), (
        args.pdbs,
        args.settings,
        "do not have the same length!",
    )

    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, 'random_seed'), 'w') as seed_file:
        seed_file.write(str(random_seed))
    for pdb_file, settings_file, subunits in zip(args.pdbs, args.settings, args.nsubs):
        atom_data = read_pdb(pdb_file)
        filament_settings = read_settings(settings_file)

        bending = CoshBend(subunits, args.mult, filament_settings, args.twist_offset)
        full_pdb = generate_pdb(atom_data, bender=bending)
        name, ext = os.path.splitext(os.path.basename(pdb_file))
        output_file = os.path.join(args.output, f"{name}_{bending.get_subunits()}_{round(args.mult,1)}{ext}")
        write_output(full_pdb, output_file)


if __name__ == "__main__":
    _main_()
