import argparse
import mrcfile
import numpy as np
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import sys

MEMBRANE_THICKNESS = 6

def average_filter(input: np.array, mask_size=(3,3,3)):

    weights = np.ones(shape=mask_size)
    weights = weights/(mask_size[0]*mask_size[1]*mask_size[2])
    template = ndimage.convolve(input,weights=weights,mode='reflect')
    return template

def gen_sphere(diameter: int, value : float, size = None) -> np.array:
    '''
    Generates the fiducial template
    :param diameter: Diameter in pixel
    :param value: Value that is used to fill it
    '''
    if size is None:
        w = int(diameter * 1.3)
        if w % 2 == 0:
            w = w + 1
    else:
        w = size

    radius = diameter / 2  # radius
    volume_size = w ##int(radius * 2) + 1
    coords = np.ogrid[: volume_size, : volume_size, : volume_size]
    distance = np.sqrt(
        (coords[0] - volume_size / 2) ** 2
        + (coords[1] - volume_size / 2) ** 2
        + (coords[2] - volume_size / 2) ** 2
    )
    template = (distance <= radius).astype(np.float32)
    template[template==1] = value

    return template

def generate_fiducial(diameter: int, value: float):
    fid = gen_sphere(diameter, value)
    fid = gaussian_filter(fid, sigma=0.5)
    return fid

def generate_vesikel(diameter: int, value: float):
    size = int(diameter * 1.3)
    if size % 2 == 0:
        size = size + 1

    fid = gen_sphere(diameter, value, size=size)
    fid2 = gen_sphere(diameter-MEMBRANE_THICKNESS, value, size=size)

    sphere = fid - fid2
    fid = None
    fid2 = None

    fid_inner_1 = gen_sphere(diameter - 2, value, size=size)
    fid_inner_2 = gen_sphere(diameter - (MEMBRANE_THICKNESS-2), value, size=size)
    sphere2 = fid_inner_1 - fid_inner_2
    del fid_inner_1
    del fid_inner_2
    # sphire2 = gaussian_filter(sphere2, sigma=1)

    fid = sphere - sphere2
    del sphere
    del sphere2

    return fid

def create_parser() -> argparse.ArgumentParser:

    def create_fiducial_parser(parser: argparse.ArgumentParser):
        parser.add_argument("-d", "--diameter", type=int, required=True, help="Diameter of the fiducial in pixel")
        parser.add_argument("--apix", type=float, required=True, help="Pixelsize of the map")
        parser.add_argument("-v", "--value", type=float, default=1000, help="Electron potential")
        parser.add_argument("-o", "--output_file", default="fiducial.mrc", help="Output fiducial file.")

    def create_vesicle_parser(parser: argparse.ArgumentParser):
        parser.add_argument("-d", "--diameter", type=int, required=True, help="Diameter in pixel")
        parser.add_argument("--apix", type=float, required=True, help="Pixelsize of the template")
        parser.add_argument("-v", "--value", type=float, default=2, help="Electron potential")
        parser.add_argument("-o", "--output_file", default="vesicle.mrc", help="Output vesicle file.")

    parser_parent = argparse.ArgumentParser(
        description="Generate map templates for TEM-Simulator"
    )
    subparsers = parser_parent.add_subparsers(help="sub-command help")

    fiducial_parser = subparsers.add_parser(
        "fiducial",
        help="Create fiducial map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_fiducial_parser(fiducial_parser)

    vesicle_parser = subparsers.add_parser(
        "vesicle",
        help="Create vesicle map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_vesicle_parser(vesicle_parser)

    return parser_parent


def _main_():
    # Read path to pdb
    parser = create_parser()

    args = parser.parse_args()

    diameter = args.diameter
    apix = args.apix
    value = args.value
    output_file = args.output_file
    if "fiducial" in sys.argv[1]:
        fid = generate_fiducial(diameter, value)
    elif "vesicle" in sys.argv[1]:
        fid = generate_vesikel(diameter, value)


    with mrcfile.new(output_file, overwrite=True) as mrc:
        mrc.set_data(fid)
        mrc._set_voxel_size(apix,apix,apix)






if __name__ == "__main__":
    _main_()
