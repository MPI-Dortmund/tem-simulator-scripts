import argparse
import os
import mrcfile
import numpy as np
# last time i run it with --size 37 and --esize 15
# MODEL=model_3; python extract_noise.py --rec ${MODEL}/reconstruction_lp60.mrc --cls ${MODEL}/class_mask.mrc --size 37 --esize 15 --number 100 --out subvolumes_lp60_noise/ --filename ${MODEL}_0XXX


def extract(vol, coords, size):
    pos0_from = (coords[0] - (size - 1) // 2)
    pos0_to = (coords[0] + (size - 1) // 2 + 1)
    pos1_from = (coords[1] - (size - 1) // 2)
    pos1_to = (coords[1] + (size - 1) // 2 + 1)
    pos2_from = (coords[2] - (size - 1) // 2)
    pos2_to = (coords[2] + (size - 1) // 2 + 1)

    sub = vol[pos0_from: pos0_to, pos1_from: pos1_to, pos2_from: pos2_to]

    return sub


def _main_() -> None:
    parser = argparse.ArgumentParser(description='Extracts noise subvolumes')
    parser.add_argument('--rec', type=str, required=True,
                        help='Reconstruction')
    parser.add_argument('--cls', type=str, required=True,
                        help='Class mask')
    parser.add_argument('--size', type=int, required=True,
                        help='Box size.')
    parser.add_argument('--esize', type=int, required=True,
                        help='empty size. the subvolume needs to without any class pixel to get selected')
    parser.add_argument('--number', type=int, required=True,
                        help='Number of noise volumes to extract')
    parser.add_argument('--out', type=str, required=True,
                        help='output path')
    parser.add_argument('--filename', type=str, default="model_x",
                        help='ID ')
    args = parser.parse_args()

    path_reconstruction = args.rec
    path_class_mask = args.cls
    path_output = args.out
    num_noise_vol = args.number
    size = args.size
    esize = args.esize
    filebasename = args.filename
    os.makedirs(path_output, exist_ok=True)

    vol_rec = mrcfile.mmap(path_reconstruction, permissive=True, mode='r').data
    vol_class = mrcfile.mmap(path_class_mask, permissive=True, mode='r').data
    free_positions = np.where(vol_class == 0)

    written = 0
    while written < num_noise_vol:
        # coord_s0 = np.random.randint(size, vol_class.shape[0]-size)
        # coord_s1 = np.random.randint(size, vol_class.shape[1] - size)
        # coord_s2 = np.random.randint(size, vol_class.shape[2] - size)
        rand_position_index = np.random.randint(low=0, high=len(free_positions[0]))
        coord_s0 = free_positions[0][rand_position_index]
        coord_s1 = free_positions[1][rand_position_index]
        coord_s2 = free_positions[2][rand_position_index]

        sub_vol_class = extract(vol_class, coords=(coord_s0, coord_s1, coord_s2), size=esize)

        if np.sum(sub_vol_class) == 0:
            sub_rec = extract(vol_rec, coords=(coord_s0, coord_s1, coord_s2), size=size)
            if sub_rec.shape == (size, size, size):
                sub_rec = -1 * sub_rec
                filename = filebasename + '_' + str(written).zfill(3) + '.mrc'
                with mrcfile.new(os.path.join(path_output, filename)) as newmrc:
                    newmrc.set_data(sub_rec)
                written = written + 1
                print(written)


if __name__ == "__main__":
    _main_()

