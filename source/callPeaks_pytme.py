# written by herman march 2023 as the default tm extract peaks was giving localizations
# that seemed to be incorrect in his data
# modified by joe to improve it

# takes .em probability maps from pytom template matching and turns into a star file
# you need to have run the template matching first
# it will extract the peaks for all the folders in your TMdir
# (template matching directory)
# you should only need to change lines 18, 19 and 20

import os
import argparse
import numpy as np
from skimage.feature import peak_local_max
from tme.matching_utils import load_pickle


def parse_args():
    parser = argparse.ArgumentParser(description="Call Peaks on pytom output.")
    parser.add_argument(
        "--outpath", type=str, required=True, help="Path to the output file."
    )
    parser.add_argument(
        "--TMdir",
        type=str,
        required=True,
        help="Path to the template matching directory.",
    )
    parser.add_argument(
        "--scoresname", type=str, required=True, help="Name of the scores file."
    )
    args = parser.parse_args()
    return args


# read_em from frosina
def read_em(path_to_emfile):
    """Function that reads .em file (form the tom format).
    :param path_to_emfile: file path to the .em file
    :return: header and value
    :rtype: dict, np.array
    """
    with open(path_to_emfile, "r") as f:
        header = dict()
        header["Machine_Coding"] = np.fromfile(f, dtype=np.byte, count=1)
        header["version"] = np.fromfile(f, dtype=np.byte, count=1)
        header["old_param"] = np.fromfile(f, dtype=np.byte, count=1)
        header["data_type_code"] = np.fromfile(f, dtype=np.byte, count=1)
        header["image_dimensions"] = np.fromfile(f, dtype=np.int32, count=3)
        header["the_rest"] = np.fromfile(f, dtype=np.byte, count=496)
        if header["data_type_code"] == 1:
            dtype = np.byte
        elif header["data_type_code"] == 2:
            dtype = np.int16
        elif header["data_type_code"] == 4:
            dtype = np.int32
        elif header["data_type_code"] == 5:
            dtype = np.float32
        elif header["data_type_code"] == 8:
            dtype = np.complex64
        elif header["data_type_code"] == 9:
            dtype = np.double
        else:
            dtype = np.double
            print("dtype was undefined, by default it wil be set to np.double")
        new_image_dim = header["image_dimensions"][::-1]
        header["image_dimensions"] = np.array(new_image_dim)
        value = np.fromfile(f, dtype=dtype)
        value = np.reshape(value, header["image_dimensions"])
        if value.shape[0] == 1 and len(value.shape) == 3:
            value = value[0, :, :]
    return header, value


if __name__ == "__main__":
    args = parse_args()

    starout_path = args.outpath
    TMdir = args.TMdir
    scoresname = args.scoresname

    # Open STAR file and write header
    starout = open(starout_path, "w", newline="")

    # x y z name
    with open(starout_path, "w", newline="") as ofile:
        scores_path = os.path.join(TMdir, scoresname)
        scores = load_pickle(scores_path)[0]

        # Find peaks in volume
        coordinates = peak_local_max(
            scores, min_distance=10, exclude_border=15, num_peaks=4000, p_norm=2
        )
        coordinates = sorted(coordinates, key=lambda x: scores[tuple(x)], reverse=True)
        lines = [
            f"{z}\t{y}\t{x}\t{scores[z,y,x]}\t{TMdir}\n" for z, y, x in coordinates
        ]
        starout.writelines(lines)
    print(f"Finished processing {TMdir}")
