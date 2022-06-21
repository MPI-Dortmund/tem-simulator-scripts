import argparse
import re
import os

def _main_():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="Path to TEM simulator input file")
    args = parser.parse_args()
    path_to_info = args.input

    # Open file and fine target lines
    with open(path_to_info, 'r') as textfile:
        filetext = textfile.read()

    matches_ntilt = re.findall("^.*ntilts.?=.*$", filetext, re.MULTILINE)
    matches_start = re.findall("^.*theta_start.?=.*$", filetext, re.MULTILINE)
    matches_incr = re.findall("^.*theta_incr.?=.*$", filetext, re.MULTILINE)

    assert len(matches_ntilt) == 1, "At least and only one line for 'ntilts' should be found"
    assert len(matches_start) == 1, "At least and only one line for 'theta_start' should be found"
    assert len(matches_incr) == 1, "At least and only one line for 'theta_incr' be found"

    with open("tiltseries.rawtlt", 'w') as f:
        ntilt = matches_ntilt[0].split()[-1]
        start = matches_start[0].split()[-1]
        incr = matches_incr[0].split()[-1]

        for i in range(int(ntilt)):
            current_tilt = int(start) + i * int(incr)
            f.write(str(current_tilt) + os.linesep)

if __name__ == "__main__":
    _main_()
