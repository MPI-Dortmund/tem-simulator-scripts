# This script takes a PDB and extracts the important lines for the "transf" files of the TEM-Simulator

import re
import os
import argparse

def _main_():
    # Read path to pdb
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=str, help="Path to pdb")
    args = parser.parse_args()
    path_to_pdb = args.pdb


    # Open file and fine target lines
    with open(path_to_pdb, 'r') as textfile:
        filetext = textfile.read()
        textfile.close()
    # https://regex101.com/r/Cucv7i/1/
    # https://regex101.com/r/4tMjXQ/1/
    # https://regex101.com/r/fBIXqf/1/
    biomolecule = "1"
    res = re.findall(
                "BIOMOLECULE: +"+biomolecule+" .*?REMARK \d+ *$",
                filetext,
                re.M | re.S
            )
    for r in res:
        if "BIOMT" in r:
            res = r
            break
    matches = []
    if res:
        matches = re.findall(
            "^.*BIOMT\d+ +\d+ (.*\d) +$",
            re.sub(
                ' +',
                ' ',
                res,
            ),
            re.M
        )

    if res is None or len(matches)==0:
        print("No assembly can be found. No trans file required. Stop")
        import sys
        sys.exit(0)

    # Write target lines to disk
    with open(os.path.splitext(os.path.basename(path_to_pdb))[0]+"_trans.txt", 'w') as f:
        f.write("#"+os.linesep)
        f.write(f"  {len(matches)} 4"+os.linesep)
        f.write("#" + os.linesep)

        for match in matches:
            line = match
            f.write(line+os.linesep)
        f.close()

if __name__ == "__main__":
    _main_()