# tem-simulator-scripts

This repository contains convenient scripts for the TEM-Simulator V1.3 (http://tem-simulator.sourceforge.net/)

## Installation
```
conda create -n temscripts -c conda-forge python=3.9 biotite
conda activate temscripts
python setup.py sdist
pip install dist/tem-simulator-scripts-0.1.tar.gz
```

Make also sure that IMOD is installed, as the `tsimscripts_pipe.sh` needs the commands `submfg` and `trimvol` from IMOD.

## Preparation

Download all PDBs you want to simulate. In case PDB is not available (only PDBx or mmcif), download
PDBX or mmcif and convert it to PDB with https://mmcif.pdbj.org/converter/index.php?l=en . Its recommended to download
directely the biological assemblies.

## Overview

Here we give a short overview of the scripts contained in this package:

 - **tsimscripts_pipe.sh**: This is probably the only script a user needs from this collection.
 - **tsimscripts_gen_coords.py**: GEn
 - **tsimscripts_gen_input.py**: This script generates the input file for the TEM-Simulator.
 - tsimscripts_gen_filaments.py
 - tsimscripts_gen_raw_tilt.py
 - tsimscripts_gen_trans_file.py 
 - tsimscripts_gen_map.py

 - tsimscripts_extract.py


## Pipeline all commands

A convenient script is fully automating the whole process incl. reconstruction. 

It requires that imod is installed. The require files for simulating filaments can be found in `resources/filament_files`

A sample command looks like this:

```bash
tsimscripts_pipe.sh --pdbs pdbs/*.pdb --npdbs 150 --output out_sim_tomo_1 --random_seed 10 --pdbs_fil filament_files/*.pdb --settings_fil filament_files/*.json --nsubs 100 --random_seed 10
```

It will generate a tomogram with all pdbs included in the folder `pdbs/` and will add 4 vesicles and 10 fiducials.

## How to simulate a tomogram with particles from a single PDB

You can also run each step individually in case you want to have more control:

0. Install the package with python setup.py install
1. Optional: Generate trans files using `tsimscripts_gen_trans_file.py`. Not necessary when using biological assemblies directly.
2. Optional: Generate fiducial using `tsimscripts_gen_map.py`: I would add 8-10 fiduicals. Patch tracking didnt work well.
3. Generate coordinate `tsimscripts_gen_coords.py`
4. Generate input file using `tsimscripts_gen_input_file.py`
5. Generate raw tilt file with `tsimscripts_gen_raw_tilt.py`
6. Run `TEM-simulator input.txt`
7. Reconstruction with etomo
   * During `Coarse Alignement` I needed to use a `low freq. rolloff sigma` of `0.01`, a `high freq. cutoff sigma` of `0.25` and a `high freq. rolloff sigma` of `0.05`. Moreover I activated the option `No Cosine Stretch`.
   * In `Final Aligned Stack Complete` I set the binning to `2`
   * In `Tomogram Generation In Progress` I set the the `thickness in Z` to `300`
8. Extract the subvolumes with `tsimscripts_extract.py`: If you have used fiducials, check if they are centered in the extracted tomograms.
If they offcentered, repeat the alignement in step 7.

## How to simulate a tomogram with particles from multiple PDBs

Here is a complete example how to simulate a tomogram. It assumes that there are several pdbs in the folder pdbs/.
Only the PDB 6x9q needed a trans_file in this example.

```bash
tsimscripts_gen_map.py fiducial -d 50 --apix 1 -v 10000
mkdir trans_files
tsimscripts_gen_trans_file.py pdbs/6x9q.pdb
mv 6x9q_trans.txt trans_files/
tsimscripts_gen_coords.py --pdbs pdbs/*.pdb --npdbs 20 -o coords/
tsimscripts_gen_input.py --pdbs pdbs/*.pdb fiducial.mrc --trans trans_files/6x9q_trans.txt --coords coords/*.txt --defocus_upper 6 --defocus_lower 7
TEM-simulator input.txt
```