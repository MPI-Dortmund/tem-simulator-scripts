#!/usr/bin/env bash

set -e

progs=(
    tsimscripts_gen_filaments.py
    tsimscripts_gen_input.py
    tsimscripts_gen_map.py
    tsimscripts_gen_coords.py
    cp
    submfg
    trimvol
    )

error=false
for prog in ${progs[@]}
do
    if ! which ${prog} &>/dev/null
    then
        echo "Program '${prog}' is missing!"
        error=true
    fi

done
if [[ ${error} == true ]]
then
    exit 1
fi


early_abort=false
defocus_lower=5
defocus_upper=2
nfiducial=10
nvesicle=4
fiducialsize=50
vesiclesize=100
thickness=125
random_seed=
pdbs=()
pdbs_fil=()
settings_fil=()
npdbs=()
nsubs=()
while [[ $# -gt 0 ]]
do
    case $1 in
        --settings_fil)
            shift
            while [[ ! ${1} =~ ^- && ! ${1} == '' ]]
            do
                settings_fil=(${settings_fil[@]} ${1})
                shift
            done
            ;;
        --pdbs_fil)
            shift
            while [[ ! ${1} =~ ^- && ! ${1} == '' ]]
            do
                pdbs_fil=(${pdbs_fil[@]} ${1})
                shift
            done
            ;;
        --nsubs) 
            shift
            while [[ ! ${1} =~ ^- && ! ${1} == '' ]]
            do
                nsubs=(${nsubs[@]} ${1})
                shift
            done
            ;;
        --pdbs)
            shift
            while [[ ! ${1} =~ ^- && ! ${1} == '' ]]
            do
                pdbs=(${pdbs[@]} ${1})
                shift
            done
            ;;
        --npdbs) 
            shift
            while [[ ! ${1} =~ ^- && ! ${1} == '' ]]
            do
                npdbs=(${npdbs[@]} ${1})
                shift
            done
            ;;
        --fiducialsize) fiducialsize=${2}; shift; shift;;
        --vesiclesize) vesiclesize=${2}; shift; shift;;
        --nvesicle) nvesicle=${2}; shift; shift;;
        --nfiducial) nfiducial=${2}; shift; shift;;
        --output) output=${2}; shift; shift;;
        --random_seed) random_seed=${2}; shift; shift;;
        --defocus_lower) defocus_lower=${2}; shift; shift;;
        --defocus_upper) defocus_upper=${2}; shift; shift;;
        --thickness) thickness=${2}; shift; shift;;
        --s1) early_abort=true; shift;;
    esac
done
mkdir -p $(dirname ${output})
output=$(realpath ${output%/}.$(date -Is))
reconstruction_dir=${output}/reconstruction
simulation_dir=${output}/simulation
mkdir -p ${simulation_dir}
ptcls_dir=${output}/ptcls
mkdir -p ${ptcls_dir}
ptcls_fil_dir=${output}/ptcls_filament
mkdir -p ${ptcls_fil_dir}
fiducial_dir=${output}/fiducials
mkdir -p ${fiducial_dir}
coords_dir=${output}/coords
coords_fil_dir=${output}/coords_filament
particle_pdb_dir=${output}/pdbs
mkdir -p ${particle_pdb_dir}
filament_pdb_dir=${output}/pdbs_filament
mkdir -p ${filament_pdb_dir}
filament_pdb_dir_in=${output}/pdbs_filament_in
mkdir -p ${filament_pdb_dir_in}

set -x
if [[ ! -z ${random_seed} ]]
then
    seed_cmd_1="--random_seed ${random_seed}"
    seed_cmd_2="--random_seed $((${random_seed}+1))"
    seed_cmd_3="--random_seed $((${random_seed}+3))"
fi

echo $(date) - Inputs
echo PDBS: ${pdbs[@]}
echo NPDBS: ${npdbs[@]}
echo PDBS FIL: ${pdbs_fil[@]}
echo SETTINGS FIL: ${settings_fil[@]}
echo NSUBS: ${nsubs[@]}

echo $(date) - Generate filaments
cp ${pdbs[@]} ${particle_pdb_dir}
cp ${pdbs_fil[@]} ${filament_pdb_dir_in}
cp ${settings_fil[@]} ${filament_pdb_dir_in}
tsimscripts_gen_filaments.py --pdbs ${filament_pdb_dir_in}/*.pdb --settings ${filament_pdb_dir_in}/*.json --nsubs ${nsubs[@]} ${seed_cmd_3} -o ${filament_pdb_dir}

echo $(date) - Generate PDB to MRC input files
tsimscripts_gen_input.py particle --pdbs ${particle_pdb_dir}/*.pdb --output_file ${ptcls_dir}/input.txt
cd ${ptcls_dir}
echo $(date) - Run PDB to MRC
TEM-simulator input.txt
rm *_map_im.mrc
cd -

echo $(date) - Generate PDB to MRC input files filament
tsimscripts_gen_input.py particle --pdbs ${filament_pdb_dir}/*.pdb  --output_file ${ptcls_fil_dir}/input.txt
cd ${ptcls_fil_dir}
echo $(date) - Run PDB to MRC filament
TEM-simulator input.txt
rm *_map_im.mrc
cd -

echo $(date) - Generate fiducial template
tsimscripts_gen_map.py fiducial -d ${fiducialsize} --apix 1 --value 10000 --output_file ${fiducial_dir}/fiducial.mrc

echo $(date) - Generate vesicle template
tsimscripts_gen_map.py vesicle -d ${vesiclesize} --apix 10 --value 2 --output_file ${fiducial_dir}/vesicle.mrc

echo $(date) - Generate coordinates
tsimscripts_gen_coords.py --pdbs ${filament_pdb_dir}/*.pdb --npdbs 3 --ptcls ${ptcls_fil_dir}/*.mrc ${seed_cmd_1} -o ${coords_fil_dir} --write_occupancy --write_raw_occupancy --tilt_range 15 --max_trials_pos 20 --max_trials_rot 20 --allow_clip --vheight ${thickness}
tsimscripts_gen_coords.py --pdbs ${particle_pdb_dir}/*.pdb --npdbs ${npdbs[@]} --maps ${fiducial_dir}/fiducial.mrc ${fiducial_dir}/vesicle.mrc --nmaps ${nfiducial} ${nvesicle} --ptcls ${ptcls_dir}/*.mrc ${seed_cmd_2} -o ${coords_dir} --write_occupancy --occupancy ${coords_fil_dir}/occupancy_raw.mrc --value_offset ${#pdbs_fil[@]} --vheight ${thickness}


[[ ${early_abort} == true ]] && exit 0
echo $(date) - Generate 3D simulation input files
tsimscripts_gen_input.py tomogram --pdbs  ${particle_pdb_dir}/*.pdb ${fiducial_dir}/fiducial.mrc ${fiducial_dir}/vesicle.mrc ${filament_pdb_dir}/*.pdb --coords ${coords_dir}/*.txt ${coords_fil_dir}/*.txt --defocus_upper ${defocus_upper} --defocus_lower ${defocus_lower} --output_file ${simulation_dir}/input.txt --thickness ${thickness}

cd ${simulation_dir}
echo $(date) - Run 3D simulation
TEM-simulator input.txt
cd -

function generate_reconstruction() {
    local suffix=${1}
    local output_dir=${2}${suffix}
    mkdir -p ${output_dir}

    echo $(date) - Run 3D reconstruction ${suffix}
    ln -rs ${simulation_dir}/tiltseries${suffix}.mrc ${output_dir}/tiltseries${suffix}.mrc

    cat << EOF > ${output_dir}/newst.com
\$setenv IMOD_OUTPUT_FORMAT MRC
\$newstack -StandardInput
AntialiasFilter	-1
InputFile	tiltseries${suffix}.mrc
OutputFile	tiltseries${suffix}_ali.mrc
TransformFile	tiltseries${suffix}.xf
TaperAtFill	1,1
AdjustOrigin	
SizeToOutputInXandY	512,512
OffsetsInXandY	0.0,0.0
#DistortionField	.idf
ImagesAreBinned	1.0
BinByFactor	2
#GradientFile	tiltseries${suffix}.maggrad
\$if (-e ./savework) ./savework
EOF

    cat << EOF > ${output_dir}/tilt.com
\$setenv IMOD_OUTPUT_FORMAT MRC
\$tilt -StandardInput
InputProjections tiltseries${suffix}_ali.mrc
OutputFile tiltseries${suffix}_full_rec.mrc
IMAGEBINNED 2
TILTFILE tiltseries${suffix}.tlt
THICKNESS 400
RADIAL 0.35 0.035
FalloffIsTrueSigma 1
XAXISTILT 0.0
SCALE 0.0 0.2
PERPENDICULAR 
MODE 2
FULLIMAGE 1024 1024
SUBSETSTART 0 0
AdjustOrigin 
ActionIfGPUFails 1,2
XTILTFILE tiltseries${suffix}.xtilt
OFFSET 0.0
SHIFT 0.0 0.0
\$if (-e ./savework) ./savework
EOF

    cat << EOF > ${output_dir}/tiltseries${suffix}.tlt
 -59.99
 -57.99
 -55.98
 -53.98
 -51.98
 -49.98
 -47.98
 -45.99
 -44.00
 -42.00
 -40.01
 -38.02
 -36.03
 -34.04
 -32.05
 -30.05
 -28.06
 -26.05
 -24.04
 -22.03
 -20.01
 -18.00
 -16.00
 -14.00
 -12.00
 -10.00
  -8.00
  -6.00
  -4.00
  -2.00
   0.00
   2.00
   4.00
   6.01
   8.01
  10.01
  12.01
  14.01
  16.00
  18.00
  20.00
  22.00
  24.00
  26.00
  28.00
  30.02
  32.04
  34.06
  36.08
  38.10
  40.09
  42.09
  44.09
  46.09
  48.08
  50.08
  52.08
  54.08
  56.07
  58.06
  60.05
EOF

    cat << EOF > ${output_dir}/tiltseries${suffix}.xf
   0.0006198  -1.0000004   1.0000004   0.0006198       0.389      -0.084
   0.0005976  -0.9999045   0.9999046   0.0005976       0.442      -0.100
   0.0006045  -0.9998404   0.9998405   0.0006045       0.401      -0.082
   0.0006297  -0.9999255   0.9999256   0.0006297       0.352      -0.096
   0.0007459  -0.9999398   0.9999398   0.0007459       0.335      -0.042
   0.0005964  -1.0001166   1.0001166   0.0005964       0.348      -0.079
   0.0005275  -1.0003042   1.0003042   0.0005275       0.305      -0.089
   0.0006434  -1.0001740   1.0001740   0.0006434       0.284      -0.072
   0.0007011  -1.0001374   1.0001374   0.0007011       0.255      -0.018
   0.0005853  -1.0000136   1.0000136   0.0005853       0.102      -0.009
   0.0008688  -1.0000852   1.0000852   0.0008688       0.126       0.033
   0.0007815  -0.9998458   0.9998458   0.0007815       0.037      -0.006
   0.0016303  -0.9977883   0.9977883   0.0016303      -1.139       0.523
   0.0006915  -0.9998064   0.9998064   0.0006915       0.050      -0.058
   0.0005593  -0.9998720   0.9998721   0.0005593      -0.043      -0.090
   0.0006001  -0.9998538   0.9998540   0.0006001      -0.064      -0.077
   0.0004840  -0.9999244   0.9999244   0.0004840      -0.000      -0.084
   0.0004050  -1.0001044   1.0001044   0.0004050      -0.011      -0.069
   0.0004546  -1.0002133   1.0002133   0.0004546      -0.033      -0.122
   0.0004199  -1.0002264   1.0002264   0.0004199       0.033      -0.142
   0.0004891  -1.0002873   1.0002873   0.0004891      -0.005      -0.139
   0.0002907  -1.0003338   1.0003338   0.0002907       0.019      -0.084
   0.0002528  -1.0004265   1.0004265   0.0002528      -0.006      -0.071
   0.0002071  -1.0003651   1.0003651   0.0002071      -0.067      -0.066
   0.0001267  -1.0003604   1.0003604   0.0001267      -0.116      -0.062
   0.0003043  -1.0004730   1.0004730   0.0003043      -0.178      -0.030
   0.0002461  -1.0005348   1.0005348   0.0002461      -0.181      -0.026
   0.0003241  -1.0004606   1.0004606   0.0003241      -0.235      -0.053
   0.0003469  -1.0003299   1.0003299   0.0003469      -0.290      -0.043
   0.0000847  -1.0003934   1.0003934   0.0000847      -0.317       0.056
  -0.0000520  -1.0003105   1.0003105  -0.0000520      -0.253       0.106
  -0.0000537  -1.0003361   1.0003361  -0.0000537      -0.382       0.126
   0.0000069  -1.0003512   1.0003512   0.0000069      -0.409       0.137
  -0.0000510  -1.0004497   1.0004497  -0.0000510      -0.381       0.111
  -0.0000801  -1.0004816   1.0004816  -0.0000801      -0.315       0.100
  -0.0000235  -1.0004929   1.0004929  -0.0000235      -0.326       0.098
   0.0000060  -1.0005163   1.0005163   0.0000060      -0.317       0.104
   0.0000239  -1.0005397   1.0005397   0.0000239      -0.335       0.097
  -0.0000553  -1.0004996   1.0004996  -0.0000553      -0.339       0.097
  -0.0000807  -1.0004818   1.0004818  -0.0000807      -0.346       0.097
  -0.0003872  -1.0005653   1.0005653  -0.0003872      -0.333       0.072
  -0.0003128  -1.0005479   1.0005479  -0.0003128      -0.353       0.049
  -0.0003706  -1.0004785   1.0004785  -0.0003706      -0.408       0.048
  -0.0004034  -1.0005858   1.0005858  -0.0004034      -0.450       0.054
  -0.0003803  -1.0006167   1.0006167  -0.0003803      -0.451       0.044
  -0.0003834  -1.0005541   1.0005541  -0.0003834      -0.480       0.086
  -0.0004030  -1.0004792   1.0004792  -0.0004030      -0.452       0.081
  -0.0003078  -1.0006343   1.0006343  -0.0003078      -0.372       0.040
  -0.0002227  -1.0004218   1.0004218  -0.0002227      -0.366      -0.064
  -0.0002175  -1.0000156   1.0000156  -0.0002175      -0.436       0.023
  -0.0002665  -0.9999346   0.9999346  -0.0002665      -0.436       0.062
  -0.0001109  -1.0000333   1.0000333  -0.0001109      -0.383       0.013
  -0.0002168  -1.0000359   1.0000359  -0.0002168      -0.402      -0.046
  -0.0001949  -0.9998775   0.9998775  -0.0001949      -0.451      -0.050
  -0.0001948  -0.9998790   0.9998790  -0.0001948      -0.431      -0.030
  -0.0003354  -0.9997974   0.9997974  -0.0003354      -0.473      -0.041
  -0.0002641  -0.9997369   0.9997370  -0.0002641      -0.495      -0.032
  -0.0001875  -0.9996973   0.9996973  -0.0001875      -0.494      -0.039
  -0.0001725  -0.9996424   0.9996424  -0.0001725      -0.636      -0.018
  -0.0000709  -0.9996196   0.9996196  -0.0000709      -0.570       0.002
  -0.0002205  -0.9997372   0.9997372  -0.0002205      -0.704      -0.033
EOF

    cat << EOF > ${output_dir}/tiltseries${suffix}.xtilt
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
   0.00
EOF


    cd ${output_dir}
    submfg newst.com
    submfg tilt.com
    trimvol ${output_dir}/tiltseries${suffix}_full_rec.mrc ${output_dir}/tiltseries${suffix}_rec.mrc -RotateX
    cd -
}
generate_reconstruction "_nonoise" "${reconstruction_dir}"
generate_reconstruction "" "${reconstruction_dir}"
echo $(date) - All done
