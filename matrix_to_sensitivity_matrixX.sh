#!/bin/bash

help="
usage:
 bash ~/ms_tools/MS_Metagene_Tools/matrix_to_sensitivity_matrixX.sh matrix_name ctrl_name KD_names outfile_name
# \n
# computes RNAseq sensitivity of samples containing "KD_names" relative to "ctrl_name"\n
# NOTE: this version uses a matrix as input
# \n
# example: \n
# \n
# bash ~/ms_tools/MS_Metagene_Tools/RNAseq_sensitivity.sh \ \n
# "matrix.gz" \ \n
# eGFP "Ars2,Cbp20,Cbp80,NCBP3,Z18" out_dir \ \n
# "matrix_rel_eGFP_sensitivity.gz" \n

"

if [ $1 == "-h" ]
then
echo ${help}
exit
fi

matrix=$1
ctrl_name=${2}
kd_names=${3}
out_name=${4}

script_dir="$(dirname -- "$0")/"

echo "using matrix ${matrix}"


## get the min value to be used as pseudocount
echo "finding minimum value to be used as pseudocount"
min_val=$(python "${script_dir}min_value_in_matrix.py" "${matrix}" greaterX=0)
echo "found min val: "$min_val
#0.02657


#make the sensitivity matrix wo dropping empty values
echo "--------------------------------"
echo "--------------------------------"
echo "computing the sensitivity matrix"
echo "--------------------------------"
echo "--------------------------------"
python ${script_dir}matrix_to_sensitivity_profile_Effie_like.py "${matrix}" $ctrl_name $kd_names ${min_val} no no $out_name




