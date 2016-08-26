#!/bin/bash


#usage:
# bash ~/Documents/MS_Metagene_Tools/RNAseq_sensitivity.sh reference-point plus-bed minus-bed plus-bw minus-bw up down ref_point title ctrl_name kd_names out_dir
#
# computes RNAseq sensitivity of bw using common control
# this requires the bed files to have strand information in field 6 !!
#
# example:
#
# bash ~/Documents/MS_Metagene_Tools/RNAseq_sensitivity.sh reference-point \
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/annotations/gencode_v23_hg19_rep_histones_plusw.bed" \
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/annotations/gencode_v23_hg19_rep_histones_minusw.bed" \
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/*plus.bw" \
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/*minus.bw" \
# 5000 5000 TES histone eGFP "Ars2,Cbp20,Cbp80,NCBP3,Z18" out_dir \
# "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/"
#

ref=$1
plus_bed=$2
minus_bed=$3
plus_bw=$4
minus_bw=$5
up=$6
dn=$7
anchor=$8
title=$9
ctrl_name=${10}
kd_names=${11}
out_dir=${12}

#title="histone"
#anchor="TES"
#plus_bed="/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/annotations/gencode_v23_hg19_rep_histones_plusw.bed"
#minus_bed="/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/annotations/gencode_v23_hg19_rep_histones_minusw.bed"
#bw_dir="/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/"
#out_dir="/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/"


###histone TES case

bash ~/Documents/MS_Metagene_Tools/computeMatrixStranded.sh reference-point \
"${plus_bed}" "${minus_bed}" \
"${plus_bw}" "${minus_bw}" \
$up $dn "${anchor}" "${out_dir}${title}_${anchor}"


## get the min value to be used as pseudocount
min_val=`python ~/Documents/MS_Metagene_Tools/min_value_in_matrix.py "${out_dir}${title}_${anchor}_joined.gz" greaterX=0`
#0.02657

#make the sensitivity matrix wo dropping empty values
python ~/Documents/MS_Metagene_Tools/matrix_to_sensitivity_profile_Effie_like.py "${out_dir}${title}_${anchor}_joined.gz" $ctrl_name $kd_names ${min_val} no

##check it out
plotHeatmap -m "${out_dir}${title}_${anchor}_joined_sensitivity.gz" \
      --sortUsing max --kmeans 1 --colorMap Blues --missingDataColor white \
      --refPointLabel "${anchor}" --plotTitle "${title}" \
      -out "${out_dir}${title}_${anchor}_joined_sensitivity_heatmap.pdf"



