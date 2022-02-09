#!bin/bash

# create lists of bgen files for samples 1-3

# read in variables

script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory containing ${cohort_name}_params.sh
source ${script_path}${cohort_name}_params.sh

# sample 1

echo ${sample1_bgen_path}hunt_sample1_mothers_chr1.bgen > ${sample1_bgen_path}hunt_sample1_mothers_bgen_list.txt
for j in {2..22}; do
  echo ${sample1_bgen_path}hunt_sample1_mothers_chr${j}.bgen >> ${sample1_bgen_path}hunt_sample1_mothers_bgen_list.txt
done

# sample 2

echo ${sample2_bgen_path}hunt_sample2_fathers_chr1.bgen > ${sample2_bgen_path}hunt_sample2_fathers_bgen_list.txt
for j in {2..22}; do
  echo ${sample2_bgen_path}hunt_sample2_fathers_chr${j}.bgen >> ${sample2_bgen_path}hunt_sample2_fathers_bgen_list.txt
done

# sample 3

echo ${sample3_bgen_path}chr1.bgen > ${sample3_bgen_path}hunt_sample3_all_bgen_list.txt # n.b. sample 3 file names may require editing 
for j in {2..22}; do
  echo ${sample3_bgen_path}chr${j}.bgen >> ${sample3_bgen_path}hunt_sample3_all_bgen_list.txt # n.b. sample 3 file names may require editing 
done
