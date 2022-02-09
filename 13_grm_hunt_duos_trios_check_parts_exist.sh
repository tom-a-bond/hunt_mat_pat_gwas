#!/bin/bash/

# check all parts of grms were successfully written

################################################################################
##### variables which require editing:
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory where scripts are kept
gen_path=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/CleanGenotypeData/ # path to plink 1 binary fileset
gen_file=ukb_cal_chr1_22_v2_cleaned # filename stem for plink 1 binary fileset; should have been cleaned already
out_file_sample1=${cohort_name}_sample1_maternal_geno # output file name for sample 1 (marginal maternal gwas)
out_file_sample2=${cohort_name}_sample2_paternal_geno # output file name for sample 2 (marginal paternal gwas)
out_file_sample3=${cohort_name}_sample3_offspring_geno # output file name for sample 3 (marginal offspring gwas)
#out_file_sample4m=${cohort_name}_sample4m_maternal_geno # output file name for sample 4m (conditional trios gwas; maternal ids); n.b. no longer required, because we are only analysing relateds
#out_file_sample4f=${cohort_name}_sample4f_paternal_geno # output file name for sample 4f (conditional trios gwas; paternal ids); n.b. no longer required, because we are only analysing relateds
#out_files=(${out_file_sample1} ${out_file_sample2} ${out_file_sample3} ${out_file_sample4m} ${out_file_sample4f})
out_files=(${out_file_sample1} ${out_file_sample2} ${out_file_sample3}) # n.b.- we now no longer need to calculate the GRMs for samples 4m and 4f, because we are only analysing relateds
out_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # output path
################################################################################

# variables which do not require editing (with possible exception of grm_parts and threads:
grm_parts=20 # number of parts to split the grm calculation into (higher number = lower time/memory requirements)
threads=20 # number of cores available
maf=0.01
source ${script_path}${cohort_name}_params.sh # read in additional variables from ${cohort_name}_params.sh
keep_files=(${keep_file_sample1} ${keep_file_sample2} ${keep_file_sample3} ${keep_file_sample4m} ${keep_file_sample4f})

# check grm parts exist
for((k = 0; k < ${#out_files[@]}; k++)); do
  unset miss
  unset nonmiss
  keep_file=${keep_files[${k}]}
  out_file=${out_files[${k}]}
  echo out_file = ${out_file}
  for((j = 1; j <= ${grm_parts}; j++)); do
    i=$(printf "%0${#grm_parts}d" ${j}) # add leading zeroes because this is how gcta names the output files
    echo i = "${i}"
    if [ -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.bin ]; then
        echo grm file present
    else 
        echo grm file not present "******************************************************"
    fi
    if [ -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.N.bin ]; then
        echo N file present
    else 
        echo N file not present "******************************************************" 
    fi
    if [ -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.id ]; then
        echo id file present
    else 
        echo id file not present "******************************************************" 
    fi
    if [ ! -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.bin ] || [ ! -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.N.bin ] || [ ! -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.id ]; then miss=(${miss[@]} "${j}"); fi
    if [ -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.bin ] && [ -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.N.bin ] && [ -f "${out_path}"${out_file}.part_${grm_parts}_"${i}".grm.id ]; then nonmiss=(${nonmiss[@]} "${j}"); fi
  done
  echo "${miss[@]}" > "${out_path}"${cohort_name}_partitioned_grm_missing_parts_${out_file}.txt
  echo "${nonmiss[@]}" > "${out_path}"${cohort_name}_partitioned_grm_nonmissing_parts_${out_file}.txt
  echo n missing parts = "${#miss[@]}" # should be 0
  echo n nonmissing parts = "${#nonmiss[@]}" # should equal number of grm parts used
done



