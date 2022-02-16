#!/bin/bash

# prune GRM to create sparse GRM
# this script is not memory or time intensive so it should be possible to run it interactively

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
source ${script_path}0_${cohort_name}_params.sh # read in additional variables from ${cohort_name}_params.sh
keep_files=(${keep_file_sample1} ${keep_file_sample2} ${keep_file_sample3} ${keep_file_sample4m} ${keep_file_sample4f})

# variables which do not require editing (with possible exception of grm_parts and threads:
grm_parts=20 # number of parts to split the grm calculation into (higher number = lower time/memory requirements)
threads=20 # number of cores available
maf=0.01
keep_files=(${keep_file_sample1} ${keep_file_sample2} ${keep_file_sample3} ${keep_file_sample4m} ${keep_file_sample4f})

# combine grms
for((i = 0; i < ${#out_files[@]}; i++)); do
  echo $i
  out_file=${out_files[${i}]}
  ${gcta} \
  		--grm ${out_path}${out_file} \
  		--make-bK-sparse 0.05 \
  		--threads 20 \
  		--out ${out_path}${out_file}_sparse_0_05
done


