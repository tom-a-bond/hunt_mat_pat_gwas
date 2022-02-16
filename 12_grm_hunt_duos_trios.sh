#!/bin/bash
#PBS -l select=1:ncpus=20:mem=2gb
#PBS -l walltime=24:00:00
#PBS -J 1-20
#PBS -j oe
#PBS -k o

# calculate GRM, for samples 1m 2m 3m 4m and 4f
# partition GRM calculation into 20 array sub-jobs (see https://cnsgenomics.com/software/gcta/#MakingaGRM)
# formula for calculating memory requirements (see https://cnsgenomics.com/software/gcta/#MakingaGRM) (n * (n + 1) / 2 * 12) / 1024^3  + 0.5; n.b. I ended up doubling this because it ran out of memory

################################################################################
##### variables which require editing:
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory where scripts are kept
gen_path=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/CleanGenotypeData/ # path to plink 1 binary fileset
gen_file_sample1=${cohort_name}_sample1_maternal_geno # filename stem for plink 1 binary fileset for sample 1; should have been cleaned already
gen_file_sample2=${cohort_name}_sample2_paternal_geno # filename stem for plink 1 binary fileset for sample 2; should have been cleaned already
gen_file_sample3=foo # filename stem for plink 1 binary fileset for sample 3 (can be the whole cohort, i.e. no subsetting required); should have been cleaned already
gen_file_sample4m=${cohort_name}_sample4m_maternal_geno # filename stem for plink 1 binary fileset for sample 4m; should have been cleaned already
gen_file_sample4f=${cohort_name}_sample4f_paternal_geno # filename stem for plink 1 binary fileset for sample 4f; should have been cleaned already
out_file_sample1=${cohort_name}_sample1_maternal_geno # output file name for sample 1 (marginal maternal gwas)
out_file_sample2=${cohort_name}_sample2_paternal_geno # output file name for sample 2 (marginal paternal gwas)
out_file_sample3=${cohort_name}_sample3_offspring_geno # output file name for sample 3 (marginal offspring gwas)
#out_file_sample4m=${cohort_name}_sample4m_maternal_geno # output file name for sample 4m (conditional trios gwas; maternal ids); n.b. no longer required, because we are only analysing relateds
#out_file_sample4f=${cohort_name}_sample4f_paternal_geno # output file name for sample 4f (conditional trios gwas; paternal ids); n.b. no longer required, because we are only analysing relateds
#out_files=(${out_file_sample1} ${out_file_sample2} ${out_file_sample3} ${out_file_sample4m} ${out_file_sample4f})
out_files=(${out_file_sample1} ${out_file_sample2} ${out_file_sample3}) # n.b.- we now no longer need to calculate the GRMs for samples 4m and 4f, because we are only analysing relateds
grm_out_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # output path
################################################################################

# variables which do not require editing (with possible exception of grm_parts and threads:
grm_parts=20 # number of parts to split the grm calculation into (higher number = lower time/memory requirements)
threads=20 # number of cores available
maf=0.01
source ${script_path}0_${cohort_name}_params.sh # read in additional variables from ${cohort_name}_params.sh
keep_files=(${keep_file_sample1} ${keep_file_sample2} ${keep_file_sample3} ${keep_file_sample4m} ${keep_file_sample4f})
gen_files=(${gen_file_sample1} ${gen_file_sample2} ${gen_file_sample3} ${gen_file_sample4m} ${gen_file_sample4f})

# calculate grm part: loop over 3 samples (n.b.- we now no longer need to calculate the GRMs for samples 4m and 4f, because we are only analysing relateds)
for((i = 0; i < ${#out_files[@]}; i++)); do
  keep_file=${keep_files[${i}]}
  gen_file=${gen_files[${i}]}
  out_file=${out_files[${i}]}
  ${gcta} \
  		--bfile ${gen_path}${gen_file} \
  		--maf ${maf} \
  		--autosome \
  		--make-grm-part ${grm_parts} ${PBS_ARRAY_INDEX} \
  		--thread-num ${threads} \
  		--out ${grm_out_path}${out_file}
done
