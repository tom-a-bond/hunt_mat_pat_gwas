#!/bin/bash
#PBS -l select=1:ncpus=20:mem=20gb
#PBS -l walltime=24:00:00
#PBS -J 1-22
#PBS -j oe
#PBS -k o

# fastgwa gwas, sample 3 (all mothers, fathers and offspring, treated as offspring)

# read in variables
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
sample=3
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory where scripts are kept
source ${script_path}${cohort_name}_params.sh

for((i = 0; i < ${#outcomes[@]}; i++)); do
  outcome=${outcomes[${i}]}
  echo $i
  echo $outcome

  ${gcta} \
		--bgen ${bgen_path}${bgen_file} \
		--sample ${sample_path}${sample_file} \
  	--keep ${keep_path}${keep_file_sample} \
		--maf ${maf} \
		--info ${info} \
		--autosome \
		--threads ${threads} \
		--load-model ${out_path}${cohort_name}_sample_${sample}_${outcome}_model_only.fastGWA \
		--out ${out_path}${cohort_name}_sample_${sample}_${outcome}_chr_${PBS_ARRAY_INDEX}
    # gzip output
    gzip -f ${out_path}${cohort_name}_sample_${sample}_${outcome}_chr_${PBS_ARRAY_INDEX}.fastGWA
done


  
  