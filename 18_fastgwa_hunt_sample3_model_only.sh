#!/bin/bash
#PBS -l select=1:ncpus=20:mem=20gb
#PBS -l walltime=04:00:00
#PBS -J 0-26
#PBS -j oe
#PBS -k o

# fit fastgwa model only, loading snps from all chrs (for grammar-gamma approximation)
# this model will then be used to run fastgwa for all chrs individually
# n.b. array sub jobs (i.e. line 4 above) should be numbered 0-(1 minus number of outcomes); there are 27 outcomes in ukb

# read in variables
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
sample=3
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/
source ${script_path}0_${cohort_name}_params.sh

# run fastgwa
#for((i = 0; i < ${#outcomes[@]}; i++)); do
  outcome=${outcomes[${PBS_ARRAY_INDEX}]}
  #echo $i
  echo $outcome
  qcovar_file_name=sample${sample}_${outcome}.qcovar
  covar_file_name=sample${sample}_${outcome}.covar
  # set name for continuous covariate file, dependent on whether we neeed to adjust this phenotype for fasting time
  if [[ ${fasting_outcomes[@]} =~ ${outcome} ]]; then qcovar_file_name=sample${sample}_fasting_${outcome}.qcovar; fi
  ${gcta} \
		--grm-sparse ${grm_path}${grm_file_sample}_sparse_0_05 \
		--mbgen ${imp_path}${cohort_name}_sample3_all_bgen_list.txt \
		--sample ${sample_path}${sample_file} \
		--pheno ${phen_path}sample${sample}_${outcome}.phen \
		--qcovar ${phen_path}${qcovar_file_name} \
		--covar ${phen_path}${covar_file_name} \
		--maf ${maf} \
		--info ${info} \
		--autosome \
		--threads ${threads} \
		--fastGWA-mlm \
		--model-only \
		--out ${out_path}${cohort_name}_sample_${sample}_${outcome}_model_only
#done


  
  