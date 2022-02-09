#!/bin/bash

# check fastgwa output for samples 1-3, all phenotypes

################################################################################
##### variables which require editing
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
out_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/results/gwas/duos_trios/
outcomes=(bmi chol crp dbp glu hba1c hdl ldl sbp tg whr bw gest lnbmi lnwhr lnglu lnhba1c lnhdl lntg lncrp)
################################################################################
samples=(1 2 3)

# check fastgwa output
unset missing_samples
unset missing_outcomes
unset missing_files
unset missing_chromosomes
for((j = 0; j < ${#samples[@]}; j++)); do
  sample=${samples[${j}]}
  echo sample = ${sample}
  for((i = 0; i < ${#outcomes[@]}; i++)); do
    outcome=${outcomes[${i}]}
    echo $outcome
    for((k = 1; k <= 22; k++)); do
      file_to_check=${out_path}${cohort_name}_sample_${sample}_${outcome}_chr_${k}.fastGWA.gz
      if [[ ! -e ${file_to_check} || ! -s ${file_to_check} ]]; then
    	   echo ${file_to_check} is missing or has filesize 0 bytes
    	   missing_samples=(${missing_samples[@]} ${sample})
    	   missing_outcomes=(${missing_outcomes[@]} ${outcome})
    	   missing_chromosomes=(${missing_chromosomes[@]} ${k})
    	   missing_files=(${missing_files[@]} ${file_to_check})
    	fi
  	done
  done
done
echo ${missing_files[@]}
echo ${missing_samples[@]}
echo ${missing_outcomes[@]}
echo ${missing_chromosomes[@]}
echo ${#missing_files[@]} # should be 0
echo ${#missing_samples[@]} # should be 0
echo ${#missing_outcomes[@]} # should be 0
echo ${#missing_chromosomes[@]} # should be 0



