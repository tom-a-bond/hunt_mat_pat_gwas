
# fastgwa gwas, sample 4m (trios, maternal genotype conditional analyses)
# merge log files from individual snps to one file, for each phenotype

cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
sample=4f # 4m or 4f
out_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/results/gwas/duos_trios/ # output path (directory containing folders written by fastgwa_hunt_sample4m.sh)
outcomes=(bmi chol crp dbp glu hba1c hdl ldl sbp tg whr bw gest lnbmi lnwhr lnglu lnhba1c lnhdl lntg lncrp)

# merge log files for all snps, write output to parent directory
for((j = 0; j < ${#outcomes[@]}; j++)); do
  outcome=${outcomes[${j}]}
  echo ${outcome}
  run=${cohort_name}_${sample}_${outcome}
  cd ${out_path}${run}
  cat *.log > ../${cohort_name}_sample_${sample}_${outcome}_all_snps.log
done

  