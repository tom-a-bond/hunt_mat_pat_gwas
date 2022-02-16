#!/bin/bash
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=24:00:00
#PBS -J 0-26
#PBS -j oe
#PBS -k o

# plink gwas, sample 4m (trios, maternal genotype conditional analyses)
# n.b. the line in the header "#PBS -J 0-26" should be 0-(1 minus number of outcomes). In UKB we have 27 outcomes, including logged variables as separate outcomes

cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
sample=4m # 4m or 4f
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory where scripts are kept
source ${script_path}0_${cohort_name}_params.sh
threads=1 # manually override the variable loaded by the params script
snps=(`cat ${snp_list_path}${snp_list}`)
outcome=${outcomes[${PBS_ARRAY_INDEX}]}

load_r_module # function defined in ${script_path}${cohort_name}_params.sh
# loop over snps; write temporary covars file, including the relevant parental genotype and offspring genotype
for((i = 0; i < ${#snps[@]}; i++)); do
  snp=${snps[${i}]}
  echo ${snp}
  run=${cohort_name}_${sample}_${outcome}
  mkdir -p ${out_path}${run} # temporary directory to store intermediate files
  echo ${snp} > ${out_path}${run}/snp.txt

  # extract each parent's and offspring's genotype
  # using --hard-call-threshold 0.499999 to ensure that we have no missing genotypes
  # fathers
  if [[ ${sample} == '4m' ]]; then
    ${plink2} \
    --bgen ${sample4_bgen_path}${sample4f_bgen_file} ref-first \
    --sample ${sample4_sample_path}${sample4f_sample_file} \
    --keep ${keep_path}${keep_file_sample4f_unrelated} \
    --snps ${snp} \
    --recode A \
    --hard-call-threshold 0.499999 \
    --out ${out_path}${run}/${cohort_name}_trios_fathers_${snp}
  fi
  
  # mothers
  if [[ ${sample} == '4f' ]]; then
    ${plink2} \
    --bgen ${sample4_bgen_path}${sample4m_bgen_file} ref-first \
    --sample ${sample4_sample_path}${sample4m_sample_file} \
    --keep ${keep_path}${keep_file_sample4m_unrelated} \
    --snps ${snp} \
    --recode A \
    --hard-call-threshold 0.499999 \
    --out ${out_path}${run}/${cohort_name}_trios_mothers_${snp}
  fi

  # offspring
  ${plink2} \
  --bgen ${sample4_bgen_path}${sample4o_bgen_file} ref-first \
  --sample ${sample4_sample_path}${sample4o_sample_file} \
  --keep ${keep_path}${keep_file_sample4o_unrelated} \
  --snps ${snp} \
  --recode A \
  --hard-call-threshold 0.499999 \
  --out ${out_path}${run}/${cohort_name}_trios_offspring_${snp}

  # set name for continuous covariate file (we will not adjust for fasting time in HUNT)
  covar_file_name=sample${sample}_${outcome}_unrelated.covar

  # append parental and offspring genotype to continous covariate file
  Rscript ${script_path}0_sample4_append_geno_covars_file_hunt.R \
  ${phen_path}${covar_file_name} \
  ${out_path}${run}/${cohort_name}_trios_mothers_${snp}.raw \
  ${out_path}${run}/${cohort_name}_trios_fathers_${snp}.raw \
  ${out_path}${run}/${cohort_name}_trios_offspring_${snp}.raw \
  ${keep_path}${linkage_file_sample4_unrelated} \
  ${sample} \
  ${snp} \
  ${out_path}${run}/ \
  ${phen_path}sample${sample}_${outcome}_unrelated.phen

  # association testing for unrelateds only, using plink
  ${plink2} \
		--bgen ${sample4_bgen_path}${sample_bgen_file} ref-first \
		--sample ${sample4_bgen_path}${sample_sample_file} \
		--keep ${keep_path}${keep_file_sample} \
		--extract ${out_path}${run}/snp.txt \
		--pheno ${phen_path}sample${sample}_${outcome}.phen.unrelated.resid \
		--covar ${out_path}${run}/covars_${snp}_unrelated.covar \
		--covar-variance-standardize \
		--threads ${threads} \
		--glm hide-covar \
		--out ${out_path}${run}/${cohort_name}_sample_${sample}_${outcome}_${snp}

  # clean up
  rm ${out_path}${run}/covars_${snp}_unrelated.qcovar
  rm ${out_path}${run}/snp.txt
  rm ${phen_path}sample${sample}_${outcome}.phen.unrelated.resid
  rm ${out_path}${run}/${cohort_name}_trios_*_${snp}.*
done



  
  