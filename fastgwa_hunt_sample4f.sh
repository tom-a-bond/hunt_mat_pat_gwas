#!/bin/bash
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=24:00:00
#PBS -J 0-26
#PBS -j oe
#PBS -k o

# fastgwa gwas, sample 4m (trios, maternal genotype conditional analyses)
# n.b. the line in the header "#PBS -J 0-26" should be 0-(1 minus number of outcomes). In UKB we have 27 outcomes, including logged variables as separate outcomes

cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
sample=4f # 4m or 4f
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory where scripts are kept
source ${script_path}${cohort_name}_params.sh
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
    --keep ${keep_path}${keep_file_sample4f} \
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
    --keep ${keep_path}${keep_file_sample4m} \
    --snps ${snp} \
    --recode A \
    --hard-call-threshold 0.499999 \
    --out ${out_path}${run}/${cohort_name}_trios_mothers_${snp}
  fi

  # offspring
  ${plink2} \
  --bgen ${sample4_bgen_path}${sample4o_bgen_file} ref-first \
  --sample ${sample4_sample_path}${sample4o_sample_file} \
  --keep ${keep_path}${keep_file_sample4o} \
  --snps ${snp} \
  --recode A \
  --hard-call-threshold 0.499999 \
  --out ${out_path}${run}/${cohort_name}_trios_offspring_${snp}

  # set name for continuous covariate file, dependent on whether we need to adjust the current phenotype for fasting time
  covar_file_name=sample${sample}_${outcome}.covar
  qcovar_file_name=sample${sample}_${outcome}.qcovar
  if [[ ${fasting_outcomes[@]} =~ ${outcome} ]]; then qcovar_file_name=sample${sample}_fasting_${outcome}.qcovar; fi
  
  # append parental and offspring genotype to continous covariate file
  Rscript ${script_path}sample4_append_geno_covars_file_hunt.R \
  ${phen_path}${qcovar_file_name} \
  ${out_path}${run}/${cohort_name}_trios_mothers_${snp}.raw \
  ${out_path}${run}/${cohort_name}_trios_fathers_${snp}.raw \
  ${out_path}${run}/${cohort_name}_trios_offspring_${snp}.raw \
  ${keep_path}${linkage_file_sample4} \
  ${sample} \
  ${snp} \
  ${out_path}${run}/ \
  ${phen_path}sample${sample}_${outcome}.phen

  # if there is cryptic relatedness in the sample, run fastgwa LMM
  if [[ ${relateds} == true ]]; then
    # set variables for chr 22 bgen (which is required to estimate model parameters; we use sample 1 or sample 2 bgens as appropriate)
    if [[ ${sample} == '4m' ]]; then
      model_bgen_path=${sample1_bgen_path}; model_bgen_file=hunt_sample1_mothers_chr22.bgen
      model_sample_path=${sample1_sample_path}; model_sample_file=hunt_sample1_mothers_chr22.sample
    fi
    if [[ ${sample} == '4f' ]]; then
      model_bgen_path=${sample2_bgen_path}; model_bgen_file=hunt_sample2_fathers_chr22.bgen
      model_sample_path=${sample2_sample_path}; model_sample_file=hunt_sample2_fathers_chr22.sample
    fi
    #  fastGWA model only; using snps from chr22 to estimate model parameters
    ${gcta} \
  		--grm-sparse ${grm_path}${grm_file_sample}_sparse_0_05 \
  		--bgen ${model_bgen_path}${model_bgen_file}
  		--sample ${model_sample_path}${model_sample_file} \
  		--keep ${keep_path}${keep_file_sample} \
  		--pheno ${phen_path}sample${sample}_${outcome}.phen.resid \
  		--joint-covar \
  		--qcovar ${out_path}${run}/qcovars_${snp}.qcovar \
  		--maf ${maf} \
  		--info ${info} \
  		--autosome \
  		--threads ${threads} \
  		--fastGWA-mlm \
  		--model-only \
  		--out ${out_path}${run}/${run}_${snp}
    # run fastgwa (load model)
    ${gcta} \
  		--bgen ${sample4_bgen_path}${sample_bgen_file} \
  		--sample ${sample_path}${sample_file} \
  		--keep ${keep_path}${keep_file_sample} \
  		--extract ${out_path}${run}/snp.txt \
  		--nofilter \
  		--threads ${threads} \
  		--load-model ${out_path}${run}/${run}_${snp} \
  		--out ${out_path}${cohort_name}_sample_${sample}_${outcome}_${snp}
  fi

  # if there is not relatedness in the sample, run fastgwa linear regression
  if [[ ${relateds} == false ]]; then
    echo Error: for HUNT, please set relateds=true
  fi
  
  # clean up
  if [[ ${relateds} == true ]]; then rm ${out_path}${run}/${run}_${snp}.fastGWA.mdl.id; rm ${out_path}${run}/${run}_${snp}.fastGWA.mdl.bin; fi
  rm ${out_path}${run}/qcovars_${snp}.qcovar
  rm ${out_path}${run}/snp.txt
  rm ${out_path}${run}/${cohort_name}_trios_*_${snp}.*
done



  
  