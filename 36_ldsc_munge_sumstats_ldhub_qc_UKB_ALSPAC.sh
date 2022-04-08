#!/bin/bash

# munge ldsc sum stats

only_run_missings=true # if true, will only run samples/outcomes which have missing or old ldsc munged sum stats files; if false, will run all samples/outcomes
ma_name=(UKB_ALSPAC) # string to identify this MA
cohorts=(ukb alspac)
samples=(1 2 3)
gwas_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/
phens=(`cat ${gwas_path}input/ph_ukb_alspac.txt`)
script_path=/user/work/tb17613/proj/mat_pat_bmi_mr/scripts/
ldsc_path=/user/home/tb17613/software/ldsc/
ldscores_path=/user/home/tb17613/software/ldsc/downloaded_ld_scores/eur_w_ld_chr/
# array over phens and samples
unset samples_to_run
for((i=0; i<${#phens[@]}; i++)); do
  for((j=0; j<${#samples[@]}; j++)); do
    echo ${phens[${i}]}
    echo ${samples[${j}]}
    phens_to_run=(${phens_to_run[@]} ${phens[${i}]})
    samples_to_run=(${samples_to_run[@]} ${samples[${j}]})
  done
done
# subset runs to only those with missing munged sum stats (if no runs have been previously run, will run all runs
if [[ ${only_run_missings} == true ]]; then
  # check date; n.b. this is a crude date check involving only day
  dates_run_on=({1..31}) # uncomment if we don't need to check whether files are old
  #dates_run_on=(17 24) # uncomment and set (as array variable) to the integer dates (i.e. days of the month) which files should have been run on, if we want to check they were written on the correct days
  unset missing_samples
  unset missing_outcomes
  unset missing_files
  unset missing_array_ids
  for((j = 0; j < ${#samples_to_run[@]}; j++)); do
    sample=${samples_to_run[${j}]}
    sample_upper=${sample^^}
    echo sample = ${sample}
    outcome=${phens_to_run[${j}]}
    echo ${outcome}
    outcome_upper=${outcome^^}
    file_to_check1=${gwas_path}output/${ma_name}/${ma_name}_sample_${sample}_${phen}_all_chr_1.sumstats.gz
    file_date1=($(date -r ${file_to_check1}))
    if [[ ! -e ${file_to_check1} || ! -s ${file_to_check1} || ! ${dates_run_on[@]} =~ ${file_date1[2]} ]]; then
  	   echo ${file_to_check1} is missing or has filesize 0 bytes or is old
  	   missing_samples=(${missing_samples[@]} ${sample})
  	   missing_outcomes=(${missing_outcomes[@]} ${outcome})
  	   #missing_files=(${missing_files[@]} ${file_to_check})
  	   missing_array_ids=(${missing_array_ids[@]} ${j})
  	fi
  done
  n_array_jobs=$((${#missing_array_ids[@]}))
fi
if [[ ${only_run_missings} == false ]]; then
  n_array_jobs=$((${#samples_to_run[@]}))
  missing_array_ids=($(seq 0 $((${n_array_jobs} - 1))))
fi

script_name=${script_path}ldsc_munge_sumstats_ldhub_qc_${ma_name}_JOB_SCRIPT.sh
cat << eof1 > ${script_name}
#!/bin/bash
#SBATCH --job-name=munge
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --array=1-${n_array_jobs}
#SBATCH --mem=10gb

phens_to_run=(${phens_to_run[@]})
samples_to_run=(${samples_to_run[@]})
missing_array_ids=(${missing_array_ids[@]})
array_id=\${missing_array_ids[\${SLURM_ARRAY_TASK_ID}]}
phen=\${phens_to_run[\${array_id}]}
sample=\${samples_to_run[\${array_id}]}
in_file_stem=${gwas_path}output/${ma_name}/${ma_name}_sample_\${sample}_\${phen}_all_chr_1

# remove snps in mhc region, maf <1%, chi sq >80, as per ld hub paper, and remove snps if not present in all cohorts meta-analysed
module load lang/r/4.0.3-bioconductor-gcc # (current default)
Rscript ${script_path}ldsc_ldhub_filters_mp_bmi_mr.R \
\${in_file_stem} \
\${phen} \
0.01 \
80

#module load lang/python/anaconda/3.9.7-2021.12-tensorflow.2.7.0 # did not work (apparently ldsc requires conda 2 not 3)
module load lang/python/anaconda/2.7-2019.10
source activate ldsc

${ldsc_path}munge_sumstats.py \
--sumstats \${in_file_stem}_tmp_ldhub_qc.tbl.gz \
--N-col total_n \
--snp SNP \
--a1 Allele1 \
--a2 Allele2 \
--p P-value \
--signed-sumstats Effect,0 \
--chunksize 500000 \
--out \${in_file_stem}_ldhub_qc \
--merge-alleles ${ldscores_path}w_hm3.snplist

rm \${in_file_stem}_tmp_ldhub_qc.tbl.gz

eof1

chmod +x ${script_name}
sbatch ${script_name}
