#!/bin/bash

# run genomic sem for all phenos by chr
# bluepebble

only_run_missings=false # if true, will only run samples/outcomes which have missing gneomic sem files; if false, will run all samples/outcomes; n/b. not yet implemented
ma_name=(UKB_ALSPAC) # string to identify this MA
chrs=({1..22})
gwas_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/
phens=(`cat ${gwas_path}input/ph_ukb_alspac.txt`)
script_path=/user/work/tb17613/proj/mat_pat_bmi_mr/scripts/
out_folder=genomic_sem/
snp_list_path=/user/work/tb17613/proj/mat_pat_bmi_mr/data/
snp_list=ukb_alspac_hunt_imputed_snps_maf_0_01_r2_0_8.txt.gz
ref_panel_path=/user/home/tb17613/software/genomic_sem/
ref_panel=reference.1000G.maf.0.005.txt
threads=20

# array over phens and samples
unset chrs_to_run
for((i=0; i<${#phens[@]}; i++)); do
  for((j=0; j<${#chrs[@]}; j++)); do
    echo ${phens[${i}]}
    echo ${chrs[${j}]}
    phens_to_run=(${phens_to_run[@]} ${phens[${i}]})
    chrs_to_run=(${chrs_to_run[@]} ${chrs[${j}]})
  done
done
# subset runs to only those with genomic sem output (if no runs have been previously run, will run all runs); n.b. not yet implemented
if [[ ${only_run_missings} == true ]]; then
  # check date; n.b. this is a crude date check involving only day
  dates_run_on=({1..31}) # uncomment if we don't need to check whether files are old
  #dates_run_on=(17 24) # uncomment and set (as array variable) to the integer dates (i.e. days of the month) which files should have been run on, if we want to check they were written on the correct days
  unset missing_chrs
  unset missing_outcomes
  unset missing_files
  unset missing_array_ids
  for((j = 0; j < ${#chrs_to_run[@]}; j++)); do
    chr=${chrs_to_run[${j}]}
    echo chr ${chr}
    outcome=${phens_to_run[${j}]}
    echo ${outcome}
    file_to_check1=${gwas_path}output/${ma_name}/genomic_sem/${ma_name}_genomic_sem_results_mat_${outcome}_chr${chr}.txt.gz
    file_date1=($(date -r ${file_to_check1}))
    file_to_check2=${gwas_path}output/${ma_name}/genomic_sem/${ma_name}_genomic_sem_results_pat_${outcome}_chr${chr}.txt.gz
    file_date1=($(date -r ${file_to_check2}))
    file_to_check3=${gwas_path}output/${ma_name}/genomic_sem/${ma_name}_genomic_sem_results_off_${outcome}_chr${chr}.txt.gz
    file_date1=($(date -r ${file_to_check3}))
    if [[ ! -e ${file_to_check1} || ! -s ${file_to_check1} || ! ${dates_run_on[@]} =~ ${file_date1[2]} || ! -e ${file_to_check2} || ! -s ${file_to_check2} || ! ${dates_run_on[@]} =~ ${file_date2[2]} || ! -e ${file_to_check3} || ! -s ${file_to_check3} || ! ${dates_run_on[@]} =~ ${file_date3[2]} ]]; then
  	   echo ${file_to_check1} or ${file_to_check2} or ${file_to_check3} is missing or has filesize 0 bytes or is old
  	   missing_chrs=(${missing_chrs[@]} ${chr})
  	   missing_outcomes=(${missing_outcomes[@]} ${outcome})
  	   missing_array_ids=(${missing_array_ids[@]} ${j})
  	fi
  done
  n_array_jobs=$((${#missing_array_ids[@]}))
fi
if [[ ${only_run_missings} == false ]]; then
  n_array_jobs=$((${#chrs_to_run[@]}))
  missing_array_ids=($(seq 0 $((${n_array_jobs} - 1))))
fi

script_name=${script_path}genomic_sem_${ma_name}_JOB_SCRIPT.sh
cat << eof1 > ${script_name}
#!/bin/bash
#SBATCH --job-name=genomic_sem
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${threads}
#SBATCH --time=48:00:00
#SBATCH --array=0-$((${n_array_jobs} - 1))
#SBATCH --mem=100gb
#SBATCH -o ${gwas_path}output/${ma_name}/${out_folder}genomic_sem%A_%a.out

# run wlm; save conditional gwas sum stats

phens_to_run=(${phens_to_run[@]})
chrs_to_run=(${chrs_to_run[@]})
missing_array_ids=(${missing_array_ids[@]})
array_id=\${missing_array_ids[\${SLURM_ARRAY_TASK_ID}]}
phen=\${phens_to_run[\${array_id}]}
chr=\${chrs_to_run[\${array_id}]}

module load lang/r/4.0.3-bioconductor-gcc # (current default)
Rscript ${script_path}run_genomic_sem_trios.R \
${gwas_path}output/${ma_name}/ \
${ma_name}_\${phen}_wlm.tsv.gz \
${out_folder} \
${threads} \
\${phen} \
${snp_list_path}${snp_list} \
${ref_panel_path}${ref_panel} \
${ma_name} \
\${chr}

eof1

chmod +x ${script_name}
sbatch ${script_name}

