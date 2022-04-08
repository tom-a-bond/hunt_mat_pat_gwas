#!/bin/bash

# run on bluepebble
# n.b. don't forget to set/change "dates_run_on" below

# set vars
only_run_missings=true # if true, will only run samples/outcomes which have missing or old html reports; if false, will run all samples/outcomes
ma_name=(UKB_ALSPAC) # string to identify this MA
cohorts=(ukb alspac)
samples=(1 2 3 4m 4f)
gwas_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/
phens=(`cat ${gwas_path}input/ph_ukb_alspac.txt`)
script_path=/user/work/tb17613/proj/mat_pat_bmi_mr/scripts/
qc_dir=${gwas_path}qc/
dir_output=${gwas_path}qc/${ma_name}
dir_references=/user/work/tb17613/software/gwasinspector/reference
header_translations=alt_headers_edit_fastgwa_chr_pos_a1_a2.txt
allele_ref_std=HRC_r1-1_GRCh37.sqlite
allele_ref_std_population=COMMON
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
# subset runs to only those with missing html reports (if no runs have been previously run, will run all runs
if [[ ${only_run_missings} == true ]]; then
  # check date; n.b. this is a crude date check involving only day
  #dates_run_on=({1..31}) # uncomment if we don't need to check whether files are old
  dates_run_on=(17 24) # uncomment and set (as array variable) to the integer dates (i.e. days of the month) which files should have been run on, if we want to check they were written on the correct days
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
    if [[ ${sample} == '1' || ${sample} == '2' || ${sample} == '3' ]]; then snps_chr=chr; fi
    if [[ ${sample} == '4m' || ${sample} == '4f' ]]; then snps_chr=snps; fi
    file_to_check1=${dir_output}/SAMPLE_${sample_upper}_${outcome_upper}_ukb_sample_${sample}_${outcome}_all_${snps_chr}.fastGWA_report.html
    file_date1=($(date -r ${file_to_check1}))
    file_to_check2=${dir_output}/SAMPLE_${sample_upper}_${outcome_upper}_alspac_sample_${sample}_${outcome}_all_${snps_chr}.fastGWA_report.html
    file_date2=($(date -r ${file_to_check2}))
    file_to_check3=${dir_output}/SAMPLE_${sample_upper}_${outcome_upper}_UKB_ALSPAC_sample_${sample}_${outcome}_all_${snps_chr}_1.tbl_report.html
    file_date3=($(date -r ${file_to_check3}))
    file_to_check4=${dir_output}/SAMPLE_${sample_upper}_${outcome_upper}_Report.html
    file_date4=($(date -r ${file_to_check4}))
    if [[ ! -e ${file_to_check1} || ! -s ${file_to_check1} || ! ${dates_run_on[@]} =~ ${file_date1[2]} || ! -e ${file_to_check2} || ! -s ${file_to_check2} || ! ${dates_run_on[@]} =~ ${file_date2[2]} ||  ! -e ${file_to_check3} || ! -s ${file_to_check3} || ! ${dates_run_on[@]} =~ ${file_date3[2]} || ! -e ${file_to_check4} || ! -s ${file_to_check4} || ! ${dates_run_on[@]} =~ ${file_date4[2]} ]]; then
  	   echo ${file_to_check1} or ${file_to_check2} or ${file_to_check3} or ${file_to_check4} is missing or has filesize 0 bytes or is old
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

script_name=${script_path}qc_${ma_name}_all_samples_JOB_SCRIPT.sh
cat << eof1 > ${script_name}
#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=0-$((${n_array_jobs[@]} - 1))
#SBATCH --mem=20gb

phens_to_run=(${phens_to_run[@]})
samples_to_run=(${samples_to_run[@]})
missing_array_ids=(${missing_array_ids[@]})
array_id=\${missing_array_ids[\${SLURM_ARRAY_TASK_ID}]}
phen=\${phens_to_run[\${array_id}]}
sample=\${samples_to_run[\${array_id}]}

# choose which cohorts are available for this pheno
module load lang/r/4.0.3-bioconductor-gcc
# save output of metal_ukb_alspac_hunt_choose_cohorts.R to bash array ${run_cohorts}; metal_ukb_alspac_hunt_choose_cohorts.R will require editing if passing more args
run_cohorts=(\$(Rscript ${script_path}metal_ukb_alspac_hunt_choose_cohorts.R \${phen} ${cohorts[@]}))
n_cohorts=\${#run_cohorts[@]}

# set input file names ending
if [ \${sample} == 1 ] || [ \${sample} == 2 ] || [ \${sample} == 3 ]; then file_end=chr; fi
if [ \${sample} == 4m ] || [ \${sample} == 4f ]; then file_end=snps; fi

# set variables for this run
dir_data=${qc_dir}${ma_name}/sample_\${sample}_\${phen}_files

# copy input files to temp dir for this run
mkdir -p \${dir_data}
cp ${gwas_path}output/${ma_name}/${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_1.tbl.gz \
\${dir_data}
for((i = 0; i < \${#run_cohorts[@]}; i++)); do
  cohort=\${run_cohorts[\${i}]}
  cp ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.gz \
  \${dir_data}
done

# write config.ini file
qc_config_file=\${dir_data}/config_file_\${phen}_sample\${sample}_all_\${file_end}.txt
cat ${script_path}gwas_inspector_config_file_part_1.txt > \${qc_config_file}
filename_output_tag=sample_\${sample}_\${phen}
echo -e 'filename_output_tag = '\${filename_output_tag^^} >> \${qc_config_file}
echo -e 'dir_data = '\${dir_data} >> \${qc_config_file}
echo -e 'dir_output = '\${dir_data} >> \${qc_config_file}
echo -e 'dir_references = '${dir_references} >> \${qc_config_file}
cat ${script_path}gwas_inspector_config_file_part_2.txt >> \${qc_config_file}
echo -e 'header_translations = '${header_translations} >> \${qc_config_file}
echo -e 'allele_ref_std = '${allele_ref_std} >> \${qc_config_file}
echo -e 'allele_ref_std_population = '${allele_ref_std_population} >> \${qc_config_file}
cat ${script_path}gwas_inspector_config_file_part_3.txt >> \${qc_config_file}
echo -e 'file_order_string = c("'${ma_name}'_sample_'\${sample}'_'\${phen}'_all_'\${file_end}'_1")' >> \${qc_config_file}
cat ${script_path}gwas_inspector_config_file_part_4.txt >> \${qc_config_file}

# run qc
Rscript ${script_path}run_gwasinspector.R \${qc_config_file}

# move reports to overall output dir
mv \${dir_data}/\${filename_output_tag^^}_Report.html ${dir_output}
mv \${dir_data}/\${filename_output_tag^^}_${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_1.tbl_report.html ${dir_output}
for((i = 0; i < \${#run_cohorts[@]}; i++)); do
  cohort=\${run_cohorts[\${i}]}
  mv \${dir_data}/\${filename_output_tag^^}_\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA_report.html \
  ${dir_output}
done

# delete redundant plots
rm \${dir_data}/*histogram.png
rm \${dir_data}/*sm.png
rm \${dir_data}/*correlation.png

# cleanup input files
rm \${dir_data}/${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_1.tbl.gz
for((i = 0; i < \${#run_cohorts[@]}; i++)); do
  cohort=\${run_cohorts[\${i}]}
  rm \${dir_data}/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.gz
done

eof1
chmod +x ${script_name}
sbatch ${script_name}

