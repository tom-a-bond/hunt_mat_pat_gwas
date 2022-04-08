#!/bin/bash

# run on bluepebble

ma_name=(UKB_ALSPAC) # string to identify this MA
cohorts=(ukb alspac)
samples=(1 2 3 4m 4f)
#samples=(4m 4f)
gwas_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/
phens=(`cat ${gwas_path}input/ph_ukb_alspac.txt`)
script_path=/user/work/tb17613/proj/mat_pat_bmi_mr/scripts/
metal=/user/home/tb17613/software/generic-metal/metal
drop_ambig=false # logical; whether to drop ambiguous snps (A/T, C/G)
chr_pos_a1_a2=false # logical; whether to make chr:pos_a1/a2 column (currently the input data already have this, so set to false)
if [[ (${sample} == 4m || ${sample} == 4f) && (${chr_pos_a1_a2} == true || ${drop_ambig} == true) ]]; then
  echo 'Error: drop_ambig and chr_pos_a1_a2 should both be false for sample 4m/4f'; exit 1
fi
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
n_array_jobs=$((${#samples_to_run[@]}))
# make output dir
out_path=${gwas_path}output/${ma_name}/
mkdir -p ${out_path}
mkdir -p ${out_path}metal_scripts

script_name=${script_path}metal_${ma_name}_all_samples_JOB_SCRIPT.sh
cat << eof1 > ${script_name}
#!/bin/bash
#SBATCH --job-name=ma
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --array=0-$((${n_array_jobs[@]} - 1))
#SBATCH --mem=10gb

phens_to_run=(${phens_to_run[@]})
samples_to_run=(${samples_to_run[@]})
phen=\${phens_to_run[\${SLURM_ARRAY_TASK_ID}]}
sample=\${samples_to_run[\${SLURM_ARRAY_TASK_ID}]}

# choose which cohorts are available for this pheno
module load lang/r/4.0.3-bioconductor-gcc
# save output of metal_ukb_alspac_hunt_choose_cohorts.R to bash array ${run_cohorts}; metal_ukb_alspac_hunt_choose_cohorts.R will require editing if passing more args
run_cohorts=(\$(Rscript ${script_path}metal_ukb_alspac_hunt_choose_cohorts.R \${phen} ${cohorts[@]}))
n_cohorts=\${#run_cohorts[@]}

# set input file names ending
if [ \${sample} == 1 ] || [ \${sample} == 2 ] || [ \${sample} == 3 ]; then file_end=chr; fi
if [ \${sample} == 4m ] || [ \${sample} == 4f ]; then file_end=snps; fi

# write metal script
metal_script=${out_path}metal_scripts/metal_script_\${phen}_sample\${sample}_all_\${file_end}.txt
# write set options section (for sample 4 we will not use the AVERAGEFREQ and MINMAXFREQ options because plink does not provide allele freq in gwas output)
if [ \${sample} == 1 ] || [ \${sample} == 2 ] || [ \${sample} == 3 ]; then
  echo -e '# set options\nSCHEME STDERR\nAVERAGEFREQ ON\nMINMAXFREQ ON\nCUSTOMVARIABLE total_n' > \${metal_script}
fi
if [ \${sample} == 4m ] || [ \${sample} == 4f ]; then
  echo -e '# set options\nSCHEME STDERR\nCUSTOMVARIABLE total_n' > \${metal_script}
fi
for((i = 0; i < \${#run_cohorts[@]}; i++)); do
  cohort=\${run_cohorts[\${i}]}
  # write describe command
  if [ \${sample} == 1 ] || [ \${sample} == 2 ] || [ \${sample} == 3 ]; then
    echo -e '# describe and process the input file for ith cohort\nMARKER CHR:POS_A1/A2\nALLELE A1 A2\nFREQ AF1\nEFFECT BETA\nSTDERR SE\nPVAL P\nLABEL total_n as N' >> \${metal_script}
  fi
  if [ \${sample} == 4m ] || [ \${sample} == 4f ]; then
    echo -e '# describe and process the input file for ith cohort\nMARKER CHR:POS_A1/A2\nALLELE A1 A2\nEFFECT BETA\nSTDERR SE\nPVAL P\nLABEL total_n as N' >> \${metal_script}
  fi
  # write process command
  echo -e 'PROCESS '${gwas_path}'input/'\${cohort}'_to_send/'\${cohort}'_sample_'\${sample}'_'\${phen}'_all_'\${file_end}'.fastGWA' >> \${metal_script}
done

# write analyze command
echo -e 'OUTFILE '${out_path}${ma_name}'_sample_'\${sample}'_'\${phen}'_all_'\${file_end}'_ .tbl\nANALYZE HETEROGENEITY\nQUIT' >> \${metal_script}

# unzip input files, remove ambiguous snps and make chr_pos_alleles variable (not all cohorts have rsid)
for((i = 0; i < \${#run_cohorts[@]}; i++)); do
  cohort=\${run_cohorts[\${i}]}
  if [[ ${drop_ambig} == true ]]; then
    gzip -cd ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.gz > \
    ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.palind
    echo N snps for sample \${sample}, all chrs, \${phen}, \${cohort} before removing palindromic snps: \$((\$(wc -l < ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.palind) - 1))
    awk '!((\$4 == "A" && \$5 == "T") || (\$4 == "T" && \$5 == "A") || (\$4 == "G" && \$5 == "C") || (\$4 == "C" && \$5 == "G")) {print}' \
    ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.palind > \
    ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp
    echo N snps for sample \${sample}, all chrs, \${phen}, \${cohort} after removing palindromic snps: \$((\$(wc -l < ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp) - 1))
  fi
  if [[ ${drop_ambig} == false ]]; then
    gzip -cd ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.gz > \
    ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp
    # for alspac sample 2, drop problematic snps (see gwas_qc_test_bmi_alspac_only_DESKTOP.R)
    if [[ \${cohort} == alspac && \${sample} == 2 ]]; then
      Rscript ${script_path}alspac_drop_problem_snps.R \
      ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp \
      /user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/input/alspac_1000G_snps_for_ma_without_problem_palind.txt.gz
    fi
    echo N snps for sample \${sample}, all chrs, \${phen}, \${cohort} '(not removing palindromic snps)': \$((\$(wc -l < ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp) - 1))
  fi
  # make chr_pos_alleles column if required
  if [[ ${chr_pos_a1_a2} == true ]]; then
    awk '{print \$1":"\$3"_"\$4"/"\$5, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp > \
    ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA
  fi
  if [[ ${chr_pos_a1_a2} == false ]]; then
    cp ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp \
    ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA
  fi
done

# remove any output from identically named previous runs (because metal does not overwrite files [it suffixes numbers to filenames instead])
rm -f ${out_path}${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_*.tbl
rm -f ${out_path}${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_*.tbl.gz
rm -f ${out_path}${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_*.tbl.info

# run metal
chmod +x \${metal_script}
${metal} \${metal_script} > ${out_path}${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}.log

# delete unzipped input files
for((i = 0; i < \${#run_cohorts[@]}; i++)); do
  cohort=\${run_cohorts[\${i}]}
  rm -f ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.palind
  rm -f ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA.tmp
  rm -f ${gwas_path}input/\${cohort}_to_send/\${cohort}_sample_\${sample}_\${phen}_all_\${file_end}.fastGWA
done

# zip ma output files
gzip ${out_path}${ma_name}_sample_\${sample}_\${phen}_all_\${file_end}_*.tbl

eof1
chmod +x ${script_name}
sbatch ${script_name}

