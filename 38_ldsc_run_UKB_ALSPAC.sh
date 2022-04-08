#!/bin/bash
#SBATCH --job-name=ldsc
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --array=0-22
#SBATCH --mem=10gb

# run bivariate ldsc for mother-child, father-child and mother-father gwas sum stats

ma_name=(UKB_ALSPAC) # string to identify this MA
gwas_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/
phens=(`cat ${gwas_path}input/ph_ukb_alspac.txt`) # n.b. there are 23 phenos, so array job should be 0-22
script_path=/user/work/tb17613/proj/mat_pat_bmi_mr/scripts/
ldsc_path=/user/home/tb17613/software/ldsc/
ldscores_path=/user/home/tb17613/software/ldsc/downloaded_ld_scores/eur_w_ld_chr/
phen=${phens[${SLURM_ARRAY_TASK_ID}]}

module load lang/python/anaconda/2.7-2019.10
source activate ldsc

# child-mother

${ldsc_path}ldsc.py \
--rg \
${gwas_path}output/${ma_name}/${ma_name}_sample_1_${phen}_all_chr_1_ldhub_qc.sumstats.gz,${gwas_path}output/${ma_name}/${ma_name}_sample_3_${phen}_all_chr_1_ldhub_qc.sumstats.gz \
--ref-ld-chr ${ldscores_path} \
--w-ld-chr ${ldscores_path} \
--out ${gwas_path}output/${ma_name}/ldsc_${phen}_child_mother

# child-father

${ldsc_path}ldsc.py \
--rg \
${gwas_path}output/${ma_name}/${ma_name}_sample_2_${phen}_all_chr_1_ldhub_qc.sumstats.gz,${gwas_path}output/${ma_name}/${ma_name}_sample_3_${phen}_all_chr_1_ldhub_qc.sumstats.gz \
--ref-ld-chr ${ldscores_path} \
--w-ld-chr ${ldscores_path} \
--out ${gwas_path}output/${ma_name}/ldsc_${phen}_child_father

# mother-father

${ldsc_path}ldsc.py \
--rg \
${gwas_path}output/${ma_name}/${ma_name}_sample_1_${phen}_all_chr_1_ldhub_qc.sumstats.gz,${gwas_path}output/${ma_name}/${ma_name}_sample_2_${phen}_all_chr_1_ldhub_qc.sumstats.gz \
--ref-ld-chr ${ldscores_path} \
--w-ld-chr ${ldscores_path} \
--out ${gwas_path}output/${ma_name}/ldsc_${phen}_mother_father


