#!/bin/bash
#SBATCH --job-name=wlm
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --array=0-22
#SBATCH --mem=50gb

# run wlm; save conditional gwas sum stats

ma_name=(UKB_ALSPAC) # string to identify this MA
gwas_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/
script_path=/user/work/tb17613/proj/mat_pat_bmi_mr/scripts/
plot_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/qc/${ma_name}/
threads=1
phens=(`cat ${gwas_path}input/ph_ukb_alspac.txt`) # n.b. there are 23 phenos, so array job should be 0-22
phen=${phens[${SLURM_ARRAY_TASK_ID}]}

module load lang/r/4.0.3-bioconductor-gcc # (current default)
Rscript ${script_path}run_wlm_trios.R \
${gwas_path}output/${ma_name}/ \
${ma_name}_sample_1_${phen}_all_chr_1.tbl.gz \
${ma_name}_sample_2_${phen}_all_chr_1.tbl.gz \
${ma_name}_sample_3_${phen}_all_chr_1.tbl.gz \
ldsc_${phen}_child_mother.log \
ldsc_${phen}_child_father.log \
ldsc_${phen}_mother_father.log \
${ma_name}_${phen}_wlm.tsv.gz \
${threads} \
${phen} \
${plot_path}


