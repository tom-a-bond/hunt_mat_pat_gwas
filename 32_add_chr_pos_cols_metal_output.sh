#!/bin/bash
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --array=0-22
#SBATCH --mem=5gb

# add chr and pos variables to metal output
# bluepebble

ma_name=UKB_ALSPAC
phens=(whr tg sbp lnwhr lntg lnhdl lnglu lnfmi lncrp lnbmi lnapob lnapoa ldl hdl glu fmi dbp crp chol bw bmi apob apoa)
phen=${phens[${SLURM_ARRAY_TASK_ID}]}
in_path=/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/output/

module load lang/r/4.0.3-bioconductor-gcc # (current default)
# sample 1
echo sample 1
Rscript /user/work/tb17613/proj/mat_pat_bmi_mr/scripts/add_chr_pos_cols_metal_output.R \
1 \
${phen} \
chr \
${in_path} \
${ma_name}
# sample 2
echo sample 2
Rscript /user/work/tb17613/proj/mat_pat_bmi_mr/scripts/add_chr_pos_cols_metal_output.R \
2 \
${phen} \
chr \
${in_path} \
${ma_name}
# sample 3
echo sample 3
Rscript /user/work/tb17613/proj/mat_pat_bmi_mr/scripts/add_chr_pos_cols_metal_output.R \
3 \
${phen} \
chr \
${in_path} \
${ma_name}
# sample 4m
echo sample 4m
Rscript /user/work/tb17613/proj/mat_pat_bmi_mr/scripts/add_chr_pos_cols_metal_output.R \
4m \
${phen} \
snps \
${in_path} \
${ma_name}
# sample 4f
echo sample 4f
Rscript /user/work/tb17613/proj/mat_pat_bmi_mr/scripts/add_chr_pos_cols_metal_output.R \
4f \
${phen} \
snps \
${in_path} \
${ma_name}



