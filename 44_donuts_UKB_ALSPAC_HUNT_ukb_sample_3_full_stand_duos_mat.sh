#!/bin/bash
#SBATCH --job-name=wlm
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --array=0-10
#SBATCH --mem=60gb
#SBATCH --account=SSCM013522

# run wlm; save conditional gwas sum stats
# bluepebble

ma_name=(UKB_ALSPAC_HUNT) # string to identify this MA
extra_text=_ukb_sample_3_full_stand # don't add _duos_mat here, but remember to add it to out_file in Rscript call below
gwas_path=foo
script_path=foo
plot_path=foo/${ma_name}/
alpha=0 # spousal genotypic correlation
threads=1
phens=(lnbmi lnwhr sbp dbp lnglu lnhba1c chol lnhdl ldl lntg lncrp) # n.b. there are 11 phenos, so array job should be 0-10
phen=${phens[${SLURM_ARRAY_TASK_ID}]}
in_file_stem1=${ma_name}_sample_1_${phen}_all_chr${extra_text}_1
in_file_stem2=${ma_name}_sample_2_${phen}_all_chr${extra_text}_1
in_file_stem3=${ma_name}_sample_3_${phen}_all_chr${extra_text}_1
ldsc_file_stem1=ldsc_${phen}_child_mother${extra_text}.log
ldsc_file_stem2=ldsc_${phen}_child_father${extra_text}.log
ldsc_file_stem3=ldsc_${phen}_mother_father${extra_text}.log
# hack to run lnbmi and lnglu with reduced nsnp ldsc output
if [[ ${phen} == lnbmi ]]; then ldsc_file_stem3=ldsc_lnbmi_mother_father${extra_text}_SUBSET_1000000_VARIANTS_SAMPLE_2.log; fi
if [[ ${phen} == lnglu ]]; then ldsc_file_stem1=ldsc_lnglu_child_mother${extra_text}_SUBSET_1000000_VARIANTS_SAMPLE_1.log; fi

module load lang/r/4.0.3-bioconductor-gcc # (current default)
Rscript ${script_path}run_donuts_duos_mat.R \
${gwas_path}output/${ma_name}/ \
${in_file_stem1}.tbl.gz \
${in_file_stem2}.tbl.gz \
${in_file_stem3}.tbl.gz \
${ldsc_file_stem1} \
${ldsc_file_stem2} \
${ldsc_file_stem3} \
${ma_name}_${phen}_donuts${extra_text}_duos_mat \
${threads} \
${phen} \
${plot_path} \
${alpha} \
${extra_text}_duos_mat # this last argument appears to be redundant now



