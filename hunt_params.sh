#!/bin/bash

# variables used for running gwas models

################################################################################
##### variables which require editing (variables will be the same for all fastgwa scripts, i.e. sample 1, sample 2 etc):
################################################################################
relateds=true # true/false; will crytically related individuals be inluded in the sample (in which case we will use fastGWA to control for this, otherwise we will use linear regression instead. Should be true for HUNT)
gcta=/home/ttbond/gcta_1.93.3beta/gcta_1.93.3_20210817_fixed # gcta executable; should be v1.93.3beta2 or greater
qctool=/home/ttbond/qctool/build/release/qctool_v2.0.8 # qctool exectuable
plink2=/home/ttbond/plink2_linux_x86_64/plink2 # plink 2 executable. Must be plink 2, not plink 1.9
load_r_module() { # function to load R module
  echo loading R module; export MODULEPATH=/sw/el7/modulefiles; module load R/4.0.5 # commands required to get Rscript.exe "in path"
}
keep_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/alspac/ # path to id lists written by MP_BMI_MR_hunt_analysis.R (i.e. "out_path" in MP_BMI_MR_hunt_analysis.R)
keep_file_sample1=dmo_ids_mothers_reorder.txt # id list for sample 1 (marginal maternal gwas), written by MP_BMI_MR_hunt_analysis.R
keep_file_sample2=dfo_ids_fathers_reorder.txt # id list for sample 2 (marginal paternal gwas), written by MP_BMI_MR_hunt_analysis.R
keep_file_sample3=dall_ids.txt # id list for sample 3 (marginal offspring gwas), written by MP_BMI_MR_hunt_analysis.R
keep_file_sample4m=dtrios_ids_mothers_reorder.txt # id list for sample 4m (conditional trios gwas; maternal ids), written by MP_BMI_MR_hunt_analysis.R
keep_file_sample4f=dtrios_ids_fathers_reorder.txt # id list for sample 4f (conditional trios gwas; paternal ids), written by MP_BMI_MR_hunt_analysis.R
keep_file_sample4o=dtrios_ids_offspring.txt # id list for sample 4m/4f, offspring ids, written by MP_BMI_MR_hunt_analysis.R
linkage_file_sample4=dtrios_id_linker.txt # file to link mother, father and offspring ids for sample 4; should be in same directory as keep_file_sample4m; should have ids: "id_mo", "id_fa" and "id_of" for mother, father and offspring id respectively; written by MP_BMI_MR_hunt_analysis.R
grm_file_sample1=${cohort_name}_sample1_maternal_geno # grm file name for sample 1 (marginal maternal gwas)
grm_file_sample2=${cohort_name}_sample2_paternal_geno # grm file name for sample 2 (marginal paternal gwas)
grm_file_sample3=${cohort_name}_sample3_offspring_geno # grm file name for sample 3 (marginal offspring gwas)
grm_file_sample4m=${cohort_name}_sample4m_maternal_geno # grm file name for sample 4m (conditional trios gwas; maternal ids)
grm_file_sample4f=${cohort_name}_sample4f_paternal_geno # grm file name for sample 4f (conditional trios gwas; paternal ids)
grm_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # path to folder containing grms
sample1_bgen_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # path to subsetted bgen files for sample 1 ("outpath" in subset_hunt_bgen_sample1.sh)
sample2_bgen_path=${sample1_bgen_path} # path to subsetted bgen files for sample 2 ("outpath" in subset_hunt_bgen_sample2.sh)
sample3_bgen_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # path to bgen files for sample 3 (= original bgen files for full hunt smaple, i.e. not subsetted)
sample4_bgen_path=${sample1_bgen_path} # path to subsetted bgen files containing just Yengo et al BMI snps for sample 4m/4f ("outpath" in subset_hunt_bgen_495_yengo_bmi_snps.sh)
imp_path=${sample1_bgen_path} # path which hunt_bgen_file_lists.sh writes lists of bgens to
sample1_bgen_file=hunt_sample1_mothers_chr${PBS_ARRAY_INDEX}.bgen # subsetted bgen files for sample 1, written by subset_hunt_bgen_sample1.sh, chromosomes are indexed by ${PBS_ARRAY_INDEX}
sample2_bgen_file=hunt_sample2_fathers_chr${PBS_ARRAY_INDEX}.bgen # subsetted bgen files for sample 2, written by subset_hunt_bgen_sample2.sh, chromosomes are indexed by ${PBS_ARRAY_INDEX}
sample3_bgen_file=chr${PBS_ARRAY_INDEX}.bgen # bgen files for sample 3 (= original bgen files for full hunt smaple, i.e. not subsetted)
n_yengo_snps=495 # number of yengo snps available with high imputaiton quality in hunt; should be 495
sample4m_bgen_file=hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen # bgen file containing just Yengo et al BMI snps, for mothers in sample 4m (one file containing all 22 autosomes); written by subset_hunt_bgen_495_yengo_bmi_snps.sh
sample4f_bgen_file=hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen # bgen file containing just Yengo et al BMI snps, for fathers in sample 4f (one file containing all 22 autosomes); written by subset_hunt_bgen_495_yengo_bmi_snps.sh
sample4o_bgen_file=hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen # bgen file containing just Yengo et al BMI snps, for offspring in sample 4m/4f (one file containing all 22 autosomes); written by subset_hunt_bgen_495_yengo_bmi_snps.sh
sample1_sample_path=${sample1_bgen_path} # path to subsetted bgen .sample files for sample 1 ("outpath" in subset_hunt_bgen_sample1.sh)
sample2_sample_path=${sample2_bgen_path} # path to subsetted bgen .sample files for sample 2 ("outpath" in subset_hunt_bgen_sample2.sh)
sample3_sample_path=${sample3_bgen_path} # path to bgen .sample files for sample 3 (= original bgen files for full hunt smaple, i.e. not subsetted)
sample4_sample_path=${sample4_bgen_path} # path to subsetted bgen .sample files containing just Yengo et al BMI snps for sample 4m/4f ("outpath" in subset_hunt_bgen_495_yengo_bmi_snps.sh)
sample1_sample_file=hunt_sample1_mothers_chr${PBS_ARRAY_INDEX}.sample # subsetted bgen files for sample 1, written by subset_hunt_bgen_sample1.sh, chromosomes are indexed by ${PBS_ARRAY_INDEX}
sample2_sample_file=hunt_sample2_mothers_chr${PBS_ARRAY_INDEX}.sample # subsetted bgen files for sample 2, written by subset_hunt_bgen_sample2.sh, chromosomes are indexed by ${PBS_ARRAY_INDEX}
sample3_sample_file=chr${PBS_ARRAY_INDEX}.sample # bgen files for sample 3 (= original bgen files for full hunt smaple, i.e. not subsetted)
sample4m_sample_file=hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample # bgen file containing just Yengo et al BMI snps, for mothers in sample 4m (one file containing all 22 autosomes); written by subset_hunt_bgen_495_yengo_bmi_snps.sh
sample4f_sample_file=hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample # bgen file containing just Yengo et al BMI snps, for fathers in sample 4f (one file containing all 22 autosomes); written by subset_hunt_bgen_495_yengo_bmi_snps.sh
sample4o_sample_file=hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample # bgen file containing just Yengo et al BMI snps, for offspring in sample 4m/4f (one file containing all 22 autosomes); written by subset_hunt_bgen_495_yengo_bmi_snps.sh
out_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/results/gwas/duos_trios/ # output path for gwas results
phen_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/phen/ # path to phenotype files
covars_path=${phen_path} # path to covariate files
snp_list_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ # path directory containing list of yengo bmi snps (yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt)

################################################################################
##### variables which do not require editing (with possible exception of threads):
################################################################################
outcomes=(bmi chol crp dbp glu hba1c hdl ldl sbp tg whr bw gest lnbmi lnwhr lnglu lnhba1c lnhdl lntg lncrp)
fasting_outcomes=()
snp_list=yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt
if [[ ${sample} == '1' ]]; then keep_file_sample=${keep_file_sample1}; grm_file_sample=${grm_file_sample1}; sample_path=${sample1_sample_path}; sample_file=${sample1_sample_file}; fi
if [[ ${sample} == '2' ]]; then keep_file_sample=${keep_file_sample2}; grm_file_sample=${grm_file_sample2}; sample_path=${sample2_sample_path}; sample_file=${sample2_sample_file}; fi
if [[ ${sample} == '3' ]]; then keep_file_sample=${keep_file_sample3}; grm_file_sample=${grm_file_sample3}; sample_path=${sample3_sample_path}; sample_file=${sample3_sample_file}; fi
if [[ ${sample} == '4m' ]]; then keep_file_sample=${keep_file_sample4m}; grm_file_sample=${grm_file_sample4m}; sample_path=${sample4m_sample_path}; sample_file=${sample4m_sample_file}; sample_bgen_file=${sample4m_bgen_file}; fi
if [[ ${sample} == '4f' ]]; then keep_file_sample=${keep_file_sample4f}; grm_file_sample=${grm_file_sample4f}; sample_path=${sample4f_sample_path}; sample_file=${sample4f_sample_file}; sample_bgen_file=${sample4f_bgen_file}; fi
threads=20 # number of cores available
maf=0.01
info=0.3
geno=0.1 # fastgwa default; n.b. when using bgen format imputed data, this filter will not remove any snps

