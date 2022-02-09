#!/bin/bash

# extract 495 BMI snps from imputed genotype data stored in BGEN format, and calculate snp summary stats
# in UKB this script only takes a few minutes to run

bgenpath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Data/Genetic_v3/ # path to bgen files
samplepath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/GeneticLinkingFiles/ # path to bgen .sample file
outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # folder to write output (subsetted bgen files) to
inpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/ # folder containing yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt
keep_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/alspac/ # path to id lists written by MP_BMI_MR_hunt_analysis.R (i.e. "out_path" in MP_BMI_MR_hunt_analysis.R)
qctool=/home/ttbond/qctool/build/release/qctool_v2.0.8 # path to qc tool executable
bgenix=/path/to/bgenix/executable # path to bgenix executable
n_yengo_snps=495 # 495 yengo snps should be available in all cohorts
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory containing add_suffixes_bgen_bmi_snps.R

# first extract 495 yengo et al. bmi snps, writing one output file per chromosome for the full hunt sample

for i in {1..22}; do
  echo chr ${i}
  ${bgenix} \
  -g ${bgenpath}chr"${i}".bgen \
  -incl-rsids ${inpath}data/yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt > \
  "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr${i}.bgen
done
  
# then subset the full hunt sample bgen files (495 snps) to the samples required for sample4m and sample4f fastgwa analysis
# j indexes the id lists for the 1st ... 6th siblings, written by MP_BMI_MR_hunt_analysis_v2.R

for j in {1..6}; do
# mothers
${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-incl-samples ${keep_path}dtrios_ids_mothers_${j}.txt \
-bgen-omit-sample-identifier-block \
-og "${outpath}"hunt_sample4_mothers_${j}_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-os "${outpath}"hunt_sample4_mothers_${j}_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample
# fathers  
${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-incl-samples ${keep_path}dtrios_ids_fathers_${j}.txt \
-bgen-omit-sample-identifier-block \
-og "${outpath}"hunt_sample4_fathers_${j}_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-os "${outpath}"hunt_sample4_fathers_${j}_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample
done
# offspring  
${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-incl-samples ${keep_path}dtrios_ids_offspring.txt \
-og "${outpath}"hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-os "${outpath}"hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample

# PLEASE CHECK THAT THE SAMPLE FILES WRITTEN ABOVE HAVE THE CORRECT NUMBER OF IDS (I.E. THE SAME AS IN THE ID LISTS [dtrios_ids_mothers_${j}.txt etc.])

# add suffixes to parental ids

export MODULEPATH=/sw/el7/modulefiles; module load R/4.0.5 # put any commands required to get Rscript.exe "in path" here
Rscript ${script_path}add_suffixes_bgen_bmi_snps.R \
${keep_path} \
${outpath} \
${n_yengo_snps}

# merge parental files for sibs 1-6

# mothers
${qctool} \
-g "${outpath}"hunt_sample4_mothers_1_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_mothers_1_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_mothers_2_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_mothers_2_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_mothers_3_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_mothers_3_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_mothers_4_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_mothers_4_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_mothers_5_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_mothers_5_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_mothers_6_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_mothers_6_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-og "${outpath}"hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-os "${outpath}"hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample
# fathers
${qctool} \
-g "${outpath}"hunt_sample4_fathers_1_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_fathers_1_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_fathers_2_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_fathers_2_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_fathers_3_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_fathers_3_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_fathers_4_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_fathers_4_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_fathers_5_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_fathers_5_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-g "${outpath}"hunt_sample4_fathers_6_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-s "${outpath}"hunt_sample4_fathers_6_yengo_${n_yengo_snps}_bmi_snps_all_chr_suffixes.sample \
-og "${outpath}"hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen \
-os "${outpath}"hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr.sample

# PLEASE CHECK THAT THE SAMPLE FILES PRODUCED IMMEDIATELY ABOVE HAVE THE CORRECT NUMBER OF IDS
# (I.E. IT SHOULD BE THE SAME NUMBER OF IDS AS IN dtrios)

# calculate stats including info scores and maf (full hunt sample)

${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-snp-stats \
-osnp "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_all_chr.sumstats




