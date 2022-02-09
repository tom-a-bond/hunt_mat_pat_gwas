#!/bin/bash
# Tom Bond, 7th Feb 2022

# this script writes pheno, covar and begn files for unrelated individuals in samples 4m/4f

# 1: write ID lists for all mothers and fathers in trios sample (4m/4f) ####

# the code below was adapted from MP_BMI_MR_hunt_analysis_v2.R
# begin R code
library(data.table)
library(dplyr)
library(stringr)
rm(list = ls())
in_path = 'C:/Users/tb17613/OneDrive - The University of Queensland/work/proj/maternal_paternal_bmi_mr/hunt/' # path to data files contianing pehnotyep data for DMO, DFO, DTRIOS etc
out_path = in_path # path to write output files (ids lists, phenotype and covariate files)
hunt_phen = 'BMI'
hunt_phen_standardised = 'bmi'
# load analysis datasets ####
dtrios = read.table(paste0(in_path, 'Cleaned_', hunt_phen, '_TRIO.txt'), header = TRUE, stringsAsFactors = FALSE)
# get original ids for parents (those present in plink/bgen files, prior to adding suffixes)
dtrios$id_m_original = unlist(lapply(str_split(dtrios$id_m, fixed('_')), function(x){
  x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
}))
dtrios$id_f_original = unlist(lapply(str_split(dtrios$id_f, fixed('_')), function(x){
  x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
}))
# get suffixes added to parents ids
dtrios$id_m_suffix = unlist(lapply(str_split(dtrios$id_m, fixed('_')), function(x){
  x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
}))
dtrios$id_f_suffix = unlist(lapply(str_split(dtrios$id_f, fixed('_')), function(x){
  x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
}))
# write id lists for subsetting plink/bgen files
# mother-father-offspring trios (sample 4m/4f)
# mothers
write.table(dtrios[c('id_m_original', 'id_m_original')], paste0(out_path, 'dtrios_ids_mothers_all.txt'),
            qu = FALSE, ro = FALSE, co = FALSE)
write.table(dtrios[c('id_m_suffix', 'id_m_suffix')], paste0(out_path, 'dtrios_id_suffixes_mothers_all.txt'),
            qu = FALSE, ro = FALSE, co = FALSE)
# fathers
write.table(dtrios[c('id_f_original', 'id_f_original')], paste0(out_path, 'dtrios_ids_fathers_all.txt'),
            qu = FALSE, ro = FALSE, co = FALSE)
write.table(dtrios[c('id_f_suffix', 'id_f_suffix')], paste0(out_path, 'dtrios_id_suffixes_fathers_all.txt'),
            qu = FALSE, ro = FALSE, co = FALSE)
# end of R code

# 2: subset plink called genotype data to mothers, fathers and offspring from genotyped trios data set ####

# the code below is adapted from hunt_subset_plink_files.sh
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
plink2=/home/ttbond/plink2_linux_x86_64/plink2 # plink 2 executable. Must be plink 2, not plink 1.9
gen_path=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/CleanGenotypeData/ # path to plink 1 binary fileset containing qc'd called genotypes for whole hunt sample (I assume there is only one set of plink files for the entire cohort)
gen_file=ukb_cal_chr1_22_v2_cleaned # filename stem for plink 1 binary fileset containing qc'd called genotypes for whole hunt sample (I assume there is only one set of plink files for the entire cohort)
keep_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/alspac/ # path to id lists written above (i.e. "out_path" in the R code above)
out_path_gen=/path/to/write/plink/files/to/ # directory that subsetted plink files will be written to
out_file_sample4m=${cohort_name}_sample4m_maternal_geno # output file name for sample 4m (conditional trios gwas; maternal ids)
out_file_sample4f=${cohort_name}_sample4f_paternal_geno # output file name for sample 4f (conditional trios gwas; paternal ids)
out_file_sample4o=${cohort_name}_sample4_offspring_geno # output file name for sample 4f (conditional trios gwas; paternal ids)
# mothers
${plink2} \
--bfile ${gen_path}${gen_file} \
--keep ${keep_path}dtrios_ids_mothers_all.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample4m}_all
# fathers
${plink2} \
--bfile ${gen_path}${gen_file} \
--keep ${keep_path}dtrios_ids_fathers_all.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample4f}_all
# offspring genotype; id list was written already by MP_BMI_MR_hunt_analysis_v2.R
${plink2} \
--bfile ${gen_path}${gen_file} \
--keep ${keep_path}dtrios_ids_offspring.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample4o}_all

# 3: get lists of unrelated individuals ####

# using 20 threads
threads=20
# mothers
${plink2} \
--bfile ${out_path_gen}${out_file_sample4m}_all \
--autosome \
--threads ${threads} \
--king-cutoff 0.04419417 \
--out ${out_path_gen}${out_file_sample4m}_all
# fathers
${plink2} \
--bfile ${out_path_gen}${out_file_sample4f}_all \
--autosome \
--threads ${threads} \
--king-cutoff 0.04419417 \
--out ${out_path_gen}${out_file_sample4f}_all
# offspring
${plink2} \
--bfile ${out_path_gen}${out_file_sample4o}_all \
--autosome \
--threads ${threads} \
--king-cutoff 0.04419417 \
--out ${out_path_gen}${out_file_sample4o}_all

# 4: write phenotype and covariate files for unrelated samples (4m/4f) ####

# the R code below was adapted from MP_BMI_MR_hunt_analysis_v2.R
# begin R code
library(data.table)
library(dplyr)
library(stringr)
rm(list = ls())
# loop over phenotypes
in_path = 'C:/Users/tb17613/OneDrive - The University of Queensland/work/proj/maternal_paternal_bmi_mr/hunt/' # path to data files
out_path = in_path # path to write output files (ids lists, phenotype and covariate files)
hunt_phens = c('BMI', 'CHO', 'CRP', 'DBP', 'GLU', 'HBA', 'HDL', 'LDL', 'SBP', 'TG', 'WHR', 'BW', 'GA',
               'BMI_log', 'WHR_log', 'GLU_log', 'HBA_log', 'HDL_log', 'TG_log', 'CRP_log')
hunt_phens_standardised = c('bmi', 'chol', 'crp', 'dbp', 'glu', 'hba1c', 'hdl', 'ldl', 'sbp', 'tg', 'whr', 'bw', 'gest',
               'lnbmi', 'lnwhr', 'lnglu', 'lnhba1c', 'lnhdl', 'lntg', 'lncrp')
for(i in 1:length(hunt_phens)){
  hunt_phen = hunt_phens[i]
  hunt_phen_standardised = hunt_phens_standardised[i]
  writeLines(hunt_phen)
  # load analysis datasets ####
  dtrios = read.table(paste0(in_path, 'Cleaned_', hunt_phen, '_TRIO.txt'), header = TRUE, stringsAsFactors = FALSE)
  # get original ids for parents (those present in plink/bgen files, prior to adding suffixes)
  dtrios$id_m_original = unlist(lapply(str_split(dtrios$id_m, fixed('_')), function(x){
    x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
  }))
  dtrios$id_f_original = unlist(lapply(str_split(dtrios$id_f, fixed('_')), function(x){
    x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
  }))
  # get suffixes added to parents ids
  dtrios$id_m_suffix = unlist(lapply(str_split(dtrios$id_m, fixed('_')), function(x){
    x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
  }))
  dtrios$id_f_suffix = unlist(lapply(str_split(dtrios$id_f, fixed('_')), function(x){
    x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
  }))
  # remove related individuals
  unrel_mother = read.table('${out_path_gen}${out_file_sample4m}_all.king.cutoff.in.id', header = FALSE) # n.b. please check whether these files have a header
  unrel_father = read.table('${out_path_gen}${out_file_sample4f}_all.king.cutoff.in.id', header = FALSE)
  unrel_offs = read.table('${out_path_gen}${out_file_sample4o}_all.king.cutoff.in.id', header = FALSE)
  unrel_mother$unrel_mother = rep(1, nrow(unrel_mother))
  unrel_father$unrel_father = rep(1, nrow(unrel_father))
  unrel_offs$unrel_offs = rep(1, nrow(unrel_offs))
  writeLines(paste0('nrow dtrios before removing relateds:', nrow(dtrios))
  dtrios = merge(dtrios, unrel_mother[c('V1', 'rel_mother')], by.x = 'id_m_original', by.y = 'V1', all.x = TRUE)
  dtrios = dtrios[which(dtrios$unrel_mother == 1), ]
  writeLines(paste0('nrow dtrios after removing related mothers:', nrow(dtrios))
  dtrios = merge(dtrios, unrel_father[c('V1', 'rel_father')], by.x = 'id_f_original', by.y = 'V1', all.x = TRUE)
  dtrios = dtrios[which(dtrios$unrel_father == 1), ]
  writeLines(paste0('nrow dtrios after removing related fathers:', nrow(dtrios))
  dtrios = merge(dtrios, unrel_offs[c('V1', 'rel_offs')], by.x = 'id_o', by.y = 'V1', all.x = TRUE)
  dtrios = dtrios[which(dtrios$unrel_offs == 1), ]
  writeLines(paste0('nrow dtrios after removing related offspring:', nrow(dtrios))
  # write phenotype files
  # sample 4m (trios, mother ids)
  write.table(dtrios[c('id_m', 'id_m', 'PHE_o')],
              paste0(out_path, 'sample4m_', hunt_phen_standardised, '_unrelated.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 4f (trios, father ids)
  write.table(dtrios[c('id_f', 'id_f', 'PHE_o')],
              paste0(out_path, 'sample4f_', hunt_phen_standardised, '_unrelated.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)
  # covariate files
  if(hunt_phen == 'BMI'){ # the covariate files should be identical for each phenotype, so only need to be written once
    cont_covars = c('Sex', 'AGE', 'AGE2', 'AGESEX', 'AGE2SEX')
    pcs = sprintf('PC%s', 1:20)
    cat_covars = c('batch', 'PR')
    # n.b. for samples 4m and 4f the offspring's and other parent's genotype will be added as covariates later
    # sample 4m (trios, mother ids)
    write.table(dtrios[c('id_m_original', 'id_m_original', paste0(cont_covars, '_o'), paste0(pcs, '_m'), 'batch_m', 'PR_o')], # n.b. using offspring covariates, parental pcs, parental genotyping batch and offspring hunt wave
                paste0(out_path, 'sample4m_', hunt_phen_standardised, '_unrelated.covar'), # all covariates (continuous and categorical)
                qu = FALSE, ro = FALSE, co = FALSE)
    # sample 4f (trios, father ids)
    write.table(dtrios[c('id_f_original', 'id_f_original', paste0(cont_covars, '_o'), paste0(pcs, '_f'), 'batch_f', 'PR_o')], # n.b. using offspring covariates, parental pcs, parental genotyping batch and offspring hunt wave
                paste0(out_path, 'sample4f_', hunt_phen_standardised, '_unrelated.covar'), # all covariates (continuous and categorical)
                qu = FALSE, ro = FALSE, co = FALSE)
    # linker file
    write.table(dtrios[c('id_m_original', 'id_f_original', 'id_o')], paste0(out_path, 'dtrios_id_linker_original_ids.txt'),
                qu = FALSE, ro = FALSE, co = TRUE)
  }
}
# end of R code

# 5: subset bgens to 495 bmi snps, for mothers, fathers and offspring samples
# n.b. we are not excluding relatedness from these bgen files, but instead removing related individuals from the phenotype and covariate files

# code adapted from subset_hunt_bgen_495_yengo_bmi_snps.sh
bgenpath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Data/Genetic_v3/ # path to bgen files
samplepath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/GeneticLinkingFiles/ # path to bgen .sample file
outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # folder to write output (subsetted bgen files) to
inpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/ # folder containing yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt
# keep_path is defined above
qctool=/home/ttbond/qctool/build/release/qctool_v2.0.8 # path to qc tool executable
bgenix=/path/to/bgenix/executable # path to bgenix executable
n_yengo_snps=495 # 495 yengo snps should be available in all cohorts
# subset the full hunt sample bgen files (495 snps) to the samples required for sample4m and sample4f fastgwa analysis
# the input bgen files for the qctool calls below were written by subset_hunt_bgen_495_yengo_bmi_snps.sh
# mothers
${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-incl-samples ${keep_path}dtrios_ids_mothers_all.txt \
-bgen-omit-sample-identifier-block \
-og "${outpath}"hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.bgen \
-os "${outpath}"hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.sample
# fathers  
${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-incl-samples ${keep_path}dtrios_ids_fathers_all.txt \
-bgen-omit-sample-identifier-block \
-og "${outpath}"hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.bgen \
-os "${outpath}"hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.sample
# offspring  
${qctool} \
-g "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr#.bgen \
-s "${outpath}"hunt_full_sample_yengo_${n_yengo_snps}_bmi_snps_chr1.sample \
-incl-samples ${keep_path}dtrios_ids_offspring.txt \
-og "${outpath}"hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.bgen \
-os "${outpath}"hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.sample
# create bgi indexes
${bgenix} -g "${outpath}"hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.bgen -index
${bgenix} -g "${outpath}"hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.bgen -index
${bgenix} -g "${outpath}"hunt_sample4_offspring_yengo_${n_yengo_snps}_bmi_snps_all_chr_unrelated.bgen -index

