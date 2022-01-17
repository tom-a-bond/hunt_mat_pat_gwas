
# fastgwa gwas, sample 4m (trios, maternal genotype conditional analyses)
# merge output from individual snps to one file, for each phenotype

library(dplyr)
library(data.table)

rm(list = ls())
cohort_name = 'hunt' # cohort name; one of "ukb", "hunt", "alspac"
sample = '4m' # 4m or 4f
snp_list_path = '/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/' # path to list of Yengo et al BMI snps (yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt)
out_path = '/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/results/gwas/duos_trios/' # output path (directory containing folders written by fastgwa_hunt_sample4m.sh)
snp_list = 'yengo_bmi_snps_manual_version_ukb_alspac_hunt_maf_0_01_r2_0_8_rsids.txt'
outcomes = c('bmi', 'chol', 'crp', 'dbp', 'glu', 'hba1c', 'hdl', 'ldl', 'sbp', 'tg', 'whr', 'bw', 'gest',
             'lnbmi', 'lnwhr', 'lnglu', 'lnhba1c', 'lnhdl', 'lntg', 'lncrp')
snps = scan(paste0(snp_list_path, snp_list), what = 'character')
n_snps = length(snps)

# check for missing output files
missings=character()
for(j in 1:length(outcomes)){
  outcome = outcomes[j]
  writeLines(outcome)
  run = paste0(cohort_name, '_', sample, '_', outcome)
  setwd(paste0(out_path, run))
  file_list = list.files()
  out_files = file_list[grep('.fastGWA', file_list)]
  if(length(out_files) != n_snps){ missings = c(missings, outcome) }
}
missings
length(missings) # should be 0

# merge fastGWA results for all snps, write output to parent directory
for(j in 1:length(outcomes)){
  outcome = outcomes[j]
  writeLines(outcome)
  run = paste0(cohort_name, '_', sample, '_', outcome)
  setwd(paste0(out_path, run))
  file_list = list.files()
  out_files = file_list[grep('.fastGWA', file_list)]
  if(length(out_files) != n_snps){ stop('Results for one or more snps appear to be missing') }
  res = lapply(as.list(out_files), fread, header = TRUE, data.table = FALSE) %>% do.call(rbind, .)
  fwrite(res, paste0('../', cohort_name, '_sample_', sample, '_', outcome, '_all_snps.fastGWA.gz'),
         row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ' ', compress = 'gzip') # writes output file to one directory up
}
  
  