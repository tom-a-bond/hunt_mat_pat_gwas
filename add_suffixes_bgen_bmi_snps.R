args = commandArgs(trailingOnly = TRUE)
keep_path = args[1]
outpath = args[2]
n_yengo_snps = args[3]
library(data.table)
for(i in 1:6){
  # mothers
  sample_m = fread(paste0(outpath, 'hunt_sample4_mothers_', i, '_yengo_', n_yengo_snps, '_bmi_snps_all_chr.sample'),
                   data.table = FALSE, header = FALSE, skip = 2)
  suffix_m = fread(paste0(keep_path, 'dtrios_id_suffixes_mothers_', i, '.txt'),
                   data.table = FALSE, header = FALSE)
  sample_m[ , 1] = paste0(sample_m[ , 1], suffix_m[ , 1])
  sample_m[ , 2] = paste0(sample_m[ , 2], suffix_m[ , 2])
  fwrite(rbind(c('ID_1', 'ID_2', 'missing', 'sex'),
               c(0, 0, 0, 'D'),
               sample_m), paste0(outpath, 'hunt_sample4_mothers_', i, '_yengo_', n_yengo_snps, '_bmi_snps_all_chr_suffixes.sample'),
         sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  # fathers
  sample_f = fread(paste0(outpath, 'hunt_sample4_fathers_', i, '_yengo_', n_yengo_snps, '_bmi_snps_all_chr.sample'),
                   data.table = FALSE, header = FALSE, skip = 2)
  suffix_f = fread(paste0(keep_path, 'dtrios_id_suffixes_fathers_', i, '.txt'),
                   data.table = FALSE, header = FALSE)
  sample_f[ , 1] = paste0(sample_f[ , 1], suffix_f[ , 1])
  sample_f[ , 2] = paste0(sample_f[ , 2], suffix_f[ , 2])
  fwrite(rbind(c('ID_1', 'ID_2', 'missing', 'sex'),
               c(0, 0, 0, 'D'),
               sample_f), paste0(outpath, 'hunt_sample4_fathers_', i, '_yengo_', n_yengo_snps, '_bmi_snps_all_chr_suffixes.sample'),
         sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
}
