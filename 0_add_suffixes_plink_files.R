args = commandArgs(trailingOnly = TRUE)
keep_path = args[1]
outpath = args[2]
out_file_sample1 = args[3]
out_file_sample2 = args[4]
# out_file_sample4m = args[5]
# out_file_sample4f = args[6]
library(data.table)
for(i in 1:6){
  # sample 1
  fam_1 = fread(paste0(outpath, out_file_sample1, '_', i, '.fam'),
                data.table = FALSE, header = FALSE, skip = 2)
  suffix_1 = fread(paste0(keep_path, 'dmo_id_suffixes_mothers_', i, '.txt'),
                   data.table = FALSE, header = TRUE)
  if(nrow(fam_1) != nrow(suffix_1)){ Stop(paste0('fam file and suffix file for sample 1 have different number of rows (i = ', i, ')')) }
  fam_1[ , 1] = paste0(fam_1[ , 1], suffix_1[ , 1])
  fam_1[ , 2] = paste0(fam_1[ , 2], suffix_1[ , 2])
  fwrite(fam_1, paste0(outpath, out_file_sample1, '_', i, '_suffixes.fam'),
         sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  # sample 2
  fam_2 = fread(paste0(outpath, out_file_sample2, '_', i, '.fam'),
                data.table = FALSE, header = FALSE, skip = 2)
  suffix_2 = fread(paste0(keep_path, 'dfo_id_suffixes_fathers_', i, '.txt'),
                   data.table = FALSE, header = TRUE)
  if(nrow(fam_2) != nrow(suffix_2)){ Stop(paste0('fam file and suffix file for sample 2 have different number of rows (i = ', i, ')')) }
  fam_2[ , 1] = paste0(fam_2[ , 1], suffix_2[ , 1])
  fam_2[ , 2] = paste0(fam_2[ , 2], suffix_2[ , 2])
  fwrite(fam_2, paste0(outpath, out_file_sample2, '_', i, '_suffixes.fam'),
         sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  # sample 3 does not need to be subsetted
  # sample 4f/4f no longer needs to be subsetted, because we will use unrelateds only
  # # sample 4m
  # fam_4m = fread(paste0(outpath, out_file_sample4m, '_', i, '.fam'),
  #                data.table = FALSE, header = FALSE, skip = 2)
  # suffix_4m = fread(paste0(keep_path, 'dtrios_id_suffixes_mothers_', i, '.txt'),
  #                   data.table = FALSE, header = TRUE)
  # fam_4m[ , 1] = paste0(fam_4m[ , 1], suffix_4m[ , 1])
  # fam_4m[ , 2] = paste0(fam_4m[ , 2], suffix_4m[ , 2])
  # fwrite(rbind(c('ID_1', 'ID_2', 'missing', 'sex'),
  #              c(0, 0, 0, 'D'),
  #              fam_4m), paste0(outpath, out_file_sample4m, '_', i, '_suffixes.fam'),
  #        sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  # # sample 4f
  # fam_4f = fread(paste0(outpath, out_file_sample4f, '_', i, '.fam'),
  #                data.table = FALSE, header = FALSE, skip = 2)
  # suffix_4f = fread(paste0(keep_path, 'dtrios_id_suffixes_fathers_', i, '.txt'),
  #                   data.table = FALSE, header = TRUE)
  # fam_4f[ , 1] = paste0(fam_4f[ , 1], suffix_4f[ , 1])
  # fam_4f[ , 2] = paste0(fam_4f[ , 2], suffix_4f[ , 2])
  # fwrite(rbind(c('ID_1', 'ID_2', 'missing', 'sex'),
  #              c(0, 0, 0, 'D'),
  #              fam_4f), paste0(outpath, out_file_sample4f, '_', i, '_suffixes.fam'),
  #        sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
}
