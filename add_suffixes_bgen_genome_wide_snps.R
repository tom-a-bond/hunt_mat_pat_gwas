args = commandArgs(trailingOnly = TRUE)
keep_path = args[1]
outpath = args[2]
bgen_file_name = args[3]
suffixes_file_name = args[4]
library(data.table)
# mothers
sample = fread(paste0(outpath, bgen_file_name, '.sample'),
               data.table = FALSE, header = FALSE, skip = 2)
suffix = fread(paste0(keep_path, suffixes_file_name),
                 data.table = FALSE, header = FALSE)
sample[ , 1] = paste0(sample[ , 1], suffix[ , 1])
sample[ , 2] = paste0(sample[ , 2], suffix[ , 2])
fwrite(rbind(c('ID_1', 'ID_2', 'missing', 'sex'), # n.b. the sample files written by qctool do not contain a 4th column (sex), so this will not be added here and will need to be added later
             c(0, 0, 0, 'D'),
             sample), paste0(outpath, bgen_file_name, '_suffixes.sample'),
       sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
