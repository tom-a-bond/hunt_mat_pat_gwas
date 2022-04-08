args = commandArgs(trailingOnly = TRUE)
sample = args[1]
phen = args[2]
file_end = args[3]
in_path = args[4]
ma_name = args[5]
library(data.table)
library(stringr)
in_file = paste0(in_path, ma_name, '/', ma_name, '_sample_', sample, '_', phen, '_all_', file_end, '_1.tbl.gz')
inp = fread(in_file, header = TRUE, data.table = FALSE)
str_m = str_split_fixed(inp$MarkerName, fixed('_'), n = 2)
chr_pos = str_split_fixed(str_m[, 1], fixed(':'), n = 2)
inp$CHR = chr_pos[ , 1]
inp$POS = chr_pos[ , 2]
fwrite(inp, file = in_file,
       sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE, compress = 'gzip')
