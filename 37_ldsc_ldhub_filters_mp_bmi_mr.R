args = commandArgs(trailingOnly = TRUE)
in_file = args[1]
phen = args[2]
maf = as.numeric(args[3])
chisq = as.numeric(args[4])
library(data.table)
ss = fread(paste0(in_file, '.tbl.gz'), data.table = FALSE)
writeLines(paste0('before applying extra ldhub qc filters, nsnp = ', nrow(ss)))
# remove mhc region
if(nrow(ss[-which(ss$CHR == 6 & ss$POS >= 26000000 & ss$POS <= 34000000), ]) > 0){
  ss = ss[-which(ss$CHR == 6 & ss$POS >= 26000000 & ss$POS <= 34000000), ]
}
writeLines(paste0('after removing MHC region, nsnp = ', nrow(ss)))
# remove if maf <1%
ss$Freq1 <- as.numeric(ss$Freq1)
ss = ss[which((ss$Freq1 >= maf & ss$Freq1 <= (1 - maf)) | is.na(ss$Freq1)), ]
writeLines(paste0('after removing snps with MAF <1%, nsnp = ', nrow(ss)))
# remove if chisq >80
ss$chisq = qchisq(ss$`P-value`, df = 1, ncp = 0, lower.tail = FALSE)
ss = ss[which(ss$chisq <= chisq), ]
writeLines(paste0('after removing snps with chisq >80, nsnp = ', nrow(ss)))
# remove if not present in all cohorts meta-analysed
ss = ss[which(ss$total_n == max(ss$total_n)), ]
writeLines(paste0('after removing snps not present in all cohorts, nsnp = ', nrow(ss)))
# drop snps which are not present in ldsc hapmap3 data, use rsids not chr_pos_a1_a2 to enable merging with ldsc hapmap3 reference
# input file written by ldsc_hapmap_snps_list_chr_pos_a1_a2.R
hm <- fread('/user/home/tb17613/software/ldsc/downloaded_ld_scores/eur_w_ld_chr/w_hm3_chr_pos_a1_a2_inc_duplicates.snplist.gz',
            header = TRUE)
ss = data.table(ss)
m <- merge(ss, hm[ , c('SNP', 'chr:pos_a1/a2')], by.x = 'MarkerName', by.y = 'chr:pos_a1/a2', all.x = TRUE)
m <- merge(m, hm[ , c('SNP', 'chr:pos_a2/a1')], by.x = 'MarkerName', by.y = 'chr:pos_a2/a1', all.x = TRUE)
m$SNP <- m$SNP.y
m[which(is.na(m$SNP.y)), 'SNP'] <- m[which(is.na(m$SNP.y)), 'SNP.x']
drop_vars <- c('chisq', 'SNP.x', 'SNP.y', 'MarkerName')
w <- m
w <- w[which(!is.na(w$SNP)), ]
w <- w[ , c('chisq', 'SNP.x', 'SNP.y', 'MarkerName'):=NULL]
fwrite(w, paste0(in_file, '_tmp_ldhub_qc.tbl.gz'), row.names = FALSE,
       col.names = TRUE, sep = '\t', quote = FALSE, compress = 'gzip')
