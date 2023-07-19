args = commandArgs(trailingOnly = TRUE)
path = args[1]
mat_file = args[2]
pat_file = args[3]
off_file = args[4]
mat_off_ldsc_file = args[5]
pat_off_ldsc_file = args[6]
mat_pat_ldsc_file = args[7]
out_file = args[8]
threads = as.numeric(args[9])
phen = args[10]
plot_path = args[11]
alpha = as.numeric(args[12])
extra_text = ifelse(length(args) > 12, args[13], '')
library(data.table)
library(DONUTS)
library(dplyr) # for some reason this is not loaded when loading DONUTS
setDTthreads(threads)
mat_all = fread(paste0(path, mat_file), header = TRUE)
#pat_all = fread(paste0(path, pat_file), header = TRUE)
off_all = fread(paste0(path, off_file), header = TRUE)
writeLines(paste0('maternal gwas nsnp = ', nrow(mat_all)))
#writeLines(paste0('paternal gwas nsnp = ', nrow(pat_all)))
writeLines(paste0('offspring gwas nsnp = ', nrow(off_all)))
# drop snps which aren't available in all cohorts (to avoid complications with ldsc sample overlap calculations)
mat = mat_all[which(mat_all$total_n == max(mat_all$total_n, na.rm = TRUE)), ]
#pat = pat_all[which(pat_all$total_n == max(pat_all$total_n, na.rm = TRUE)), ]
off = off_all[which(off_all$total_n == max(off_all$total_n, na.rm = TRUE)), ]
writeLines(paste0('after keeping only snps available in all cohorts, maternal gwas nsnp = ', nrow(mat)))
#writeLines(paste0('after keeping only snps available in all cohorts, paternal gwas nsnp = ', nrow(pat)))
writeLines(paste0('after keeping only snps available in all cohorts, offspring gwas nsnp = ', nrow(off)))
# drop irrelevant columns, make new chr_pos_a1_a2 column accurately reflecting the tested allele, rename columns for merge
format_ss <- function(x, mat_pat_off){ # 'mat_pat_off' = 'mat'/'pat'/'off'
  x$A1 <- toupper(x$Allele1)
  x$A2 <- toupper(x$Allele2)
  x$BP <- x$POS
  x$N <- x$total_n
  x$SNP <- paste0(x$CHR, ':', x$POS, '_', x$A1, '/', x$A2)
  x[ , c('MinFreq', 'MaxFreq', 'Direction', 'Allele1', 'Allele2')] <- NULL
  names(x)[1] <- 'MarkerName_old'
  names(x)[which(names(x) == 'Effect')] <- 'BETA'
  names(x)[which(names(x) == 'StdErr')] <- 'SE'
  names(x)[which(names(x) == 'P-value')] <- 'P'
  return(x)
}
mat <- format_ss(mat, 'mat')
#pat <- format_ss(pat, 'pat')
off <- format_ss(off, 'off')
# the check for discrepant alleles has already been done when running the wlm previously:
# # there should not be any discrepant alleles:
# if(nrow(m[which(m$A1_mat != m$A1_pat | m$A1_mat != m$A1_off | m$A1_pat != m$A1_off), ]) > 0){
#   stop('alleles appear to be inconsistent between maternal, paternal and/or offspring gwas')
# }
# read in ldsc intercepts
mat_off_ldsc_read <- readLines(paste0(path, mat_off_ldsc_file))
mat_off_ldsc <- scan(text = mat_off_ldsc_read[length(mat_off_ldsc_read) - 3], what = 'character')[]
cm <- as.numeric(mat_off_ldsc[11])
cm_se <- as.numeric(mat_off_ldsc[12])
# pat_off_ldsc_read <- readLines(paste0(path, pat_off_ldsc_file))
# pat_off_ldsc <- scan(text = pat_off_ldsc_read[length(pat_off_ldsc_read) - 3], what = 'character')[]
# cp <- as.numeric(pat_off_ldsc[11])
# cp_se <- as.numeric(pat_off_ldsc[12])
# mat_pat_ldsc_read <- readLines(paste0(path, mat_pat_ldsc_file))
# mat_pat_ldsc <- scan(text = mat_pat_ldsc_read[length(mat_pat_ldsc_read) - 3], what = 'character')[]
# mp <- as.numeric(mat_pat_ldsc[11])
# mp_se <- as.numeric(mat_pat_ldsc[12])
m_est <- mat_off_ldsc_read[grep('Heritability of phenotype 1', mat_off_ldsc_read) + 5]
o_est <- mat_off_ldsc_read[grep('Heritability of phenotype 2', mat_off_ldsc_read) + 5]
#p_est <- pat_off_ldsc_read[grep('Heritability of phenotype 1', pat_off_ldsc_read) + 5]
writeLines(paste0('child-mother LDSC intercept: ', cm, ', SE: ', cm_se))
#writeLines(paste0('child-father LDSC intercept: ', cp, ', SE: ', cp_se))
#writeLines(paste0('mother-father LDSC intercept: ', mp, ', SE: ', mp_se))
writeLines(paste0('child LDSC ', o_est))
writeLines(paste0('mother LDSC ', m_est))
#writeLines(paste0('father LDSC ', p_est))
# ldsc intercepts files have already been written previously:
# sink(paste0(path, 'wlm_ldsc_intercepts_', phen, extra_text, '.txt'))
# writeLines(paste0('child-mother LDSC intercept: ', cm, ', SE: ', cm_se))
# writeLines(paste0('child-father LDSC intercept: ', cp, ', SE: ', cp_se))
# writeLines(paste0('mother-father LDSC intercept: ', mp, ', SE: ', mp_se))
# writeLines(paste0('child LDSC ', o_est))
# writeLines(paste0('mother LDSC ', m_est))
# writeLines(paste0('father LDSC ', p_est))
# sink()
# throw error if any LDSC intercepts are missing
if(is.na(cm)){
  stop('one or more LDSC intercepts appear to be missing')
}

# apply donuts
donuts <- donuts(ss.own = off, ss2 = mat, mode = 3, 
                 l12 = cm, alpha = alpha,
                 OutDir = NULL)  

# add rsids from ukb pvar file
ukb <- fread('/foo/UkBbF-ImpRaw-Chr_all.pvar.chr_pos_ref_alt.gz',
             header = TRUE)
donuts <- merge(donuts, ukb[ , c('chr_pos_ref_alt', 'ID')], by.x = 'SNP', by.y = 'chr_pos_ref_alt', all.x = TRUE)
donuts <- merge(donuts, ukb[ , c('chr_pos_alt_ref', 'ID')], by.x = 'SNP', by.y = 'chr_pos_alt_ref', all.x = TRUE)
donuts$rsid <- donuts$ID.y
donuts[which(is.na(donuts$rsid)), 'rsid'] <- donuts[which(is.na(donuts$rsid)), 'ID.x']
donuts[ , c('ID.y', 'ID.x')] <- NULL
# tidy up output
donuts <- donuts[ , c("rsid", "SNP", "CHR", "BP", "A1", "A2",
                      "beta.ss2", "se.ss2", "p.ss2", "n.ss2",
                      "beta.own", "se.own", "p.own", "n.own",
                      "beta.ind.ss2", "se.ind.ss2", "p.ind.ss2", "n.ind.ss2",
                      "beta.dir", "se.dir", "p.dir", "n.dir")]
names(donuts) <- c("rsid", "chr_pos_a1_a2", "chr", "pos", "a1", "a2",
                   "beta_mat", "se_mat", "p_mat", "n_mat",
                   "beta_off", "se_off", "p_off", "n_off",
                   "beta_mat_donuts", "se_mat_donuts", "p_mat_donuts", "n_mat_donuts",
                   "beta_off_donuts", "se_off_donuts", "p_off_donuts", "n_off_donuts")
# merge in other relevant cols
cols <- c("SNP", "Freq1", "FreqSE", "HetISq", "HetChiSq",
          "HetDf", "HetPVal", "total_n", 'A1', 'A2')
mat_m <- data.frame(mat)[cols]
#pat_m <- data.frame(pat)[cols]
off_m <- data.frame(off)[cols]
names(mat_m) <- paste0(names(mat_m), '_mat')
#names(pat_m) <- paste0(names(pat_m), '_pat')
names(off_m) <- paste0(names(off_m), '_off')
m <- merge(mat_m, off_m, by.x = 'SNP_mat', by.y = 'SNP_off', all.x = TRUE)
#m <- merge(m, off_m, by.x = 'SNP_mat', by.y = 'SNP_off', all.x = TRUE)
m <- m[which(!is.na(m$Freq1_mat) & !is.na(m$Freq1_off)), ]
# there should not be any discrepant alleles:
if(nrow(m[which(m$A1_mat != m$A1_off), ]) > 0){
  stop('alleles appear to be inconsistent between maternal, paternal and/or offspring gwas')
}
names(m) <- c('SNP',
              paste0(c('freq_a1', 'freq_a1_se', 'het_isq', 'het_chisq', 'het_df', 'het_p', 'n', 'a1', 'a2'), '_mat'),
              paste0(c('freq_a1', 'freq_a1_se', 'het_isq', 'het_chisq', 'het_df', 'het_p', 'n', 'a1', 'a2'), '_off'))
donuts$n_mat <- NULL
#donuts$n_pat <- NULL
donuts$n_off <- NULL
donuts <- merge(donuts, m, by.x = 'chr_pos_a1_a2', by.y = 'SNP', all.x = TRUE)
donuts$freq_a1 <- donuts$freq_a1_mat
donuts$freq_a1_se <- donuts$freq_a1_se_mat
# # allele freqs are correct:      
# nrow(donuts[which(donuts$a1 != donuts$a1_mat | donuts$a1 != donuts$a1_pat | donuts$a1 != donuts$a1_off |
#                     donuts$a2 != donuts$a2_mat | donuts$a2 != donuts$a2_pat | donuts$a2 != donuts$a2_off), ]) # 0
donuts <- donuts[ , c("rsid", "chr_pos_a1_a2", "chr", "pos", "a1", "a2", 'freq_a1', 'freq_a1_se',
                      "beta_mat", "se_mat", "p_mat", paste0(c('het_chisq', 'het_df', 'het_p', 'n'), '_mat'),
                      "beta_off", "se_off", "p_off", paste0(c('het_chisq', 'het_df', 'het_p', 'n'), '_off'),
                      "beta_mat_donuts", "se_mat_donuts", "p_mat_donuts", "n_mat_donuts",
                      "beta_off_donuts", "se_off_donuts", "p_off_donuts", "n_off_donuts")]
# write output file
alpha_name <- gsub('.', '_', alpha, fixed = TRUE) # swap . for _ if required
fwrite(donuts, file = paste0(path, 'donuts/', out_file, '_alpha_', alpha_name, '.tsv.gz'), sep = '\t', col.names = TRUE, quote = FALSE, compress = 'gzip')



