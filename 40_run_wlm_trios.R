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
library(data.table)
setDTthreads(threads)
mat_all = fread(paste0(path, mat_file), header = TRUE)
pat_all = fread(paste0(path, pat_file), header = TRUE)
off_all = fread(paste0(path, off_file), header = TRUE)
writeLines(paste0('maternal gwas nsnp = ', nrow(mat_all)))
writeLines(paste0('paternal gwas nsnp = ', nrow(pat_all)))
writeLines(paste0('offspring gwas nsnp = ', nrow(off_all)))
# drop snps which aren't available in all cohorts (to avoid complications with ldsc sample overlap calculations)
mat = mat_all[which(mat_all$total_n == max(mat_all$total_n, na.rm = TRUE)), ]
pat = pat_all[which(pat_all$total_n == max(pat_all$total_n, na.rm = TRUE)), ]
off = off_all[which(off_all$total_n == max(off_all$total_n, na.rm = TRUE)), ]
writeLines(paste0('after keeping only snps available in all cohorts, maternal gwas nsnp = ', nrow(mat)))
writeLines(paste0('after keeping only snps available in all cohorts, paternal gwas nsnp = ', nrow(pat)))
writeLines(paste0('after keeping only snps available in all cohorts, offspring gwas nsnp = ', nrow(off)))
# drop irrelevant columns, make new chr_pos_a1_a2 column accurately reflecting the tested allele, rename columns for merge
format_ss <- function(x, mat_pat_off){ # 'mat_pat_off' = 'mat'/'pat'/'off'
  x$A1 <- toupper(x$Allele1)
  x$A2 <- toupper(x$Allele2)
  x[ , c('MinFreq', 'MaxFreq', 'Direction', 'Allele1', 'Allele2')] <- NULL
  x$chr_pos_a1_a2 <- paste0(x$CHR, ':', x$POS, '_', x$A1, '/', x$A2)
  names(x)[1] <- 'MarkerName_old'
  names(x)[1:15] <- paste0(names(x)[1:15], '_', mat_pat_off)
  return(x)
}
mat <- format_ss(mat, 'mat')
pat <- format_ss(pat, 'pat')
off <- format_ss(off, 'off')
# merge on chr_pos_a1_a2
# (merging on rsid results in a handful more successfully merged variants than merging on chr_pos_a1_a2 or MarkerName_old,
# but at the cost of introducing some allele discrepancies between maternal/paternal/offspring gwas)
m <- merge(mat, pat, by = 'chr_pos_a1_a2', all.x = TRUE)
m <- merge(m, off, by = 'chr_pos_a1_a2', all.x = TRUE)
# there should not be any discrepant alleles:
if(nrow(m[which(m$A1_mat != m$A1_pat | m$A1_mat != m$A1_off | m$A1_pat != m$A1_off), ]) > 0){
  stop('alleles appear to be inconsistent between maternal, paternal and/or offspring gwas')
}
# read in ldsc intercepts
mat_off_ldsc_read <- readLines(paste0(path, mat_off_ldsc_file))
mat_off_ldsc <- scan(text = mat_off_ldsc_read[length(mat_off_ldsc_read) - 3], what = 'character')[]
cm <- as.numeric(mat_off_ldsc[11])
pat_off_ldsc_read <- readLines(paste0(path, pat_off_ldsc_file))
pat_off_ldsc <- scan(text = pat_off_ldsc_read[length(pat_off_ldsc_read) - 3], what = 'character')[]
cp <- as.numeric(pat_off_ldsc[11])
mat_pat_ldsc_read <- readLines(paste0(path, mat_pat_ldsc_file))
mat_pat_ldsc <- scan(text = mat_pat_ldsc_read[length(mat_pat_ldsc_read) - 3], what = 'character')[]
mp <- as.numeric(mat_pat_ldsc[11])
writeLines(paste0('child-mother LDSC intercept: ', cm))
writeLines(paste0('child-father LDSC intercept: ', cp))
writeLines(paste0('mother-father LDSC intercept: ', mp))
sink(paste0(path, 'wlm_ldsc_intercepts_', phen, '.txt'))
writeLines(paste0('child-mother LDSC intercept: ', cm))
writeLines(paste0('child-father LDSC intercept: ', cp))
writeLines(paste0('mother-father LDSC intercept: ', mp))
sink()
# throw error if any LDSC intercepts are missing
if(!is.numeric(cm) | !is.numeric(cp) | !is.numeric(mp)){
  stop('one or more LDSC intercepts appear to be missing')
}
# define wlm function; adapted from Rob Beaumont's apply_wlm.R, with fixed typos
apply_wlm_trio<-function(data,cm,cp,mp){
  # trio
  data$c_wlm_beta<-2*data$c_beta-data$m_beta-data$p_beta
  data$c_wlm_se<-sqrt(4*data$c_se**2+data$m_se**2+data$p_se**2-4*cm*sqrt((data$c_se**2)*(data$m_se**2))-4*cp*sqrt((data$c_se**2)*(data$p_se**2)))+2*mp*sqrt((data$m_se**2)*(data$p_se**2))
  data$c_wlm_p<-2*pnorm(-abs(data$c_wlm_beta/data$c_wlm_se))
  data$m_wlm_beta<-(3*data$m_beta-2*data$c_beta+data$p_beta)/2
  data$m_wlm_se<-sqrt((9/4)*data$m_se**2+data$c_se**2+(1/4)*data$p_se**2-3*cm*sqrt((data$c_se**2)*(data$m_se**2))-cp*sqrt((data$c_se**2)*(data$p_se**2)))+1.5*mp*sqrt((data$m_se**2)*(data$p_se**2))
  data$m_wlm_p<-2*pnorm(-abs(data$m_wlm_beta/data$m_wlm_se))
  data$p_wlm_beta<-(3*data$p_beta-2*data$c_beta+data$m_beta)/2
  data$p_wlm_se<-sqrt((9/4)*data$p_se**2+data$c_se**2+(1/4)*data$m_se**2-3*cp*sqrt((data$c_se**2)*(data$p_se**2))-cm*sqrt((data$c_se**2)*(data$m_se**2)))+1.5*mp*sqrt((data$m_se**2)*(data$p_se**2))
  data$p_wlm_p<-2*pnorm(-abs(data$p_wlm_beta/data$p_wlm_se))
  # Return data frame
  data
}
# rename columns for wlm function
wlm <- m
names(wlm)[which(names(wlm) == 'Effect_off')] <- 'c_beta'
names(wlm)[which(names(wlm) == 'Effect_mat')] <- 'm_beta'
names(wlm)[which(names(wlm) == 'Effect_pat')] <- 'p_beta'
names(wlm)[which(names(wlm) == 'StdErr_off')] <- 'c_se'
names(wlm)[which(names(wlm) == 'StdErr_mat')] <- 'm_se'
names(wlm)[which(names(wlm) == 'StdErr_pat')] <- 'p_se'
wlm <- apply_wlm_trio(wlm, cm = cm, cp = cp, mp = mp)
# add rsids from ukb pvar file
# ukb <- fread('/user/work/tb17613/data/ukb/UkBbF-ImpRaw-Chr_all.pvar.gz',
#              header = TRUE)
# ukb$chr_pos_ref_alt <- paste0(ukb$`#CHROM`, ':', ukb$POS, '_', ukb$REF, '/', ukb$ALT)
# ukb$chr_pos_alt_ref <- paste0(ukb$`#CHROM`, ':', ukb$POS, '_', ukb$ALT, '/', ukb$REF)
# fwrite(ukb, file = '/user/work/tb17613/data/ukb/UkBbF-ImpRaw-Chr_all.pvar.chr_pos_ref_alt.gz',
#        col.names = TRUE, sep = '\t', quote = FALSE, compress = 'gzip')
ukb <- fread('/user/work/tb17613/data/ukb/UkBbF-ImpRaw-Chr_all.pvar.chr_pos_ref_alt.gz',
             header = TRUE)
wlm <- merge(wlm, ukb[ , c('chr_pos_ref_alt', 'ID')], by.x = 'chr_pos_a1_a2', by.y = 'chr_pos_ref_alt', all.x = TRUE)
wlm <- merge(wlm, ukb[ , c('chr_pos_alt_ref', 'ID')], by.x = 'chr_pos_a1_a2', by.y = 'chr_pos_alt_ref', all.x = TRUE)
wlm$rsid <- wlm$ID.y
wlm[which(is.na(wlm$rsid)), 'rsid'] <- wlm[which(is.na(wlm$rsid)), 'ID.x']
wlm[ , c('ID.y', 'ID.x')] <- NULL
# tidy up output
wlm <- wlm[ , c("rsid", "chr_pos_a1_a2", "CHR_mat", "POS_mat", "A1_mat", "A2_mat",
             "Freq1_mat", "FreqSE_mat", "m_beta", "m_se", "P-value_mat",
             "HetISq_mat", "HetChiSq_mat", "HetDf_mat", "HetPVal_mat", "total_n_mat",
             "Freq1_pat", "FreqSE_pat", "p_beta", "p_se", "P-value_pat",
             "HetISq_pat", "HetChiSq_pat", "HetDf_pat", "HetPVal_pat", "total_n_pat",
             "Freq1_off", "FreqSE_off", "c_beta", "c_se", "P-value_off",
             "HetISq_off", "HetChiSq_off", "HetDf_off", "HetPVal_off", "total_n_off",
             "m_wlm_beta", "m_wlm_se", "m_wlm_p",
             "p_wlm_beta", "p_wlm_se", "p_wlm_p",
             "c_wlm_beta", "c_wlm_se", "c_wlm_p")]
names(wlm) <- c("rsid", "chr_pos_a1_a2", "chr", "pos", "a1", "a2",
                "freq_a1_mat", "freq_a1_se_mat", "beta_mat", "se_mat", "p_mat",
                "het_i_sq_mat", "het_chi_sq_mat", "het_df_mat", "het_p_mat", "n_mat",
                "freq_a1_pat", "freq_a1_se_pat", "beta_pat", "se_pat", "p_pat",
                "het_i_sq_pat", "het_chi_sq_pat", "het_df_pat", "het_p_pat", "n_pat",
                "freq_a1_off", "freq_a1_se_off", "beta_off", "se_off", "p_off",
                "het_i_sq_off", "het_chi_sq_off", "het_df_off", "het_p_off", "n_off",
                "beta_mat_wlm", "se_mat_wlm", "p_mat_wlm",
                "beta_pat_wlm", "se_pat_wlm", "p_pat_wlm",
                "beta_off_wlm", "se_off_wlm", "p_off_wlm")
# plots comparing allele freqs for mat/pat/off alleles
freq_plot <- function(x, y, xlab, ylab, plot_path, phen){
  name <- paste0('SCATTER_', ylab, '_VS_', xlab, '_', phen)
  title <- paste0(ylab, '_vs_', xlab)
  writeLines(name)
  png(paste0(plot_path, name, '.png'), width = 20, height = 20, units = 'cm', res = 200)
  plot(x, y, col = rgb(0, 0, 0, alpha = 0.1), cex = 0.1, pch = 1,
       xlab = xlab, ylab = ylab,
       main = title)
  dev.off()
}
freq_plot(wlm$freq_a1_mat, wlm$freq_a1_pat, xlab = 'maternal_EAF', ylab = 'paternal_EAF', plot_path = plot_path, phen = phen)
freq_plot(wlm$freq_a1_mat, wlm$freq_a1_off, xlab = 'maternal_EAF', ylab = 'offspring_EAF', plot_path = plot_path, phen = phen)
freq_plot(wlm$freq_a1_pat, wlm$freq_a1_off, xlab = 'paternal_EAF', ylab = 'offspring_EAF', plot_path = plot_path, phen = phen)
# write output file
fwrite(wlm, file = paste0(path, out_file), sep = '\t', col.names = TRUE, quote = FALSE, compress = 'gzip')



