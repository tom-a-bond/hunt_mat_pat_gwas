args = commandArgs(trailingOnly = TRUE)
covar_file = args[1]
mother_geno_file = args[2]
father_geno_file = args[3]
offs_geno_file = args[4]
link_file = args[5]
sample = args[6]
snp = args[7]
out_path = args[8]
pheno_file = args[9]
# residualise phenotype on all covariates except for other parent's and offspring's genotype, save new pheno file
writeLines('Residualising phenotype on all covariates (except parental/offpsring genotype)')
pheno = read.table(pheno_file, header = FALSE)
names(pheno)[3] = 'pheno'
covars = read.table(covar_file, header = FALSE)
pheno_m = merge(covars, pheno[c('V1', 'pheno')], by = 'V1', all.x = TRUE)
f = as.formula(paste0('pheno ~ ', paste(sprintf('V%s', 3:ncol(covars)), collapse = '+')))
model = lm(f, data = pheno_m)
writeLines('Summary of residualisation regression model:')
print(summary(model))
# remove individuals with missing values before merging residuals
resid = pheno_m
resid = resid[which(!is.na(resid$pheno)), ]
for(i in 3:ncol(covars)){
  resid = resid[which(!is.na(resid[ , paste0('V', i)])), ]
}
resid$resid = residuals(model)
write.table(resid[c('V1', 'V2', 'resid')], file = paste0(pheno_file, '.unrelated.resid'), ro = FALSE, co = FALSE, qu = FALSE)
# write covar file with other parent's and offspring's genotype
offs_geno = read.table(offs_geno_file, header = TRUE)
link = read.table(link_file, header = TRUE)
if(sample == '4m'){
  other_parent_geno = read.table(father_geno_file, header = TRUE)
  index_parent_id = 'id_m_original'
  other_parent_id = 'id_f_original'
}
if(sample == '4f'){
  other_parent_geno = read.table(mother_geno_file, header = TRUE)
  index_parent_id = 'id_f_original'
  other_parent_id = 'id_m_original'
}
other_parent_geno_m = merge(other_parent_geno, link, by.x = 'FID', by.y = other_parent_id, all.x = TRUE)
geno_varname_other_parent = names(other_parent_geno_m)[grep(gsub(':|/', '.', snp), names(other_parent_geno_m))]
geno_m = merge(covars[ , 1:2], other_parent_geno_m[, c(geno_varname_other_parent, index_parent_id)], by.x = 'V1', by.y = index_parent_id, all.x = TRUE)
offs_geno_m = merge(offs_geno, link, by.x = 'FID', by.y = 'id_o', all.x = TRUE)
geno_varname_offs = names(offs_geno_m)[grep(gsub(':|/', '.', snp), names(offs_geno_m))]
if(geno_varname_other_parent != geno_varname_offs){ stop('Effect/non-effect alleles are different in .raw files for offspring and the other (non-index) parent; is there a problem?') }
geno_m2 = merge(geno_m, offs_geno_m[, c(geno_varname_offs, index_parent_id)], by.x = 'V1', by.y = index_parent_id, all.x = TRUE)
write.table(geno_m2, paste0(out_path, 'covars_', snp, '_unrelated.covar'), ro = FALSE, co = FALSE, qu = FALSE)
