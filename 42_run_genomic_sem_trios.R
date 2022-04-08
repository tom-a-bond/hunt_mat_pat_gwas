args = commandArgs(trailingOnly = TRUE)
writeLines(paste0('\n\narguments passed to run_genomic_sem_trios.R (n.b. the warnings are benign):\n\n'))
print(args)
path = args[1]
in_file = args[2]
out_path = args[3]
threads = as.numeric(args[4])
phen = args[5]
snp_list = args[6]
ref_panel = args[7]
ma_name = args[8]
chr = as.integer(args[9])
library(data.table)
#library(devtools)
#install_github("MichelNivard/GenomicSEM")
library(GenomicSEM)
setDTthreads(threads)
options(width = 150)
setwd(paste0(path))

# code adapted from step1and2.R from Helen

# STEP 0: load in GenomicSEM ####
# Create Output 
#filename <- "step1and2"
#sink(paste(filename, ".Ro", sep = ""), append = FALSE, split = TRUE)

# # STEP 1: MUNGE THE FILES ####
# n.b. I have already done this previously using the ldsc python script
# n.b. the genomic SEM authors note that their munge() function produces slightly different results to the
# the ldsc python script (and gives slightly fewer snps- see https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects)
# but this seems unlikely to affect the results so I will not remunge
# 
# files<-c(mat_file, pat_file, off_file)
# 
# #2. hm3 = the name of the reference file to use for alligning effects to same ref allele across traits
# #can be found on our github
# hm3 <- "/user/home/tb17613/software/genomic_sem/w_hm3.noMHC.snplist"
# 
# #3. trait.names = names used to create the .sumstats.gz output files
# trait.names <- paste0(phen, c('_mat', '_pat', '_off'))
# 
#4. N = total sample size for traits (will be used later)
# n_mat <- fread(paste0(path, mat_file), select = 16)
# n_pat <- fread(paste0(path, pat_file), select = 16)
# n_off <- fread(paste0(path, off_file), select = 16)
# N <- c(max(n_mat[ , 1]), max(n_pat[ , 1]), max(n_off[ , 1]))
# 
# #Run the munge function. This will create three .sumstats.gz files
# munge(files = files, hm3 = hm3, trait.names = trait.names, N = N)

# STEP 2: RUN LD-SCORE REGRESSION ####
#1. traits = the name of the .sumstats.gz traits
traits <- c(paste0(ma_name, '_sample_1_', phen, '_all_chr_1_ldhub_qc.sumstats.gz'),
            paste0(ma_name, '_sample_2_', phen, '_all_chr_1_ldhub_qc.sumstats.gz'),
            paste0(ma_name, '_sample_3_', phen, '_all_chr_1_ldhub_qc.sumstats.gz'))

#2. sample.prev = the proportion of cases to total sample size. For quantitative traits list NA
sample.prev <- c(NA, NA, NA)

#3. population.prev = the population lifetime prevalence of the traits. For quantitative traits list NA
population.prev <- c(NA, NA, NA)

##4. ld = folder of LD scores 
#can be found on our github
ld <- "/user/home/tb17613/software/genomic_sem/eur_w_ld_chr"

#5. wld = folder of LD scores
#can be found on our github
wld <- "/user/home/tb17613/software/genomic_sem/eur_w_ld_chr"

#6. trait.names = optional sixth  argument to list trait names so they can be named in your model
trait.names <- paste0(phen, c('_mat', '_pat', '_off'))

#Run the ldsc function 
writeLines(paste0('\n\nrunning multivariate LDSC:\n\n'))
TRIOS <- ldsc(traits = traits, sample.prev = sample.prev, population.prev = population.prev, ld = ld, wld = wld, trait.names = trait.names)

# we only need to save the ldsc results for one chromosome (because they should be identical for all chrs)
if(chr == 1){
  save(TRIOS, file = paste0(path, out_path, ma_name, '_genomic_sem_ldsc_results_', phen, '.rda'))
}

# code from here onwards adapted from PW_gSEM_Partition_TRIOS_step3and4_chr1.R from Helen ####

#Setup 

# filename    <- "chr1"
# sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# load("TRIOS.RData")  ## loading files created in STEP1 and STEP2

# load gwas sum stats, write temp input files for mat/pat/off ####
# subset to snps with info >=0.8 and maf >=1% in all cohorts
# n.b. the genomic sem authors recommend subsetting info >=0.6 and maf >=1%, but we wil only use info >=0.8 for MR

d <- fread(paste0(path, in_file), header = TRUE)
keep_snps <- fread(paste0(snp_list), header = TRUE) # list of snps with info >=0.8 and maf >=1% in all cohorts
d_sub <- d[which(d$rsid %in% keep_snps$rsid), ]
ind <- which(d_sub$chr == chr)
d_sub <- d_sub[ind] # subset by chr
d_m <- d_sub[which(!is.na(d_sub$beta_mat)), c('rsid', 'chr_pos_a1_a2', 'a1', 'a2', 'freq_a1_mat', 'beta_mat', 'se_mat',  'p_mat', 'n_mat')]
d_p <- d_sub[which(!is.na(d_sub$beta_pat)), c('rsid', 'chr_pos_a1_a2', 'a1', 'a2', 'freq_a1_pat', 'beta_pat', 'se_pat',  'p_pat', 'n_pat')]
d_o <- d_sub[which(!is.na(d_sub$beta_off)), c('rsid', 'chr_pos_a1_a2', 'a1', 'a2', 'freq_a1_off', 'beta_off', 'se_off',  'p_off', 'n_off')]
names(d_m) <- c('rsid', 'chr_pos_a1_a2', 'a1', 'a2', 'maf', 'effect', 'se', 'P', 'n')
names(d_p) <- c('rsid', 'chr_pos_a1_a2', 'a1', 'a2', 'maf', 'effect', 'se', 'P', 'n')
names(d_o) <- c('rsid', 'chr_pos_a1_a2', 'a1', 'a2', 'maf', 'effect', 'se', 'P', 'n')
fwrite(d_m, file = paste0(path, out_path, ma_name, '_tmp_mat_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'),
       compress = 'gzip', quote = FALSE, sep = '\t')
fwrite(d_p, file = paste0(path, out_path, ma_name, '_tmp_pat_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'),
       compress = 'gzip', quote = FALSE, sep = '\t')
fwrite(d_o, file = paste0(path, out_path, ma_name, '_tmp_off_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'),
       compress = 'gzip', quote = FALSE, sep = '\t')

# get n snps
N <- c(max(d_m$n), max(d_p$n), max(d_o$n))

# specify gneomic SEM model ####
# Running SUNMSTATS per CHR
#Takes four necessary arguments:
files <- c(paste0(path, out_path, ma_name, '_tmp_mat_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'),
           paste0(path, out_path, ma_name, '_tmp_pat_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'),
           paste0(path, out_path, ma_name, '_tmp_off_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'))

#2. ref = the name of the reference file used to obtain SNP MAF
ref = ref_panel # downloaded from genomic sem github

#3. trait.names = the name of the files to be used in (defined above)

#4. se.logit = whether the standard errors are on an logistic scale
se.logit <- c(FALSE, FALSE, FALSE)

#run the sumstats function below
writeLines(paste0('\n\ncalling GenomicSEM:sumstats() to prepare summary stats:\n\n'))
ss <- sumstats(files = files, ref = ref, trait.names = trait.names, se.logit = se.logit,
                       OLS = c(TRUE, TRUE, TRUE), linprob = NULL, N = N,
                       keep.indel = TRUE, maf.filter = 0.005, cores = threads)

# save(ss, file="ss_chr1.RData")
# write.table(ss, file = "ss_chr1.txt", row.names = FALSE, quote = FALSE)

## Check model
covstruc <- TRIOS
SNPs <- ss

model <- paste0('
Offspring =~ .5*', phen, '_mat # make a latent offspring variable that loads on Off (with 1) and MA (with .5)
Offspring =~ .5*', phen, '_pat # make a latent offspring variable that loads on Off (with 1) and PA (with .5)
Offspring =~ 1*', phen, '_off

Maternal =~ 1*', phen, '_mat # make a latent maternal variable that loads on Maternal with 1 and on Off with.5
Maternal =~ .5*', phen, '_off

Paternal =~ 1*', phen, '_pat
Paternal =~ .5*', phen, '_off

', phen, '_mat ~~ 0*', phen, '_mat + 0*', phen, '_off + 0*', phen, '_pat # no residual (co)variance
', phen, '_pat ~~ 0*', phen, '_off + 0*', phen, '_pat
', phen, '_off ~~ 0*', phen, '_off

Maternal ~~ Maternal # free variances
Paternal ~~ Paternal
Offspring ~~ Offspring

Maternal ~~ Offspring 
Paternal ~~ Offspring
Maternal ~~ Paternal 

')

#estimation = an optional third argument specifying the estimation method to use
estimation <- "DWLS"

#std.lv = optional fourth argument specifying whether variances of latent variables should be set to 1
std.lv = FALSE

#Run your model

writeLines(paste0('\n\ncalling GenomicSEM:usermodel() to fit SEM:\n\n'))
YourModel <- usermodel(covstruc = covstruc, model = model, estimation = estimation, std.lv = std.lv)

# we only need to save the model for one chromosome (because it should be identical for all chrs)
if(chr == 1){
  save(YourModel, file = paste0(path, out_path, ma_name, '_genomic_sem_usermodel_', phen, '.rda'))
}

writeLines(paste0('\n\nsummary of fitted SEM (usermodel):\n\n'))
YourModel  ## I like to check this before running the whole GWAS

## For running full GWAS 
GWAS <- paste0('
Offspring =~ .5*', phen, '_mat # make a latent offspring variable that loads on Off (with 1) and MA (with .5)
Offspring =~ .5*', phen, '_pat # make a latent offspring variable that loads on Off (with 1) and PA (with .5)
Offspring =~ 1*', phen, '_off

Maternal =~ 1*', phen, '_mat # make a latent maternal variable that loads on Maternal with 1 and on Off with.5
Maternal =~ .5*', phen, '_off

Paternal =~ 1*', phen, '_pat
Paternal =~ .5*', phen, '_off

', phen, '_mat ~~ 0*', phen, '_mat + 0*', phen, '_off + 0*', phen, '_pat # no residual (co)variance
', phen, '_pat ~~ 0*', phen, '_off + 0*', phen, '_pat
', phen, '_off ~~ 0*', phen, '_off

Maternal ~~ Maternal # free variances
Paternal ~~ Paternal
Offspring ~~ Offspring

Maternal ~~ Offspring 
Paternal ~~ Offspring
Maternal ~~ Paternal 

Maternal ~ SNP
Paternal ~ SNP
Offspring ~ SNP # regress the latent variables on the SNP
')

#estimation = an optional third argument specifying the estimation method to use
estimation <- "DWLS"
parallel <- TRUE
SNPSE <- FALSE

# run genomic SEM model ####

writeLines(paste0('\n\ncalling GenomicSEM:userGWAS() to estimate snp effects:\n\n'))
YourModel2 <- userGWAS(covstruc = covstruc, SNPs = SNPs, model = GWAS, estimation = estimation, toler=FALSE,
                       cores = threads, parallel = parallel, SNPSE = SNPSE, sub = c("Maternal~SNP", "Paternal~SNP", "Offspring~SNP"))

#save(YourModel2, file = paste0(path, ma_name, '_genomic_sem_YourModel2_', phen, '.rda'))

res_mat <- YourModel2[[1]]
res_pat <- YourModel2[[2]]
res_off <- YourModel2[[3]]

# write model info file

sink(paste0(path, out_path, ma_name, '_genomic_sem_model_summary_', phen, '_chr', chr, '.txt'))
writeLines('\nYourmodel (estimated using all chromosomes):\n')
print(YourModel)
writeLines('\ndimensions of YourModel2:\n')
writeLines('\nmaternal:\n')
dim(YourModel2[[1]])
writeLines('\npaternal:\n')
dim(YourModel2[[2]])
writeLines('\noffspring:\n')
dim(YourModel2[[3]])
writeLines('\n\ntables of warnings for YourModel2:\n')
table(res_mat$warning)
table(res_pat$warning)
table(res_off$warning)
sink()

# write gwas output

fwrite(res_mat, file = paste0(path, out_path, ma_name, '_genomic_sem_results_mat_', phen, '_chr', chr, '.txt.gz'),
       row.names = FALSE, quote = FALSE, compress = 'gzip', sep = '\t')
fwrite(res_pat, file = paste0(path, out_path, ma_name, '_genomic_sem_results_pat_', phen, '_chr', chr, '.txt.gz'),
       row.names = FALSE, quote = FALSE, compress = 'gzip', sep = '\t')
fwrite(res_off, file = paste0(path, out_path, ma_name, '_genomic_sem_results_off_', phen, '_chr', chr, '.txt.gz'),
       row.names = FALSE, quote = FALSE, compress = 'gzip', sep = '\t')

# delete tmp input files

writeLines(paste0('\n\ndeleting temporary input files\n\n'))
file.remove(paste0(path, out_path, ma_name, '_tmp_mat_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'))
file.remove(paste0(path, out_path, ma_name, '_tmp_pat_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'))
file.remove(paste0(path, out_path, ma_name, '_tmp_off_genomic_sem_input_', phen, '_chr', chr, '.tsv.gz'))



