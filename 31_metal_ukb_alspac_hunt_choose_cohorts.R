args = commandArgs(trailingOnly = TRUE)
phen = args[1]
cohorts = character()
i = 2
j = 1
while(i <= length(args)){ cohorts[j] = args[i]; j = j + 1; i = i + 1 } # read in 2nd to last elements of args to cohorts vector
library(readxl)
m = as.matrix(read_xlsx('/user/work/tb17613/proj/mat_pat_bmi_mr/results/gwas_ma/input/compare_ma_input_phenos.xlsx'))
m[is.na(m)] = 'FALSE'
s = m[which(m[ , 1] == phen), ]
cohorts_all = names(s)[2:4][as.logical(s[2:4])] # all cohorts with this pheno available
cohorts_sub = intersect(cohorts_all, cohorts) # only cohorts present in the current MA
cat(cohorts_sub) # print cohorts to be run to std out, enabling saving in bash array
