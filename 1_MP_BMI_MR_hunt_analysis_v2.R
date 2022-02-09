# Tom Bond 19/10/21
# write hunt phenotype, covariate and id list files

# set up ####

library(data.table)
library(dplyr)
library(stringr)

rm(list = ls())

# loop over phenotypes

in_path = 'C:/Users/tb17613/OneDrive - The University of Queensland/work/proj/maternal_paternal_bmi_mr/hunt/' # path to data files
out_path = in_path # path to write output files (ids lists, phenotype and covariate files)
hunt_phens = c('BMI', 'CHO', 'CRP', 'DBP', 'GLU', 'HBA', 'HDL', 'LDL', 'SBP', 'TG', 'WHR', 'BW', 'GA',
               'BMI_log', 'WHR_log', 'GLU_log', 'HBA_log', 'HDL_log', 'TG_log', 'CRP_log')
hunt_phens_standardised = c('bmi', 'chol', 'crp', 'dbp', 'glu', 'hba1c', 'hdl', 'ldl', 'sbp', 'tg', 'whr', 'bw', 'gest',
               'lnbmi', 'lnwhr', 'lnglu', 'lnhba1c', 'lnhdl', 'lntg', 'lncrp')

for(i in 1:length(hunt_phens)){
  hunt_phen = hunt_phens[i]
  hunt_phen_standardised = hunt_phens_standardised[i]
  writeLines(hunt_phen)
  
  # load analysis datasets ####
  
  dmo = read.table(paste0(in_path, 'Cleaned_', hunt_phen, '_MO.txt'), header = TRUE, stringsAsFactors = FALSE)
  dfo = read.table(paste0(in_path, 'Cleaned_', hunt_phen, '_FO.txt'), header = TRUE, stringsAsFactors = FALSE)
  dtrios = read.table(paste0(in_path, 'Cleaned_', hunt_phen, '_TRIO.txt'), header = TRUE, stringsAsFactors = FALSE)
  dall = read.table(paste0(in_path, 'Cleaned_', hunt_phen, '.txt'), header = TRUE, stringsAsFactors = FALSE)
  
  # exclude 7th and above sibs
  
  # get sib_no variable
  dmo$sib_no = as.numeric(substr(dmo$id_m, nchar(dmo$id_m), nchar(dmo$id_m)))
  dfo$sib_no = as.numeric(substr(dfo$id_f, nchar(dfo$id_f), nchar(dfo$id_f)))
  dtrios$sib_no_m = as.numeric(substr(dtrios$id_m, nchar(dtrios$id_m), nchar(dtrios$id_m)))
  dtrios$sib_no_f = as.numeric(substr(dtrios$id_f, nchar(dtrios$id_f), nchar(dtrios$id_f)))
  
  # exclude if sib_no >6
  dmo = dmo[which(dmo$sib_no <= 6), ]
  dfo = dfo[which(dfo$sib_no <= 6), ]
  dtrios = dtrios[which(dtrios$sib_no_m <= 6), ]
  dtrios = dtrios[which(dtrios$sib_no_f <= 6), ]
  
  # get original ids for parents (those present in plink/bgen files, prior to adding suffixes)
  
  dmo$id_m_original = unlist(lapply(str_split(dmo$id_m, fixed('_')), function(x){
    x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
  }))
  dfo$id_f_original = unlist(lapply(str_split(dfo$id_f, fixed('_')), function(x){
    x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
  }))
  dtrios$id_m_original = unlist(lapply(str_split(dtrios$id_m, fixed('_')), function(x){
    x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
  }))
  dtrios$id_f_original = unlist(lapply(str_split(dtrios$id_f, fixed('_')), function(x){
    x_1_3 = x[1:3]; return(paste(x_1_3, collapse = '_'))
  }))
  
  # get suffixes added to parents ids
  
  dmo$id_m_suffix = unlist(lapply(str_split(dmo$id_m, fixed('_')), function(x){
    x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
  }))
  dfo$id_f_suffix = unlist(lapply(str_split(dfo$id_f, fixed('_')), function(x){
    x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
  }))
  dtrios$id_m_suffix = unlist(lapply(str_split(dtrios$id_m, fixed('_')), function(x){
    x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
  }))
  dtrios$id_f_suffix = unlist(lapply(str_split(dtrios$id_f, fixed('_')), function(x){
    x_4_5 = x[4:5]; return(paste0('_', paste(x_4_5, collapse = '_')))
  }))

  # write id lists for subsetting plink/bgen files ####

  if(hunt_phen == 'BMI'){ # these id list files should be identical for each phenotype, so only need to be written once
    
    # id lists for creating parental genotype files
    for(j in 1:6){ # loop over sibs 1-6
      # mother-offspring pairs (sample 1)
      tmp_m = dmo[which(dmo$sib_no == j), ]
      if(nrow(tmp_m) > 0){
        write.table(tmp_m[c('id_m_original', 'id_m_original')], paste0(out_path, 'dmo_ids_mothers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
        write.table(tmp_m[c('id_m_suffix', 'id_m_suffix')], paste0(out_path, 'dmo_id_suffixes_mothers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
      }
      # father-offspring pairs (sample 2)
      tmp_f = dfo[which(dfo$sib_no == j), ]
      if(nrow(tmp_f) > 0){
        write.table(tmp_f[c('id_f_original', 'id_f_original')], paste0(out_path, 'dfo_ids_fathers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
        write.table(tmp_f[c('id_f_suffix', 'id_f_suffix')], paste0(out_path, 'dfo_id_suffixes_fathers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
      }
      # mother-father-offspring trios (sample 4m/4f)
      # mothers
      tmp_tm = dtrios[which(dtrios$sib_no_m == j), ]
      if(nrow(tmp_tm) > 0){
        write.table(tmp_tm[c('id_m_original', 'id_m_original')], paste0(out_path, 'dtrios_ids_mothers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
        write.table(tmp_tm[c('id_m_suffix', 'id_m_suffix')], paste0(out_path, 'dtrios_id_suffixes_mothers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
      }
      # fathers
      tmp_tf = dtrios[which(dtrios$sib_no_f == j), ]
      if(nrow(tmp_tf) > 0){
        write.table(tmp_tf[c('id_f_original', 'id_f_original')], paste0(out_path, 'dtrios_ids_fathers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
        write.table(tmp_tf[c('id_f_suffix', 'id_f_suffix')], paste0(out_path, 'dtrios_id_suffixes_fathers_', j, '.txt'),
                    qu = FALSE, ro = FALSE, co = FALSE)
      }
      
      # # now write id lists for reordering the plink/bgen files once files 1-6 have been merged (NO LONGER REQUIRED)
      # # mother-offspring pairs (sample 1)
      # tmp_m = dmo
      # write.table(tmp_m[c('id_m', 'id_m')], paste0(out_path, 'dmo_ids_mothers_reorder.txt'),
      #             qu = FALSE, ro = FALSE, co = FALSE)
      # # father-offspring pairs (sample 2)
      # tmp_f = dfo
      # write.table(tmp_f[c('id_f', 'id_f')], paste0(out_path, 'dfo_ids_fathers_reorder.txt'),
      #             qu = FALSE, ro = FALSE, co = FALSE)
      # # mother-father-offspring trios (sample 4m/4f)
      # # mothers
      # tmp_tm = dtrios
      # write.table(tmp_tm[c('id_m', 'id_m')], paste0(out_path, 'dtrios_ids_mothers_reorder.txt'),
      #             qu = FALSE, ro = FALSE, co = FALSE)
      # # fathers
      # tmp_tf = dtrios
      # write.table(tmp_tf[c('id_f', 'id_f')], paste0(out_path, 'dtrios_ids_fathers_reorder.txt'),
      #             qu = FALSE, ro = FALSE, co = FALSE)

    }
    # id lists for creating offspring genotype files
    # mother-father-offspring trios (sample 4m/4f)
    dtrios_ids_offspring = data.frame(ID_1 = dtrios[ , 'id_o'],
                                      ID_2 = dtrios[ , 'id_o'])
    write.table(dtrios_ids_offspring, paste0(out_path, 'dtrios_ids_offspring.txt'),
                qu = FALSE, ro = FALSE, co = FALSE)
    # full hunt sample (sample 3)
    dall_ids = data.frame(ID_1 = dall[ , 'IID'],
                          ID_2 = dall[ , 'IID'])
    write.table(dall_ids, paste0(out_path, 'dall_ids.txt'),
                qu = FALSE, ro = FALSE, co = FALSE)
    
    # linker file for mothers/fathers/offspring in sample 4m/4f
    write.table(dtrios[c('id_m', 'id_f', 'id_o')], paste0(out_path, 'dtrios_id_linker.txt'),
                qu = FALSE, ro = FALSE, co = TRUE)
    
  }
  
  # write pheno and covar files for fastgwa ####
  
  cont_covars = c('Sex', 'AGE', 'AGE2', 'AGESEX', 'AGE2SEX')
  pcs = sprintf('PC%s', 1:20)
  cat_covars = c('batch', 'PR')

  # phenotype files

  # sample 1 (mother offspring duos)
  write.table(dmo[c('id_m', 'id_m', 'PHE_o')],
              paste0(out_path, 'sample1_', hunt_phen_standardised, '.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 2 (father offspring duos)
  write.table(dfo[c('id_f', 'id_f', 'PHE_o')],
              paste0(out_path, 'sample2_', hunt_phen_standardised, '.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 3 (full sample [all treated as offspring])
  write.table(dall[c('IID', 'IID', 'PHE')],
              paste0(out_path, 'sample3_', hunt_phen_standardised, '.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 4m (trios, mother ids)
  write.table(dtrios[c('id_m', 'id_m', 'PHE_o')],
              paste0(out_path, 'sample4m_', hunt_phen_standardised, '.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 4f (trios, father ids)
  write.table(dtrios[c('id_f', 'id_f', 'PHE_o')],
              paste0(out_path, 'sample4f_', hunt_phen_standardised, '.phen'),
              qu = FALSE, ro = FALSE, co = FALSE)

  # covariate files
  # n.b. for samples 4m and 4f the offspring's and other parent's genotype will be added as covariates later
  # sample 1 (mother offspring duos)
  write.table(dmo[c('id_m', 'id_m', paste0(cont_covars, '_o'), paste0(pcs, '_m'))], # n.b. using offspring covariates but parental pcs
              paste0(out_path, 'sample1_', hunt_phen_standardised, '.qcovar'), # continuous covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  write.table(dmo[c('id_m', 'id_m', 'batch_m', 'PR_o')], # n.b. using parental genotyping batch and offspring hunt wave
              paste0(out_path, 'sample1_', hunt_phen_standardised, '.covar'), # categorical covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 2 (father offspring duos)
  write.table(dfo[c('id_f', 'id_f', paste0(cont_covars, '_o'), paste0(pcs, '_f'))], # n.b. using offspring covariates but parental pcs
              paste0(out_path, 'sample2_', hunt_phen_standardised, '.qcovar'), # continuous covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  write.table(dfo[c('id_f', 'id_f', 'batch_f', 'PR_o')], # n.b. using parental genotyping batch and offspring hunt wave
              paste0(out_path, 'sample2_', hunt_phen_standardised, '.covar'), # categorical covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 3 (full sample, offspring, mothers, fathers [all treated as offspring])
  write.table(dall[c('IID', 'IID', cont_covars, pcs)],
              paste0(out_path, 'sample3_', hunt_phen_standardised, '.qcovar'), # continuous covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  write.table(dall[c('IID', 'IID', 'batch', 'PR')],
              paste0(out_path, 'sample3_', hunt_phen_standardised, '.covar'), # categorical covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 4m (trios, mother ids)
  write.table(dtrios[c('id_m', 'id_m', paste0(cont_covars, '_o'), paste0(pcs, '_m'))], # n.b. using offspring covariates but parental pcs
              paste0(out_path, 'sample4m_', hunt_phen_standardised, '.qcovar'), # continuous covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  write.table(dtrios[c('id_m', 'id_m', 'batch_m', 'PR_o')], # n.b. using parental genotyping batch and offspring hunt wave
              paste0(out_path, 'sample4m_', hunt_phen_standardised, '.covar'), # categorical covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  # sample 4f (trios, father ids)
  write.table(dtrios[c('id_f', 'id_f', paste0(cont_covars, '_o'), paste0(pcs, '_f'))], # n.b. using offspring covariates but parental pcs
              paste0(out_path, 'sample4f_', hunt_phen_standardised, '.qcovar'), # continuous covariates
              qu = FALSE, ro = FALSE, co = FALSE)
  write.table(dtrios[c('id_f', 'id_f', 'batch_f', 'PR_o')], # n.b. using parental genotyping batch and offspring hunt wave
              paste0(out_path, 'sample4f_', hunt_phen_standardised, '.covar'), # categorical covariates
              qu = FALSE, ro = FALSE, co = FALSE)

}

