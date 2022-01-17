#!/bin/bash
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -k o

# create plink binary files for input into grm calculation

################################################################################
##### variables which require editing:
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory containing add_suffixes_plink_files.R
plink2=/home/ttbond/plink2_linux_x86_64/plink2 # plink 2 executable. Must be plink 2, not plink 1.9
plink1_9=/location/of/plink # plink 1.9 executable. Must be plink 1.9, not plink 2
gen_path=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/CleanGenotypeData/ # path to plink 1 binary fileset containing qc'd called genotypes for whole hunt sample (I assume there is only one set of plink files for the entire cohort)
gen_file=ukb_cal_chr1_22_v2_cleaned # filename stem for plink 1 binary fileset containing qc'd called genotypes for whole hunt sample (I assume there is only one set of plink files for the entire cohort)
keep_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/alspac/ # path to id lists written by MP_BMI_MR_hunt_analysis.R (i.e. "out_path" in MP_BMI_MR_hunt_analysis.R)
out_path_gen=/path/to/write/plink/files/to/ # path to directory that subsetted plink files will be written to
out_file_sample1=${cohort_name}_sample1_maternal_geno # output file name for sample 1 (marginal maternal gwas)
out_file_sample2=${cohort_name}_sample2_paternal_geno # output file name for sample 2 (marginal paternal gwas)
out_file_sample4m=${cohort_name}_sample4m_maternal_geno # output file name for sample 4m (conditional trios gwas; maternal ids)
out_file_sample4f=${cohort_name}_sample4f_paternal_geno # output file name for sample 4f (conditional trios gwas; paternal ids)
################################################################################

# subset the full hunt sample plink files (all qc'd called genotypes) to the samples required for grm calculation
# j indexes the id lists for the 1st ... 6th siblings, written by MP_BMI_MR_hunt_analysis_v2.R

for j in {1..6}; do
  # sample 1 (maternal genotype)
  ${plink2} \
  --bfile ${gen_path}${gen_file} \
  --keep ${keep_path}dmo_ids_mothers_${j}.txt \
  --make-bed \
  --out ${out_path_gen}${out_file_sample1}_${j}
  # sample 2 (paternal genotype)
  ${plink2} \
  --bfile ${gen_path}${gen_file} \
  --keep ${keep_path}dfo_ids_fathers_${j}.txt \
  --make-bed \
  --out ${out_path_gen}${out_file_sample2}_${j}
  # sample 3 does not need to be subsetted, because it will be the full hunt sample
  # sample 4m (maternal genotype)
  ${plink2} \
  --bfile ${gen_path}${gen_file} \
  --keep ${keep_path}dtrios_ids_mothers_${j}.txt \
  --make-bed \
  --out ${out_path_gen}${out_file_sample4m}_${j}
  # sample 4f (paternal genotype)
  ${plink2} \
  --bfile ${gen_path}${gen_file} \
  --keep ${keep_path}dtrios_ids_fathers_${j}.txt \
  --make-bed \
  --out ${out_path_gen}${out_file_sample4f}_${j}
done

# PLEASE CHECK THAT THE .fam FILES WRITTEN ABOVE HAVE THE CORRECT NUMBER OF IDS (I.E. THE SAME AS IN THE ID LISTS [dmo_ids_mothers_${j}.txt etc.])

# add suffixes to parental ids

export MODULEPATH=/sw/el7/modulefiles; module load R/4.0.5 # put any commands required to get Rscript.exe "in path" here
Rscript ${script_path}add_suffixes_plink_files.R \
${keep_path} \
${out_path_gen} \
${out_file_sample1} \
${out_file_sample2} \
${out_file_sample4m} \
${out_file_sample4f}

# merge parental files for sibs 1-6, and reorder parental files so that they match up with the offspring files

# sample 1
# write list of files to merge
echo ${out_path_gen}${out_file_sample1}_1.bed ${out_path_gen}${out_file_sample1}_1.bim ${out_path_gen}${out_file_sample1}_1_suffixes.fam > \
${out_path_gen}plink_files_list_${out_file_sample1}.txt
for i in {2..6}; do
  echo ${out_path_gen}${out_file_sample1}_${i}.bed ${out_path_gen}${out_file_sample1}_${i}.bim ${out_path_gen}${out_file_sample1}_${i}_suffixes.fam >> \
  ${out_path_gen}plink_files_list_${out_file_sample1}.txt
done
# merge
${plink1_9} \
--merge-list ${out_path_gen}plink_files_list_${out_file_sample1}.txt \
--indiv-sort file ${keep_path}dmo_ids_mothers_reorder.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample1}

# sample 2
# write list of files to merge
echo ${out_path_gen}${out_file_sample2}_1.bed ${out_path_gen}${out_file_sample2}_1.bim ${out_path_gen}${out_file_sample2}_1_suffixes.fam > \
${out_path_gen}plink_files_list_${out_file_sample2}.txt
for i in {2..6}; do
  echo ${out_path_gen}${out_file_sample2}_${i}.bed ${out_path_gen}${out_file_sample2}_${i}.bim ${out_path_gen}${out_file_sample2}_${i}_suffixes.fam >> \
  ${out_path_gen}plink_files_list_${out_file_sample2}.txt
done
# merge
${plink1_9} \
--merge-list ${out_path_gen}plink_files_list_${out_file_sample2}.txt \
--indiv-sort file ${keep_path}dfo_ids_fathers_reorder.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample2}

# sample 3 does not need to be subsetted

# sample 4m
# write list of files to merge
echo ${out_path_gen}${out_file_sample4m}_1.bed ${out_path_gen}${out_file_sample4m}_1.bim ${out_path_gen}${out_file_sample4m}_1_suffixes.fam > \
${out_path_gen}plink_files_list_${out_file_sample4m}.txt
for i in {2..6}; do
  echo ${out_path_gen}${out_file_sample4m}_${i}.bed ${out_path_gen}${out_file_sample4m}_${i}.bim ${out_path_gen}${out_file_sample4m}_${i}_suffixes.fam >> \
  ${out_path_gen}plink_files_list_${out_file_sample4m}.txt
done
# merge
${plink1_9} \
--merge-list ${out_path_gen}plink_files_list_${out_file_sample4m}.txt \
--indiv-sort file ${keep_path}dtrios_ids_mothers_reorder.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample4m}

# sample 4f
# write list of files to merge
echo ${out_path_gen}${out_file_sample4f}_1.bed ${out_path_gen}${out_file_sample4f}_1.bim ${out_path_gen}${out_file_sample4f}_1_suffixes.fam > \
${out_path_gen}plink_files_list_${out_file_sample4f}.txt
for i in {2..6}; do
  echo ${out_path_gen}${out_file_sample4f}_${i}.bed ${out_path_gen}${out_file_sample4f}_${i}.bim ${out_path_gen}${out_file_sample4f}_${i}_suffixes.fam >> \
  ${out_path_gen}plink_files_list_${out_file_sample4f}.txt
done
# merge
${plink1_9} \
--merge-list ${out_path_gen}plink_files_list_${out_file_sample4f}.txt \
--indiv-sort file ${keep_path}dtrios_ids_fathers_reorder.txt \
--make-bed \
--out ${out_path_gen}${out_file_sample4f}

# PLEASE CHECK THAT THE .fam FILES PRODUCED IMMEDIATELY ABOVE HAVE THE CORRECT NUMBER OF IDS
# (I.E. IT SHOULD BE THE SAME NUMBER OF IDS AS IN dmo, dfo AND dtrios RESPECTIVELY)


