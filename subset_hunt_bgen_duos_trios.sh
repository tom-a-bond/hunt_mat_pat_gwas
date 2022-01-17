#!/bin/bash
#PBS -l select=1:ncpus=20:mem=200gb
#PBS -l walltime=08:00:00
#PBS -J 1-22
#PBS -j oe
#PBS -k o

bgenpath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Data/Genetic_v3/
samplepath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/GeneticLinkingFiles/
outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/
inpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/
qctool=/home/ttbond/qctool/build/release/qctool_v2.0.8
plink2=/home/ttbond/plink2_linux_x86_64/plink2 # plink 2 executable

# n.b. the cryptic relatedness (only a handful of individuals) has been removed previously, as have individuals who failed genotype qc and non-white brit individuals, in MP_BMI_MR_ukb_define_samples_clean_data.R)
# subset individuals using qctool
# write to 1 set of bgens for duos + trios first, having excluded 2021 withdrawal of consent; n.b. exclusion of individuals who have withdrawn consent is now also being done in MP_BMI_MR_ukb_define_samples_clean_data.R, so the exclusion done here is redundant

echo chr ${PBS_ARRAY_INDEX}
# write to 1 set of bgens per chr for duos + trios first
# n.b. exclusion of withdrawn consent is now being done in MP_BMI_MR_ukb_define_samples_clean_data.R, so the exclusion here is redundant
# n.b. qctool was much slower than plink; plink only preserves dosages, not genotype probabilities, but this doesn't matter
# n.b. these files are not being used now so have been deleted
# n.b. if using this again, need to sort out id-paste modifier for -export bgen; the current version produces sample files with ids that aren't consistent with the bgen files
${plink2} \
--bgen ${bgenpath}ukb_imp_chr"${PBS_ARRAY_INDEX}"_v3.bgen ref-first \
--sample "${samplepath}"ukb53641_imp_chr1_v3_s487298.sample \
--keep ${inpath}dmo_dfo_dtrios_ids.txt \
--remove /dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/ExclusionLists/w53641_20210201.csv \
-threads 20 \
--export bgen-1.3 \
--out "${outpath}"ukb_duos_trios_chr${PBS_ARRAY_INDEX}



