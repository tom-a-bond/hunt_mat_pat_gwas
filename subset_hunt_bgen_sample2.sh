#!/bin/bash
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=72:00:00
#PBS -J 1-22
#PBS -j oe
#PBS -k o

# subset bgen files with genome wide snps, write bgen files for sample 2 gwas
# n.b. I haven't run this in ukb/alspac, so I'm not sure how long/how much memory it will take

bgenpath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Data/Genetic_v3/ # path to bgen files
samplepath=/dmf/uqdi/Genomic_Medicine/Evans_Group/_BIOBANK/Biobank_53641/GeneticLinkingFiles/ # path to bgen .sample file
outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # folder to write output (subsetted bgen files) to
keep_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/alspac/ # path to id lists written by MP_BMI_MR_hunt_analysis.R (i.e. "out_path" in MP_BMI_MR_hunt_analysis.R)
qctool=/home/ttbond/qctool/build/release/qctool_v2.0.8 # path to qc tool executable
script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory containing add_suffixes_bgen_genome_wide_snps.R

# subset the full hunt sample bgen files to the samples required for sample1/2 fastgwa analysis
# j indexes the id lists for the 1st ... 6th siblings, written by MP_BMI_MR_hunt_analysis_v2.R

for j in {1..6}; do
  # sample 1
  ${qctool} \
  -g ${bgenpath}chr${PBS_ARRAY_INDEX}.bgen \
  -s ${bgenpath}chr${PBS_ARRAY_INDEX}.sample \
  -incl-samples ${keep_path}dfo_ids_fathers_${j}.txt \
  -bgen-omit-sample-identifier-block \
  -og ${outpath}hunt_sample2_fathers_${j}_chr${PBS_ARRAY_INDEX}.bgen \
  -os ${outpath}hunt_sample2_fathers_${j}_chr${PBS_ARRAY_INDEX}.sample
done

# PLEASE CHECK THAT THE SAMPLE FILES WRITTEN ABOVE HAVE THE CORRECT NUMBER OF IDS (I.E. THE SAME AS IN THE ID LISTS [dtrios_ids_mothers_${j}.txt etc.])

# add suffixes to parental ids

export MODULEPATH=/sw/el7/modulefiles; module load R/4.0.5 # put any commands required to get Rscript.exe "in path" here
for j in {1..6}; do
  echo adding suffixes to fam file ${j}
  Rscript ${script_path}add_suffixes_bgen_genome_wide_snps.R \
  ${keep_path} \
  ${outpath} \
  hunt_sample2_fathers_${j}_chr${PBS_ARRAY_INDEX} \
  dfo_id_suffixes_fathers_${j}.txt
done

# merge parental files for sibs 1-6

# mothers
${qctool} \
-g ${outpath}hunt_sample2_fathers_1_chr${PBS_ARRAY_INDEX}.bgen \
-s ${outpath}hunt_sample2_fathers_1_chr${PBS_ARRAY_INDEX}_suffixes.sample \
-g ${outpath}hunt_sample2_fathers_2_chr${PBS_ARRAY_INDEX}.bgen \
-s ${outpath}hunt_sample2_fathers_2_chr${PBS_ARRAY_INDEX}_suffixes.sample \
-g ${outpath}hunt_sample2_fathers_3_chr${PBS_ARRAY_INDEX}.bgen \
-s ${outpath}hunt_sample2_fathers_3_chr${PBS_ARRAY_INDEX}_suffixes.sample \
-g ${outpath}hunt_sample2_fathers_4_chr${PBS_ARRAY_INDEX}.bgen \
-s ${outpath}hunt_sample2_fathers_4_chr${PBS_ARRAY_INDEX}_suffixes.sample \
-g ${outpath}hunt_sample2_fathers_5_chr${PBS_ARRAY_INDEX}.bgen \
-s ${outpath}hunt_sample2_fathers_5_chr${PBS_ARRAY_INDEX}_suffixes.sample \
-g ${outpath}hunt_sample2_fathers_6_chr${PBS_ARRAY_INDEX}.bgen \
-s ${outpath}hunt_sample2_fathers_6_chr${PBS_ARRAY_INDEX}_suffixes.sample \
-og ${outpath}hunt_sample2_fathers_chr${PBS_ARRAY_INDEX}.bgen \
-os ${outpath}hunt_sample2_fathers_chr${PBS_ARRAY_INDEX}.sample

# PLEASE CHECK THAT THE SAMPLE FILES PRODUCED IMMEDIATELY ABOVE HAVE THE CORRECT NUMBER OF IDS
# (I.E. IT SHOULD BE THE SAME NUMBER OF IDS AS IN dmo)

