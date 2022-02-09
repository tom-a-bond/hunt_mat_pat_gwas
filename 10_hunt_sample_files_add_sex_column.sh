#!/bin/bash

# gcta expects the sample files to have a 4th column ("sex"), so we need to add a dummy column to all sample files
# n.b. I haven't been able to test this script, so it may require some debugging

script_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/scripts/ # path to directory containing hunt_params.sh
cohort_name=hunt

# samples 1 and 2

outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # folder containing final bgens written by subset_hunt_bgen_sample1.sh and subset_hunt_bgen_sample2.sh
for i in {1..22}; do
  # sample 1
  mv ${outpath}hunt_sample1_mothers_chr${i}.sample ${outpath}hunt_sample1_mothers_chr${i}.sample.nosexcol
  awk '{print $1, $2, $3, $3}' ${outpath}hunt_sample1_mothers_chr${i}.sample.nosexcol | \
  awk 'NR == 1 {$4 = "sex"}1' | \
  awk 'NR == 2 {$4 = "D"}1' > \
  ${outpath}hunt_sample1_mothers_chr${i}.sample
  # sample 2
  mv ${outpath}hunt_sample2_fathers_chr${i}.sample ${outpath}hunt_sample2_fathers_chr${i}.sample.nosexcol
  awk '{print $1, $2, $3, $3}' ${outpath}hunt_sample2_fathers_chr${i}.sample.nosexcol | \
  awk 'NR == 1 {$4 = "sex"}1' | \
  awk 'NR == 2 {$4 = "D"}1' > \
  ${outpath}hunt_sample2_fathers_chr${i}.sample
done

# sample 3

sample=3
source ${script_path}${cohort_name}_params.sh
for i in {1..22}; do
  sample_file=chr${i}.sample # N.B. THIS MAY NEED TO BE EDITED
  mv ${sample_path}${sample_file} ${sample_path}${sample_file}.nosexcol
  awk '{print $1, $2, $3, $3}' ${sample_path}${sample_file}.nosexcol | \
  awk 'NR == 1 {$4 = "sex"}1' | \
  awk 'NR == 2 {$4 = "D"}1' > \
  ${sample_path}${sample_file}
done

# n.b. samples 4m and 4f probably do not require editing now we are only analysing unrelateds, because plink will probably read them (i.e. the sample files written in 9_hunt_sample4_subset_unrelated.sh) fine

# sample 4m

# sample=4m
# source ${script_path}${cohort_name}_params.sh
# mv ${sample_path}${sample_file} ${sample_path}${sample_file}.nosexcol
# awk '{print $1, $2, $3, $3}' ${sample_path}${sample_file}.nosexcol | \
# awk 'NR == 1 {$4 = "sex"}1' | \
# awk 'NR == 2 {$4 = "D"}1' > \
# ${sample_path}${sample_file}

# sample 4f

# sample=4f
# source ${script_path}${cohort_name}_params.sh
# mv ${sample_path}${sample_file} ${sample_path}${sample_file}.nosexcol
# awk '{print $1, $2, $3, $3}' ${sample_path}${sample_file}.nosexcol | \
# awk 'NR == 1 {$4 = "sex"}1' | \
# awk 'NR == 2 {$4 = "D"}1' > \
# ${sample_path}${sample_file}

# sample 4o (used for samples 4m and 4f, but probably does not require editing, assuing plink is not as fussy as gcta)

# sample=4f
# source ${script_path}${cohort_name}_params.sh
# mv ${sample4_sample_path}${sample4o_sample_file} ${sample4_sample_path}${sample4o_sample_file}.nosexcol
# awk '{print $1, $2, $3, $3}' ${sample4_sample_path}${sample4o_sample_file}.nosexcol | \
# awk 'NR == 1 {$4 = "sex"}1' | \
# awk 'NR == 2 {$4 = "D"}1' > \
# ${sample4_sample_path}${sample4o_sample_file}




