#!bin/bash

# create .bgi index files for sample 1/2 subsetted bgen files
# I am not sure how long this script will take to run; it may need to be submitted as a script rather than run interactively
# n.b. bgenix must be installed prior to running this script; see section "Obtaining and installing BGEN" at https://enkre.net/cgi-bin/code/bgen/dir?ci=tip

outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # path to bgen files for sample 1/2 ("outpath" in subset_hunt_bgen_sample1.sh and subset_hunt_bgen_sample2.sh)

for i in {1..22}; do
  bgenix -g ${outpath}hunt_sample1_mothers_chr${i}.bgen -index
  bgenix -g ${outpath}hunt_sample2_fathers_chr${i}.bgen -index
done
