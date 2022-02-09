#!bin/bash

# create .bgi index files

outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # path to bgen files
bgenix=/path/to/bgenix/executable # path to bgenix executable

for i in {1..22}; do
  ${bgenix} -g "${outpath}"chr${i}.bgen -index
done
