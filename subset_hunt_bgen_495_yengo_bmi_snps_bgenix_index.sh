#!bin/bash

# create .bgi index files for duos/trios yengo snps bgen files
# this script should be quick to run
# n.b. bgenix must be installed prior to running this script; see section "Obtaining and installing BGEN" at https://enkre.net/cgi-bin/code/bgen/dir?ci=tip

outpath=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/data/ukb/ # path to bgen files for duos/trios (created by subset_hunt_bgen_495_yengo_bmi_snps.sh)
n_yengo_snps=495 # 495 yengo snps should be available in all cohorts

bgenix -g "${outpath}"hunt_sample4_fathers_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen -index
bgenix -g "${outpath}"hunt_sample4_mothers_yengo_${n_yengo_snps}_bmi_snps_all_chr.bgen -index
