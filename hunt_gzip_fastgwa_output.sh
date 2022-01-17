#!bin/bash

# add fastGWA output and log files to gzipped tarball
# n.b. I am unsure of runtime, so this may need to be queued rather than run interactively

out_path=/dmf/uqdi/HPC/PBSHOME/ttbond/proj/mat_pat_bmi_mr/results/gwas/duos_trios/ # output path for gwas results
cohort_name=hunt # cohort name; one of "ukb", "hunt", "alspac"
mkdir ${out_path}${cohort_name}_to_send
mv ${out_path}*.fastGWA.gz ${out_path}${cohort_name}_to_send # move output files to a new directory named "to_send"
mv ${out_path}*.log ${out_path}${cohort_name}_to_send # move output files to a new directory named "to_send"
cd ${out_path}
tar -C ${out_path} -cvzf ${out_path}${cohort_name}_to_send.tar.gz ./${cohort_name}_to_send

# N.B. THERE SHOULD BE 2780 FILES IN ${out_path}${cohort_name}_to_send AS SHOWN BY THE FOLLOWING R CODE; IF THERE ARE NOT, PLEASE CONTACT TOM SO WE CAN DIAGNOSE WHY
# n_phens = 20
# n_chrs = 22
# n_samples_gwas_fastgwa = 3 # samples 1, 2, 3
# n_samples_bmi_snps_fastgwa = 2 # samples 4m, 4f
# n_phens * n_chrs * n_samples_gwas_fastgwa * 2 + n_phens * n_samples_gwas_fastgwa + n_phens * n_samples_bmi_snps_fastgwa * 2


