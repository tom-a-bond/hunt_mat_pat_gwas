args = commandArgs(trailingOnly = TRUE)
qc_config_file = args[1]
library(GWASinspector)
# read in params
job <- setup.inspector(qc_config_file)
print(job)
# run qc pipeline
run.inspector(job)