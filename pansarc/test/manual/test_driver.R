# test driver
setwd("/home/samuelc/Documents/workspace/R/pansarc/PANSARC_RSFv2")

source("../pansarc/R/parse_rcc_annot.R")
source("../pansarc/R/parse_rcc_expr.R")
source("../pansarc/R/pos_ctl_norm_factors.R")
source("../pansarc/R/analyze.R")

input_fname <- "/mnt/vm_shared/mapcore/pansarc testing/20191213_julie/rcc/20191114_20191113 PanSarcoma Validation-2_RCC"
probe_med_finfo <- c("Nanostring probe list with medians.xlsx","probe medians")
output_dir <- "/home/samuelc/Desktop"
pos_ctl_bound <- c(0.3,3)
sample_ratio_threshold <- "Nanostring custom threshold.xlsx"
hk_qc_threshold <- 0.3
do_main <- TRUE

result <- analyze(
  input_fname,
  probe_med_finfo,
  output_dir,
  pos_ctl_bound,
  sample_ratio_threshold,
  hk_qc_threshold,
  do_main
)
