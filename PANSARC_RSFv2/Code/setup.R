# load libraries and setup data
library(rmarkdown)
library(pansarc)
library(rChoiceDialogs)

# ask for input folder
input_folder_name <- jchoose.dir(
  default="D:\\Users\\INSTR-ADMIN\\Desktop\\PANSARC_RSFv2\\RCC data",
  caption="Please select folder containing RCC files")
if (is.na(input_folder_name)) {
  stop("No RCC folder selected")
}

# ask for output folder
output_folder_name <- choose.dir(
  default="D:\\Users\\INSTR-ADMIN\\Desktop\\PANSARC_RSFv2\\Output", 
  caption="Please select folder to store output files")

# do the analysis!
result <- analyze(
  input_folder_name,
  c("Nanostring probe list with medians.xlsx","probe medians"),
  output_folder_name,
  pos_ctl_bound=c(0.3,3),
  sample_ratio_threshold=5,
  hk_qc_threshold=0.3,
  do_main=TRUE
)