# Pan Sarcoma Fusion Gene Nanostring Assay - Report of Statistical Findings

# clean up
rm(list = ls(all = TRUE))  # Clean up

source("Code/setup.R")
SUCCEED <- length(result)>2 # length(result)==2 indicate some erros

output_file <- paste0("PANSARC_RSF_",format(Sys.Date(),format="%Y-%m-%d"),".pdf")
output_dir <- getwd()

if (SUCCEED) {
  rmarkdown::render(
    input="comp/succeed.Rmd",
    output_format="pdf_document",
    output_file=output_file,
    output_dir=output_dir)
} else {
  rmarkdown::render(
    input="comp/failed.Rmd",
    output_format="pdf_document",
    output_file=output_file,
    output_dir=output_dir)
}

