# TODO: Add comment
# 
# Author: samle
###############################################################################
require(assertthat)

# calculate slopes between two points ------------------------------------------
#
# assume input_d is a 2-col data frame or matrix with number of 
# rows = number of points
#
# return NA if there is only one point
slopes_btw_pts <- function(input_d) {
  assertthat::assert_that(ncol(input_d)==2)
  num_rows <- nrow(input_d)
  if (num_rows < 2) {
    # cannot calculate slot for only one point
    return(NA)
  } 
  for (i in ncol(input_d)) {
    input_d[,i] <- as.numeric(input_d[,i])
  }
  input_d <- input_d[order(input_d[,1]),]
  return(sapply(2:nrow(input_d),function(x){
            y <- input_d[x,]-input_d[x-1,]
            y[2]/y[1]
          }))
}

#' Calculate normalization factor for positive control
#' 
#' calculate normalization factor for positive control
#' 
#' @param input_d: data frame of nanostring expression
#' @param norm_func: function for calculating positive control normalization factor; usually either sum, mean or geometric mean
#' @param rcc_codeclass_col_name: column name for the 'CodeClass'
#' @param rcc_name_col_name: column name for 'Name'
#' @param rcc_non_sample_col_names: names of column that do not corresponds to expression data e.g. Name, CodeClass
#' @return a list with the following elements
#' 	"values"=positive control normalization factor values
#'  "norm_d"=data frame of normalized expression data for endogenous probes
#'  "hk_d"=data frame of normalized expression data for housekeepers
#'  "pctl_d"=data frame of normalized expression data for positive controls
#'  "nctl_d"=data frame of normalized expression data for negative controls
#' @author Samuel Leung
#' @export
#' @examples
#' pos_ctl_norm_factors(r_rcc,norm_func=function(a){prod(a)^(1/length(a))})
pos_ctl_norm_factors <- function(
    input_d, 
    norm_func=sum,
    rcc_codeclass_col_name="CodeClass", 
    rcc_name_col_name="Name", 
    rcc_non_sample_col_names=c("CodeClass", "Name", "Accession")) {
  pos_ctl_row_indexes <- which(input_d[,rcc_codeclass_col_name]=="Positive")
  neg_ctl_row_indexes <- which(input_d[,rcc_codeclass_col_name]=="Negative") 
  hk_row_indexes      <- which(input_d[,rcc_codeclass_col_name]=="Housekeeping") 
  
  sample_col_indexes <- which(!names(input_d)%in%rcc_non_sample_col_names)
  sum_pos_ctls <- apply(input_d[pos_ctl_row_indexes,sample_col_indexes],2,norm_func)
  avg_sum_pos_ctls <- mean(sum_pos_ctls)
  pos_ctl_norm_factors <- avg_sum_pos_ctls/sum_pos_ctls
  
  # QC: make sure the slop is more or less linear
  # - currently just make sure no slope between any two point is < 0
  pos_ctl_values <- sapply(input_d$Name[which(input_d$CodeClass=="Positive")],function(x){as.numeric(strsplit(x,"\\(|\\)")[[1]][2])})
  assertthat::assert_that(min(
    sapply(sample_col_indexes,function(x){
      result <- slopes_btw_pts(cbind(pos_ctl_values,input_d[pos_ctl_row_indexes,x]))
      if (min(result)<=0) {
        stop(paste("positive controls failed:",names(input_d)[x]))
      }
      return(result)
    }))>0)
  
  # normalize with positive control
  output_d <- input_d
  output_d[,sample_col_indexes] <- input_d[,sample_col_indexes] * matrix(rep(pos_ctl_norm_factors,times=nrow(input_d)),ncol=length(sample_col_indexes),byrow=TRUE)
  # remove housekeepers pos/neg controls
  norm_d <- output_d[-c(pos_ctl_row_indexes,neg_ctl_row_indexes,hk_row_indexes),]
  rownames(norm_d) <- norm_d[,rcc_name_col_name]
  norm_d <- norm_d[,sample_col_indexes]
  # housekeepers expression matrix
  hk_d <- output_d[hk_row_indexes,]
  rownames(hk_d) <- hk_d[,rcc_name_col_name]
  hk_d <- hk_d[,sample_col_indexes]
  # positive control expression matrix
  pctl_d <- output_d[pos_ctl_row_indexes,]
  rownames(pctl_d) <- pctl_d[,rcc_name_col_name]
  pctl_d <- pctl_d[,sample_col_indexes]
  # negative control expression matrix
  nctl_d <- output_d[neg_ctl_row_indexes,]
  rownames(nctl_d) <- nctl_d[,rcc_name_col_name]
  nctl_d <- nctl_d[,sample_col_indexes]
  return(list(
          "values"=pos_ctl_norm_factors,
          "norm_d"=norm_d,
          "hk_d"=hk_d,
          "pctl_d"=pctl_d,
          "nctl_d"=nctl_d)
  )
}



