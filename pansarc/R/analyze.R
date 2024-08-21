# function to do the main nanostring data analysis
# 
# Author: samle
###############################################################################
require(xlsx)
require(assertthat)

#' Analyze nanostring expression data to detect fusion probe
#' 
#' Analyze nanostring expression data to detect fusion probe
#' 
#' @param input_fname: input file/folder name(s)
#' @param probe_med_finfo: info regarding probe median of initial cases; assume 
#'  this an array of 2 elements: excel file name, excel worksheet name
#' @param output_dir: output folder name for data export, if NA, no data will be exported
#' @param pos_ctl_bound: bound for positive control normalization data check; default 0.3-3
#' @param sample_ratio_threshold: threshold for fusion probe detection; if a file is provided, this file would contain custom threshold for individual genes; default: 5
#' @param hk_qc_threshold: housekeeper QC threshold; default 0.3
#' @param do_main: if TRUE do analysis with main codeset
#' @return a list with the following elements
#' 	n_obj=processed expression matrixes;    
#' 	p_ratios=probe ratios; 
#' 	s_ratios=sample ratio; 
#' 	hk_ratios=housekeeper ratios; 
#' 	do_report=function to do report (e.g. barplot) on individual sample;  
#' @author Samuel Leung
#' @export
#' @examples
#' pos_ctl_norm_factors(r_rcc,norm_func=function(a){prod(a)^(1/length(a))})
analyze <- function(input_fname,
    probe_med_finfo,
    output_dir,
    pos_ctl_bound=c(0.3,3),
    sample_ratio_threshold=NA,
    hk_qc_threshold=0.3,
    do_main=TRUE) {
  if (!do_main) {
    stop("currently only support main codeset.")
  }
  # force load package
  library(xlsx)

  # 1. parse RCC files go get raw counts -----------------------------------------
  if (do_main) {
    r_rcc <- parse_rcc_expr(input_fname,function(x){x})
  }
  
  # 2. positive control normalization --------------------------------------------
  n_obj <- pos_ctl_norm_factors(r_rcc,norm_func=function(a){prod(a)^(1/length(a))}) 
  n_rcc <- n_obj$norm_d
  n_rcc_hk <- n_obj$hk_d
  n_rcc_pctl <- n_obj$pctl_d
  n_rcc_nctl <- n_obj$nctl_d
  # make sure positive control normalization factor is between 0.3-3
  for (x in names(n_obj$values)) {
    if (n_obj$values[x]<=pos_ctl_bound[1] | n_obj$values[x]>=pos_ctl_bound[2]) {
      warning(paste0("Positive control normalization factor <=",pos_ctl_bound[1]," or >=",pos_ctl_bound[2],": ",x," (",n_obj$values[x],")"))
    }
  }
  if(min(n_obj$values)<=pos_ctl_bound[1] | max(n_obj$values)>=pos_ctl_bound[2]) {
    # return error message
    return(list(
      "ERR_MSG"="Positive control normalization factor failed.",
      "n_obj"=n_obj
    ))
  }
  
  # 3. calculate probe ratios = normalized count / median probe expression -------
  # read probe median file
  probe_med <- xlsx::read.xlsx2(probe_med_finfo[1],sheetName=probe_med_finfo[2],stringsAsFactors=FALSE)
  if (do_main) {
    probe_med_col_name <- "Main.CodeSet"
  }
  cat("NOTE: parsing probe median file (please make sure info are correct) ...\n")
  assertthat::assert_that(probe_med_col_name%in%names(probe_med))
  probe_med_col_value <- names(probe_med)[which(names(probe_med)==probe_med_col_name)+1]
  cat("  looking for column:",paste0("'",probe_med_col_name,"'"),"assuming values in column:",paste0("'",probe_med_col_value ,"'"),"\n")
  probe_med_blank_row_indexes <- sort(which(is.na(probe_med[,probe_med_col_name])|stringr::str_trim(probe_med[,probe_med_col_name])==""))
  probe_med_nhk <- probe_med[1:(probe_med_blank_row_indexes[1]-1),c(probe_med_col_name,probe_med_col_value)]
  probe_med_hk  <- probe_med[(grep("Control Genes - average normalized counts",probe_med[,probe_med_col_name])+1):(probe_med_blank_row_indexes[2]-1),c(probe_med_col_name,probe_med_col_value)]
  cat("  number of non-house keeper genes (exclude pos/neg ctls):",nrow(probe_med_nhk),"\n")
  cat("  number of house keeper genes:",nrow(probe_med_hk),"\n")
  # make sure all probe_med genes are in rcc file
  assertthat::assert_that(sum(probe_med_nhk[,probe_med_col_name]%in%rownames(n_rcc))==nrow(probe_med_nhk))
  # make sure probe_med genes are unique - withwise matching will fail!
  assertthat::assert_that(length(unique(probe_med_nhk[,probe_med_col_name]))==length(probe_med_nhk[,probe_med_col_name]))
  assertthat::assert_that(length(unique(probe_med_hk[, probe_med_col_name]))==length(probe_med_hk[, probe_med_col_name]))
  
  # calculate probe ratios
  p_ratios <- n_rcc / matrix(as.numeric(rep(probe_med_nhk[match(rownames(n_rcc),probe_med_nhk[,probe_med_col_name]),probe_med_col_value],times=ncol(n_rcc))),ncol=ncol(n_rcc),byrow=FALSE)
  
  # calculate PRR = PR / median (all probe_ratios within sample except housekeeper)
  s_ratios = p_ratios / matrix(rep(apply(p_ratios,2,median),times=nrow(n_rcc)),ncol=ncol(n_rcc),byrow=TRUE)
  
  # figure out threshold
  DEFAULT_SAMPLE_RATIO_THRESHOLD <- 5
  if (is.na(sample_ratio_threshold)) {
    # nothing is specified for threshold, set threshold to default threshold
    sample_ratio_threshold <- DEFAULT_SAMPLE_RATIO_THRESHOLD
  }
  # check if sample_ratio_threshold is in fact a filename
  sr_th_default <- NA
  if (file.exists(as.character(sample_ratio_threshold))) {
    sr_th_default <- DEFAULT_SAMPLE_RATIO_THRESHOLD
    sr_th_custom_d <- xlsx::read.xlsx2(sample_ratio_threshold,sheetIndex=1)
    assertthat::assert_that(
      sum(sr_th_custom_d$gene %in% rownames(s_ratios))==nrow(sr_th_custom_d),
      msg=paste0("unknown genes in '",sample_ratio_threshold,"': ",
        paste(sr_th_custom_d$gene[!sr_th_custom_d$gene %in% rownames(s_ratios)], collapse = ", ")
      )
    )
  } else {
    # assume this corresponds to a single threshold (number)
    sr_th_default <- tryCatch(
      as.numeric(sample_ratio_threshold),
      warning=function(cond) {
        stop(paste0("incorrect input to 'sample_ratio_threshold' ... this needs to be either a number or an Excel file."))
      }
    )
    sr_th_custom_d <- data.frame(
      gene=character(0),
      threshold=integer(0)
    )
  }

  # 4. QC by housekeepers --------------------------------------------------------
  #    a. for each HK, calculate normalized HK count / median HK count for initial 36 samples
  #    b. for each sample, calculate the average of HK values in a)
  p_ratios_hk <- n_rcc_hk / matrix(as.numeric(rep(probe_med_hk[match(rownames(n_rcc_hk),probe_med_hk[,probe_med_col_name]),probe_med_col_value],times=ncol(n_rcc_hk))),ncol=ncol(n_rcc_hk),byrow=FALSE)
  hk_ratios <- apply(p_ratios_hk,2,mean)
  hk_failed <- hk_ratios[which(hk_ratios<hk_qc_threshold)]
  
  # 5. generate plots and indicate fusion genes ----------------------------------
  do_report <- function(sample_id) {
    sample_ratio_threshold_arr <- rep(sr_th_default,nrow(s_ratios))
    if (nrow(sr_th_custom_d)>0) {
      sample_ratio_threshold_arr[match(sr_th_custom_d$gene,rownames(s_ratios))] <- sr_th_custom_d$threshold
      sample_ratio_threshold_arr <- as.numeric(sample_ratio_threshold_arr)
    }

    fg_indexes <- which(s_ratios[,sample_id]>sample_ratio_threshold_arr)
    fg_detected <- length(fg_indexes)>0
    
    # function to generate a table of detected fusion probe and sample ratio
    get_table <- function() {
      data.frame(
          "probe"=rownames(s_ratios)[fg_indexes],
          "sample ratio"=round(s_ratios[fg_indexes,sample_id],digits=1)
      )
    }
    
    # function to generate barplot
    get_plot <- function() {
      barplot(
          s_ratios[,sample_id],
          ylab="Sample ratio",
          xlab="Nanostring probe",
          ylim=c(0,max((sr_th_default*1.5),max(s_ratios[,sample_id]))),
          main=paste0(sample_id)
      )
      abline(h=sr_th_default)
      if (fg_detected) {
        result_text <- paste0(length(fg_indexes)," fusion probe",ifelse(length(fg_indexes)>1,"s","")," detected")
        legend("topright",title=result_text,apply(get_table(),1,paste,collapse=": "))
      } else {
        mtext(paste0("no fusion probe detected (","HK ratio:",round(hk_ratios[sample_id],digits=2),")"))
      }
    }    
    
    return(list(
            "fg_detected"=length(fg_indexes)>0,
            "get_plot"=get_plot,
            "get_table"=get_table
        ))
  }
  
  # 6. output data ---------------------------------------------------------------
  if (!is.na(output_dir)) {
    out_fname <- file.path(output_dir,paste0("pan_sarc_nstring_",format(Sys.Date(),format="%Y-%m-%d"),".xlsx"))
    cat("# writing output data to",out_fname,"...")
    out_d <- xlsx::write.xlsx2(r_rcc,file=out_fname,   sheetName="raw count",                   col.names=TRUE, row.name=FALSE)
    out_d <- xlsx::write.xlsx2(n_rcc,file=out_fname,   sheetName="normalized endogenous count", col.names=TRUE, row.name=TRUE, append=TRUE)
    out_d <- xlsx::write.xlsx2(n_rcc_hk,file=out_fname,sheetName="normalized housekeeper count",col.names=TRUE, row.name=TRUE, append=TRUE)
    out_d <- xlsx::write.xlsx2(
        cbind(
            "CodeClass"=c(rep("endogenous",times=nrow(probe_med_nhk)),rep("housekeeper",times=nrow(probe_med_hk))),
            rbind(probe_med_nhk,probe_med_hk)
        ),
        file=out_fname,sheetName="probe median (initial samples)",col.names=TRUE, row.name=FALSE, append=TRUE)
    out_d <- xlsx::write.xlsx2(p_ratios,file=out_fname,sheetName="probe ratio", col.names=TRUE, row.name=TRUE, append=TRUE)
    out_d <- xlsx::write.xlsx2(s_ratios,file=out_fname,sheetName="sample ratio",col.names=TRUE, row.name=TRUE, append=TRUE)
  }
  
  cat("done.  bye.\n")
  cat(date(),"\n")
  
  return(list(
          "n_obj"=n_obj,   
          "p_ratios"=p_ratios,
          "p_ratios_hk"=p_ratios_hk,
          "s_ratios"=s_ratios,
          "hk_ratios"=hk_ratios,
          "do_report"=do_report        
      ))
}

