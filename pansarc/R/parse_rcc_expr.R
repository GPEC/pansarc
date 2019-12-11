# utility functions dealing with RCC files
#
# Author: samuelc
###############################################################################
require(reshape2)

#############################
### some constants        ###
RCC.CODE.SUMMARY.COLNAME.CODECLASS <- "CodeClass"
RCC.CODE.SUMMARY.COLNAME.NAME      <- "Name"
RCC.CODE.SUMMARY.COLNAME.ACCESSION <- "Accession"
RCC.CODE.SUMMARY.COLNAME.COUNT     <- "Count"
NANOSTRING.SITE.VAN <- "Vancouver"
NANOSTRING.SITE.USC <- "USC"
NANOSTRING.SITE.AOC <- "AOC"
HYB.TIME.LONG <- "Long-Hyb"
HYB.TIME.SHORT <- "Short-Hyb"
### end of some constants ###
#############################

# parse a single rcc file for expression data ----------------------------------
#
# filename: file name of RCC file
#
# return a data matrix with the following columns:
#        "Name","Accession","Code.Class", [column.name]
parse_rcc_expr_single_file <- function(filename, column.name="Count") {
	op <- options("stringsAsFactors")
	options(stringsAsFactors=FALSE)
	con <- file(filename,"r",blocking=FALSE)
	rcc.lines <- readLines(con)
	close(con)
	cs.start.end <- which(rcc.lines %in% c("<Code_Summary>","</Code_Summary>"))
	if (length(cs.start.end)!=2) {
		cat("ERROR encountered when parsing",filename," ... failed to find Code summary tag\n")	
	}
	rcc.expr.header <- strsplit(rcc.lines[cs.start.end[1]+1],",")[[1]]
	rcc.expr.d      <- reshape2::colsplit(
        rcc.lines[c((cs.start.end[1]+2):(cs.start.end[2]-1))],
        pattern=",",#split=",",
        names=rcc.expr.header)
	# don't know how to set colsplit to NOT generate factor from strings ... therefore, need to do manually
	rcc.expr.d <- cbind(
		as.data.frame(apply(rcc.expr.d[,names(rcc.expr.d)!=RCC.CODE.SUMMARY.COLNAME.COUNT],2,as.character)),
		"TEMP_VAR_NAME_Count"=rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.COUNT]
	)
	names(rcc.expr.d)[names(rcc.expr.d)=="TEMP_VAR_NAME_Count"] <- column.name
	options(op)
	return(rcc.expr.d)
}

#' Parse and extract expression data from multiple rcc file(s) in a folder
#' 
#' Recursively go through all the input folder(s), searches for RCC files, parse
#' the expression data and put in a single data frame.  The function will generate
#' warning messages if encounter files with different number of probes from the 
#' first RCC file it parsed.
#' 
#' @param older.name: to folder to look for RCC files, will look in subfolder as well; can be a list of folders
#' @param sample.name.format.function: function to format the file name into sample name
#' @param exclude: file(s) to exclude
#' @param get.excluded.files.instead: if set to true, get excluded files instead only
#' @return a data frame of the merged expression data; row=probe, column=sample
#' @author Samuel Leung
#' @export
#' @examples
#' parse_rcc_expr(c("folder1","folder2"),function(x){x})
parse_rcc_expr <- function(folder.name, sample.name.format.function=NULL, exclude=NULL, get.excluded.files.instead=FALSE) {
  if (is.null(sample.name.format.function)) {
    # a function that does nothing
    sample.name.format.function <- function(x){x}
  }
	rcc.expr.d <- NULL
	folder.length.gt.1 <- length(folder.name)>1
	folder.name.is.not.dir <- ifelse(folder.length.gt.1,FALSE,!file.info(folder.name)$isdir)
	if (folder.name.is.not.dir) {
		# this is not a folder ... must be a file, try to parse it
		rcc.file.name <- basename(folder.name)
		if (get.excluded.files.instead == (rcc.file.name %in% exclude)) {
			if (endsWith(rcc.file.name,".rcc") | endsWith(rcc.file.name,".RCC")) {
				#cat("adding new sample:",sample.name.format.function(basename(folder.name)),"\n")
				rcc.expr.d <- parse_rcc_expr_single_file(folder.name, column.name=sample.name.format.function(rcc.file.name))
			}
		}
	} else {
		if (length(folder.name)>1) {
			folders.to.iter <- folder.name			
		} else {
			folders.to.iter <- dir(folder.name, full.names=TRUE)
		}
		for (sub.folder.name in folders.to.iter) {
			temp <- parse_rcc_expr(sub.folder.name,sample.name.format.function,exclude=exclude,get.excluded.files.instead=get.excluded.files.instead)
			if (!is.null(temp)) {
				if (!is.null(rcc.expr.d)) {
					# we require ALL genes to be the same ... i.e. cannot merge RCC files with different number of genes
					check.1 <- match(
						paste(rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.NAME],rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.ACCESSION],sep="___"),
						paste(      temp[,RCC.CODE.SUMMARY.COLNAME.NAME],      temp[,RCC.CODE.SUMMARY.COLNAME.ACCESSION],sep="___")				
					)
					check.2 <- match(
						paste(      temp[,RCC.CODE.SUMMARY.COLNAME.NAME],      temp[,RCC.CODE.SUMMARY.COLNAME.ACCESSION],sep="___"),	
						paste(rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.NAME],rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.ACCESSION],sep="___")			
					)
					if (sum(is.na(check.1))>0 | sum(is.na(check.2))) {
						warning(paste0("(merging RCCs): different set of probes ... SKIPPING ",(ncol(temp)-3)," RCCs in:",sub.folder.name))
					} else {			
						rcc.expr.d <- cbind(
							rcc.expr.d,
							temp[match(
								paste(rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.NAME],rcc.expr.d[,RCC.CODE.SUMMARY.COLNAME.ACCESSION],sep="___"),
								paste(      temp[,RCC.CODE.SUMMARY.COLNAME.NAME],      temp[,RCC.CODE.SUMMARY.COLNAME.ACCESSION],sep="___")				
							),],
							stringsAsFactors=FALSE
						)
						# assume first 3 column is ALWAYS CodeClass/Accession/Name ... this is guarantee by parse_rcc_expr_single_file
						rcc.expr.d <- rcc.expr.d[,-(c(1,2,3)+(ncol(rcc.expr.d)-ncol(temp)))]
					}
				} else {t
					rcc.expr.d <- temp # no need to cbind since this is the first parsed rcc file
				}
			}
		}	
	}
	return(rcc.expr.d)
}