# utility functions dealing with RCC files
#
# Author: samuelc
###############################################################################
require(reshape2)

# helper function to get value from matrix of key/value ------------------------
get.value.from.dictionary <- function(dictionary,key,to.character=TRUE) {
  return(
    ifelse(
      to.character,
      as.character(dictionary$value[dictionary$key==key]),
      dictionary$value[dictionary$key==key]
    )
  )
}


# parse a single rcc file for annotation ---------------------------------------
#
# filename: file name of RCC file
#
# return list of the following
#    sampleID
#    file.version
parse_rcc_annot_single_file <- function(filename, sample.name.format.function) {
  op <- options("stringsAsFactors")
  options(stringsAsFactors=FALSE)
  con <- file(filename,"r",blocking=FALSE)
  rcc.lines <- readLines(con)
  close(con)
  # header
  header.start.end <- which(rcc.lines %in% c("<Header>","</Header>"))
  if (length(header.start.end)!=2) {
    cat("ERROR encountered when parsing",filename," ... failed to find Header tag\n")	
  }
  header <- reshape2::colsplit(
      rcc.lines[(header.start.end[1]+1):(header.start.end[2]-1)],
      pattern=",",#split=",",
      names=c("key","value"))
  file.version <- get.value.from.dictionary(header,"FileVersion")
  
  # sample attributes
  sample.attributes.start.end <- which(rcc.lines %in% c("<Sample_Attributes>","</Sample_Attributes>"))
  if (length(sample.attributes.start.end)!=2) {
    cat("ERROR encountered when parsing",filename," ... failed to find Sample Attributes tag\n")	
  }
  sample.attributes <- reshape2::colsplit(
      rcc.lines[(sample.attributes.start.end[1]+1):(sample.attributes.start.end[2]-1)],
      pattern=",",#split=",",
      names=c("key","value"))
  sampleID            <- get.value.from.dictionary(sample.attributes,"ID")
  nanostring.operator <- get.value.from.dictionary(sample.attributes,"Owner")
  hyb.time            <- get.value.from.dictionary(sample.attributes,"Comments")
  nanostring.date     <- get.value.from.dictionary(sample.attributes,"Date")
  geneRLF             <- get.value.from.dictionary(sample.attributes,"GeneRLF")
  
  # lane attributes
  lane.attributes.start.end <- which(rcc.lines %in% c("<Lane_Attributes>","</Lane_Attributes>"))
  if (length(lane.attributes.start.end)!=2) {
    cat("ERROR encountered when parsing",filename," ... failed to find Sample Attributes tag\n")	
  }
  lane.attributes <- reshape2::colsplit(
      rcc.lines[(lane.attributes.start.end[1]+1):(lane.attributes.start.end[2]-1)],
      pattern=",",#split=",",
      names=c("key","value"))
  lane.number     <- get.value.from.dictionary(lane.attributes,"ID")
  fov.count       <- get.value.from.dictionary(lane.attributes,"FovCount")
  fov.counted     <- get.value.from.dictionary(lane.attributes,"FovCounted")
  scannerID       <- get.value.from.dictionary(lane.attributes,"ScannerID")
  stage.position  <- get.value.from.dictionary(lane.attributes,"StagePosition")
  binding.density <- get.value.from.dictionary(lane.attributes,"BindingDensity")
  cartridgeID     <- get.value.from.dictionary(lane.attributes,"CartridgeID")	
  
  # attributes that are derived from various attributes WITHIN the RCC file
  ottaID <- sampleID # what about non-otta study ??? 
  
  options(op)
  return(list(
          "File.Name.Original"=basename(filename), # original file name
          "File.Name"=sample.name.format.function(basename(filename)),
          "sampleID"=sampleID,
          "ottaID"=ottaID,
          "nanostring.site"=ifelse(scannerID=="DA46",NANOSTRING.SITE.VAN,ifelse(scannerID=="1303C0088",NANOSTRING.SITE.USC,NANOSTRING.SITE.AOC)),
          "nanostring.date"=nanostring.date,
          "nanostring.operator"=nanostring.operator,
          "scannerID"=scannerID,
          "file.version"=file.version,
          "geneRLF"=geneRLF,
          "cartridgeID"=cartridgeID,
          "lane.number"=lane.number,
          "hyb.time"=ifelse(
              hyb.time%in%c("long-hyb", "Long-Hyb"),
              HYB.TIME.LONG,
              ifelse(
                  hyb.time%in%c("short-hyb","Short-Hyb"),
                  HYB.TIME.SHORT,
                  paste("comment:",stringr::str_trim(hyb.time))
              )
          ),
          "stage.position"=stage.position,
          "fov.count"=fov.count,
          "fov.counted"=fov.counted,
          "binding.density"=binding.density
      ))
  
}

#' Parse and extract annotation from multiple rcc file(s) in a folder
#' 
#' Recursively go through all the input folder(s), searches for RCC files, parse
#' the annotation data and put in a single data frame.  The function will generate
#' warning messages if encounter files with different number of probes from the 
#' first RCC file it parsed.
#' 
#' @param older.name: to folder to look for RCC files, will look in subfolder as well; can be a list of folders
#' @param sample.name.format.function: function to format the file name into sample name
#' @param exclude: file(s) to exclude
#' @param get.excluded.files.instead: if set to true, get excluded files instead only
#' @return a data frame of the merged annotation data; row=probe, column=sample
#' @author Samuel Leung
#' @export
#' @examples
#' parse_rcc_expr(c("folder1","folder2"),function(x){x})
parse_rcc_annot <- function(folder.name, sample.name.format.function=NULL, exclude=NULL, get.excluded.files.instead=FALSE) {
  if (is.null(sample.name.format.function)) {
    # a function that does nothing
    sample.name.format.function <- function(x){x}
  }
  rcc.annot.d <- NULL
  folder.length.gt.1 <- length(folder.name)>1
  folder.name.is.not.dir <- ifelse(folder.length.gt.1,FALSE,!file.info(folder.name)$isdir)
  if (folder.name.is.not.dir) {
    # this is not a folder ... must be a file, try to parse it
    rcc.file.name <- basename(folder.name)
    if (get.excluded.files.instead == (rcc.file.name %in% exclude)) {
      if (endsWith(rcc.file.name,".rcc") | endsWith(rcc.file.name,".RCC")) {
        #cat("adding new sample:",sample.name.format.function(basename(folder.name)),"\n")
        rcc.annot.d <- parse_rcc_annot_single_file(folder.name, sample.name.format.function=sample.name.format.function)
        if (!is.null(rcc.annot.d)) {
          rcc.annot.d <- as.data.frame(t(unlist(rcc.annot.d)),stringsAsFactors=FALSE)
        }
      }
    }
  } else {
    if (length(folder.name)>1) {
      folders.to.iter <- folder.name			
    } else {
      folders.to.iter <- dir(folder.name, full.names=TRUE)
    }
    for (sub.folder.name in folders.to.iter) {
      temp <- parse_rcc_annot(sub.folder.name,sample.name.format.function,exclude=exclude, get.excluded.files.instead=get.excluded.files.instead)
      if (!is.null(temp)) {
        if (!is.null(rcc.annot.d)) {
          rcc.annot.d <- rbind(rcc.annot.d,temp)
        } else {
          rcc.annot.d <- temp # no need to rbind since this is the first parsed rcc file
        }
      }
    }	
  }
  return(rcc.annot.d)
}
