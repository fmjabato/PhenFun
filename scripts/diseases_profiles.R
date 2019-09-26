#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

##############################################################################
##                           CONFIGURE PROGRAM                              ##
##############################################################################

# Load necessary packages
suppressPackageStartupMessages(require(optparse))        # Parse script inputs

# Prepare input commands
option_list <- list(
  make_option(c("-f", "--file"), action="store", type="character",
              dest="file", help="File with Diseases HPO profiles (one entry per each Disease-HPO tuple)"),
  make_option(c("-i",	"--indexes"), action="store",type="character", default = "1,2,3,6,5",
              dest="indexes", help="String with Source, SourceID, DiseaseID, DiseaseName, and HPO column indexes, separated by commas [Default: %default]"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store clusters calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))

##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################

# Prepare indexes
# Obtain wanted indeexes
indx <- unlist(strsplit(opt$indexes,","))
indx <- list(Source   = as.integer(indx[1]),
             SourceID = as.integer(indx[2]),
             DisName  = as.integer(indx[3]),
             DisID    = as.integer(indx[4]),
             HPO      = as.integer(indx[5]))

# Load file
diseases <- read.table(opt$file, sep = "\t", header = F, quote = "", stringsAsFactors = F, comment.char = "")


##############################################################################
##                             TRANSFORM DATA                               ##
##############################################################################

# Group by tuples <Source,DiseaseName>
invisible(lapply(unique(diseases[,indx$Source]),function(source){
  indx_s <- which(diseases[,indx$Source] == source)
  # Study each disease
  invisible(lapply(unique(diseases[indx_s,indx$SourceID]),function(dis){
    # Obtain disease indexes
    indx_d <- intersect(indx_s,which(diseases[,indx$SourceID] == dis))
    # Update
    diseases[indx_d,indx$DisID] <<- paste(source,dis,sep=":")
  }))
}))

##############################################################################
##                              WRITE OUTPUT                                ##
##############################################################################

write.table(diseases[,c(indx$DisID,indx$DisName,indx$HPO)], file = opt$output, sep = "\t", quote = F, col.names = F, row.names = F)