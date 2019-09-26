#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
require(optparse)

option_list <- list(
  make_option(c("-i", "--input"), action="store", type="character",
              dest="input", help="Relationships file"),
  make_option(c("-t", "--target"), action="store", type="integer", default = 1,
              dest="target", help="Column where targets will be found. [Default: %default]"),
  make_option(c("-a","--all"), action="store_true", default = FALSE,
              dest="all", help="Flag used to indicate that maximum frequence must be added as the number of patients and linked to HP:0000001 (All)"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store ratios calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
 # Load filtered signal
network <- read.table(file = opt$input, header = F, sep = "\t", stringsAsFactors = F)
targets <- network[which(grepl("HP:",network[,opt$target])),opt$target]
targets <- table(targets)

# Transform table to dataframe
dftargets <- data.frame(HPO = names(targets), Freq = as.vector(targets), stringsAsFactors = F)

# Add patients frequencies
if(opt$all){
	# Calc max freq
	pats_freq <- length(unique(network[-which(grepl("HP:",network[,opt$target])),opt$target]))
	# Add
	if("HP:0000001" %in% dftargets$HPO){
		dftargets$Freq[which(dftargets$HPO == "HP:0000001")] <- pats_freq
	}else{
		dftargets <- rbind(data.frame(HPO = "HP:0000001", Freq = pats_freq, stringsAsFactors = F), dftargets)
	}
}

################################################################
##                          WRITE DATA                        ##
################################################################
write.table(dftargets, file = opt$output, col.names = F, row.names = F, sep = "\t", quote = F)