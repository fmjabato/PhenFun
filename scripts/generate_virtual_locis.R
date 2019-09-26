#! /usr/bin/env Rscript

#' @description randomizate locis position into the genome
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

##############################################################################
##                           CONFIGURE PROGRAM                              ##
##############################################################################

# Load necessary packages
suppressPackageStartupMessages(require(optparse))        # Parse script inputs

# Prepare input commands
option_list <- list(
  make_option(c("-l", "--locis"), action="store", type="character",
              dest="locis", help="Network with locis. Need 4 columns: Chr,Start,End,NodeID"),
  make_option(c("-g", "--genome"), action="store", type="character",
              dest="genome", help="Genome's chromosomes file. A two columns file with Chr_ID nad Chr_Size"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Destiny file")
)

opt <- parse_args(OptionParser(option_list=option_list))

##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################
# Load locis
locis <- read.table(file = opt$locis, sep = "\t", header = F, quote = "",stringsAsFactors = F)
colnames(locis) <- c("Chr","Start","End","Term","Value","NodeID")
# Load sizes
chr_sizes <- read.table(file = opt$genome, sep = "\t", header = F, quote = "",stringsAsFactors = F)
colnames(chr_sizes) <- c("Chr","Size")

##############################################################################
##                                RANDOMIZE                                 ##
##############################################################################

# Obtain unique list of NodeIDs
nodes <- unique(locis$NodeID)

# Per each node, generate a virtual loci of same length
vlocis <- as.data.frame(do.call(rbind,lapply(nodes,function(id){
	# Find instance of it loci
	indexes <- which(locis$NodeID == id)
	# Obtain loci size
	loci_size <- locis$End[indexes[1]] - locis$Start[indexes[1]]
	# Take a random chromosome
	correctNotFound = TRUE
	while(correctNotFound){
		rdm <- sample(1:nrow(chr_sizes),1)
		if(chr_sizes$Size[rdm] > loci_size){
			correctNotFound <- FALSE
		}
	}
	# Take a random position
	pos <- sample(1:(chr_sizes$Size[rdm] - loci_size),1)
	# Generate random loci
	return(data.frame(NodeID = id, 
		              Chr = chr_sizes$Chr[rdm],
		              Start = pos,
		              End = pos + loci_size,
	                  stringsAsFactors = F))
})))

# Substitute locis
invisible(lapply(seq(nrow(vlocis)),function(i){
	indexes <- which(locis$NodeID == vlocis$NodeID[i])
	# Substitute
	locis$Chr[indexes]   <<- vlocis$Chr[i]
	locis$Start[indexes] <<- vlocis$Start[i]
	locis$End[indexes]   <<- vlocis$End[i]
	#locis$NodeID         <<- paste("v",vlocis$NodeID[i],sep=".")
}))

# Write new network
write.table(locis, file = opt$output, row.names = F, col.names = F, sep = "\t", quote = F)