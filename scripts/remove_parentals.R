#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(GO.db))

option_list <- list(
  make_option(c("-f", "--input"), action="store", type="character",
              dest="input", help="Input table with necessary info"),
  make_option(c("-i", "--indexes"), action="store",type="character",
              dest="indexes", help="String with Signal(GO), GroupBy1(Model) and GroupBy2(HPO) column indexes, separated by commas"),
  make_option(c("-r", "--remove"), action="store",type="character", default = NULL,
              dest="remove", help="Model tags to be removed separated by commas. [Default: '']"),
  make_option(c("-H", "--header"), action="store_true",type="logical", default = FALSE,
              dest="header", help="Flag which indicates that file has header. [Default: %default]"),
  make_option(c("-v", "--verbose"), action="store_true",type="logical", default = FALSE,
              dest="verbose", help="Flag which activates verbose mode. [Default: %default]"),  
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store ratios calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################

# Load info
info <- read.table(file = opt$input, header = opt$header, sep = "\t", stringsAsFactors = F, quote = "")

# Prepare indexes
indx <- as.double(unlist(strsplit(opt$indexes,",")))

# remove
if(!is.null(opt$remove)){
	# Find
	to_remove <- which(info[,indx[2]] %in% unlist(strsplit(opt$remove,",")))
	if(length(to_remove) > 0){
		info <- info[-to_remove,]
		if(opt$verbose){
			message(paste("Enrichments removed by model:",length(to_remove)))
		}
	}
}

################################################################
##                        TRANSFORM DATA                      ##
################################################################


# Prepare groupBy IDs
groupBy1 <- unique(info[,indx[2]])
groupBy2 <- unique(info[,indx[3]])

# Check
if(length(groupBy2) > length(groupBy1)){
	# Swap ids
	aux <- groupBy1
	groupBy1 <- groupBy2
	groupBy2 <- aux
	# Swap indexes
	aux <- indx[3]
	indx[3] <- indx[2]
	indx[2] <- aux
}

if(opt$verbose){
	suppressPackageStartupMessages(require(pbapply))
	pboptions(type="timer")
	apf <- pblapply
}else{
	apf <- lapply
}


if(opt$verbose){
	message("Calculating speedup variables")
}

# Prepare speedup indexes
#indx_gb1 <- lapply(groupBy1,function(gb){return(which(info[,indx[2]] == gb))})
#indx_gb2 <- lapply(groupBy2,function(gb){return(which(info[,indx[3]] == gb))})
indx_gb1 <- apf(groupBy1,function(gb){return(which(info[,indx[2]] == gb))})
indx_gb2 <- apf(groupBy2,function(gb){return(which(info[,indx[3]] == gb))})


# Prepare GO terms allowed
go_terms <- names(Term(GOTERM))

# Per each GO term into our set, find parentals
go_targets <- unique(info[,indx[1]])

if(opt$verbose){
	message(paste("GO terms enriched:",length(go_targets)))
	message("Preparing GO terms ancestors map")
}

#go_parentals <- as.data.frame(do.call(rbind,lapply(go_targets,function(go){
go_parentals <- as.data.frame(do.call(rbind,apf(go_targets,function(go){
	ont <- Ontology(GOTERM[[go]])
	if(ont == "BP"){
		p <- GOBPANCESTOR[[go]]
	}else if(ont == "MF"){
		p <- GOMFANCESTOR[[go]]
	}else{ # CC
		p <- GOCCANCESTOR[[go]]
	}
	df <- data.frame(GO = go, Parentals = "", stringsAsFactors = F)
	df$Parentals <- list(p)
	return(df)
})))

if(opt$verbose){
	message("Finding parentals to be removed")
}


# Per each tuple, remove parental signals
#to_remove <- unlist(lapply(seq_along(groupBy1),function(i){
to_remove <- unlist(apf(seq_along(groupBy1),function(i){
	# Check with it's tuples
	to_remove <- unlist(lapply(seq_along(groupBy2),function(j){
		# Obtain indexes
		indx_t <- intersect(indx_gb1[[i]],indx_gb2[[j]])
		# Check
		if(length(indx_t) > 0){
			# Find parentals
			parentals <- unique(unlist(lapply(indx_t,function(z){
				# Check special case
				if(!info[z,indx[1]] %in% go_terms){
					return(z)
				}else{
					# Obtain term
					return(indx_t[which(info[indx_t,indx[1]] %in% go_parentals$Parentals[[which(go_parentals$GO == info[z,indx[1]])]])])
				}
			})))
			return(parentals)
		}
		# ELSE: nothing to do here
		return(integer())
	}))
	return(to_remove)
}))

if(opt$verbose){
	message(paste("Rows removed:",length(to_remove),"/",nrow(info)))	
}

# Remove lines
if(length(to_remove) > 0){
	info <- info[-to_remove,]
}

################################################################
##                          WRITE DATA                        ##
################################################################

write.table(info, file = opt$output, sep = "\t", col.names = opt$header, row.names = F, quote = F)