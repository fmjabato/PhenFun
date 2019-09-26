#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
require(optparse)
require(ontologyIndex)


option_list <- list(
  make_option(c("-e", "--enrich"), action="store", type="character",
              dest="enrich", help="File with filtered touples Pheno-System and rdm model used to filter"),
  make_option(c("-p", "--pheno_onto"), action="store", type="character", default = NULL,
              dest="pheno_onto", help="OBO file to load phenotypes ontology. Default = ontologyIndex HPO ontology"),
  make_option(c("-f", "--funsys_onto"), action="store", type="character",default = NULL,
              dest="funsys_onto", help="OBO file to load used FunSys ontology. Default = ontologyIndex GO ontology"),
  make_option(c("-i", "--indexes"), action="store", type="character",default = "1,2,4",
              dest="indexes", help="Column indexes separated by commas <Phenotype,IC,FunSys>. [Default = %default]"),
  make_option(c("-s", "--sep"), action="store", type="character",default = ";",
              dest="sep", help="FunSys signal separator into enrichment file. [Default = %default]"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file")
)

opt <- parse_args(OptionParser(option_list=option_list))



################################################################
##                        LOAD AND PARSE                      ##
################################################################
# Load enrichments
enrichments <- read.table(file = opt$enrich, sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)

# Obtain columns indexes
indexes <- as.numeric(unlist(strsplit(opt$indexes,",")))
indexes <- list(Pheno = indexes[1], IC = indexes[2], FunSys = indexes[3])

# Load Phenotypes onology
if(!is.null(opt$pheno_onto)){
	pheno_onto <- get_ontology(file = opt$pheno_onto, propagate_relationships = "is_a")
}else{
	data(hpo)
	pheno_onto <- hpo
}

# Load FunSys ontology
if(!is.null(opt$funsys_onto)){
	funsys_onto <- get_ontology(file = opt$funsys_onto, propagate_relationships = "is_a")
}else{
	data(go)
	funsys_onto <- go
}

# Obtain HPO-GO pairs
pairs_enr <- as.data.frame(do.call(rbind,lapply(seq(nrow(enrichments)),function(i){
	# Obtain Phenotype and it's IC
	pheno <- enrichments[i,indexes$Pheno]
	ic    <- enrichments[i,indexes$IC]
	# Obtain signal
	funsys <- unlist(strsplit(enrichments[i,indexes$FunSys],opt$sep))
	# Return Info
	return(data.frame(Pheno    = rep(pheno,length(funsys)),
					  IC       = rep(ic,length(funsys)),
					  FunSys   = funsys,
					  ToRemove = rep(FALSE,length(funsys)),
					  stringsAsFactors = FALSE))
})))

# Sort by IC (most specific first)
pairs_enr <- pairs_enr[order(-pairs_enr$IC),]




################################################################
##                         CHECK PAIRS                        ##
################################################################
# Check all pairs starting by most specific
phenotypes <- unique(pairs_enr$Pheno)
invisible(lapply(phenotypes,function(pheno){
	# Obtain phenotype parentals
	parentals <- head(get_ancestors(pheno_onto,pheno),-1)
	
	# Check
	if(any(parentals %in% phenotypes)){
		parentals <- parentals[which(parentals %in% phenotypes)]
	}else{ # Any parental to be checked
		return(NULL)
	}
	
	# Obtain active pairs
	funsys <- pairs_enr$FunSys[which(pairs_enr$Pheno == pheno & !pairs_enr$ToRemove)]
	
	# Check
	if(length(funsys) <= 0){ # Any funsys to be checked
		return(NULL)
	}

	# Obtain parental pairs
	indx <- which(pairs_enr$Pheno %in% parentals & !pairs_enr$ToRemove)
	
	# Check
	if(length(indx) > 0){
		pairs_parental <- pairs_enr$FunSys[indx]		
	}else{
		return(NULL)
	}

	# Per each funsys check if parentals have repetitive info
	to_remove <- unique(unlist(lapply(funsys,function(fs){
		# Obtain FunSys parentals
		fs_family <- get_ancestors(funsys_onto,fs)
		# Check if already exists into phenotype parentals
		redundant_pairs <- which(pairs_parental %in% fs_family)
		if(length(redundant_pairs) > 0){
			return(redundant_pairs)
		}else{
			return(NA)
		}
	})))
	if(any(is.na(to_remove))){
		to_remove <- to_remove[-which(is.na(to_remove))]
	}

	# Check and update
	if(length(to_remove) > 0){
		pairs_enr$ToRemove[indx[to_remove]] <<- TRUE
	}
}))


# Check init stats
total_phenos <- length(unique(pairs_enr$Pheno))
total_funsys <- length(unique(pairs_enr$FunSys))
total_pairs <- nrow(pairs_enr)

# Update set
pairs_enr <- pairs_enr[-which(pairs_enr$ToRemove),]
pairs_enr <- pairs_enr[,-ncol(pairs_enr)]

# Check final stats
final_phenos <- length(unique(pairs_enr$Pheno))
final_funsys <- length(unique(pairs_enr$FunSys))
final_pairs <- nrow(pairs_enr)


# Verbose point
message(paste("STATS: AFTER/BEFORE\nDifferent phenotypes: ",final_phenos,"/",total_phenos,
			  "\nDifferent FunSys: ",final_funsys,"/",total_funsys,
			  "\nPheno-FunSys pairs: ",final_pairs,"/",total_pairs,sep=""))


# Return to collapsed format
pairs_enr <- as.data.frame(do.call(rbind,lapply(unique(pairs_enr$Pheno),function(pheno){
	# Find
	indx <- which(pairs_enr$Pheno == pheno)
	# Collapse info
	return(data.frame(Pheno        = pheno,
					  IC           = pairs_enr$IC[indx[1]],
					  SignalAmount = length(indx),
					  Signal       = paste(pairs_enr$FunSys[indx],collapse = opt$sep),
					  stringsAsFactors = FALSE))
})))


################################################################
##                          WRITE DATA                        ##
################################################################
write.table(pairs_enr, file = opt$output, col.names = T, row.names = F, sep = "\t", quote = F)