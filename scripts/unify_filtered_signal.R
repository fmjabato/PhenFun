#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
require(optparse)

option_list <- list(
  make_option(c("-e", "--enrich"), action="store", type="character",
              dest="enrich", help="File with filtered touples Pheno-System and rdm model used to filter"),
  make_option(c("-r", "--random"), action="store", type="character",
              dest="random", help="Random models table info"),
  make_option(c("-t", "--threshold"), action="store", type="integer", default = 3,
              dest="threshold", help="Number of random models which must include the tuple to be accepted. [Default: %default]"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
# Load filtered signal
enrichs <- read.table(file = opt$enrich, sep = "\t", quote = "", header = T, stringsAsFactors = F)

# Load RDM info
rdms_info <- read.table(file = opt$enrich, sep = "\t", quote = "", header  = T, stringsAsFactors = F)

# Filter enrichments by rdm formats
enrichs <- enrichs[which(unlist(lapply(enrichs$Model,function(rdm){return(rdm %in% rdms_info$Format)}))),]

# Prepare rdm models identifiers
rdm_models <- unique(enrichs$Model)

################################################################
##                         COMBINE DATA                       ##
################################################################
# Create tuples membership matrix
tuples_by_model <- lapply(rdm_models,function(model){
  # Per term, obtain tuples
  tuples <- unlist(lapply(which(enrichs$Model == model),function(i){
    # Generate tuples
    return(paste(enrichs$HPO[i],enrichs$IC[i],unlist(strsplit(enrichs$Signal[i],";")),sep=";"))
  }))
  tuples
})

# Obtain unique tuples
tuples <- unique(unlist(tuples_by_model))

# Prepare membership matrix
membership <- data.frame(Tuple = tuples, stringsAsFactors = F)
for(i in seq_along(rdm_models)){
  membership[,as.character(rdm_models[i])] <- membership$Tuple %in% tuples_by_model[[i]]
}

# Calculate membership ratio
membership$Validations <- unlist(lapply(seq(nrow(membership)),function(i){return(sum(as.numeric(membership[i,-1])))}))

# Take only which have enough validations
indx <- which(membership$Validations >= opt$threshold)

# Prepare output format
info <- as.data.frame(do.call(rbind,lapply(membership$Tuple[indx],function(tuple){
  items <- unlist(strsplit(tuple,";"))
  return(data.frame(HPO = items[1], IC = items[2], Signal = items[3], stringsAsFactors = F))
})))
info <- as.data.frame(do.call(rbind,lapply(unique(info$HPO),function(hp){
  # Take data
  signals <- info$Signal[which(info$HPO == hp)]
  ic <- info$IC[which(info$HPO == hp)][1]
  # Return
  return(data.frame(HPO          = hp, 
                    IC           = ic, 
                    SignalAmount = length(signals), 
                    Signal       = paste(signals,collapse = ";"), 
                    stringsAsFactors = F))
})))

# Message
message(paste("\tTuples (Phenotype-System): Before / After :",length(tuples),"/",length(indx)))
message(paste("\tPhenotypes with signals: Before / After :",length(unique(enrichs$HPO)),"/",nrow(info)))

################################################################
##                          WRITE DATA                        ##
################################################################
write.table(info, file = opt$output, col.names = T, row.names = F, sep = "\t", quote = F)