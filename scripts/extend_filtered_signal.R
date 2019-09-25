#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
require(optparse)

option_list <- list(
  make_option(c("-s", "--signal"), action="store", type="character",
              dest="signal", help="Filtered signal file"),
  make_option(c("-r", "--ratios"), action="store", type="character", default = NULL,
              dest="ratios", help="Ratios file. If not specified, ratio value will be NA"),
  make_option(c("-f", "--freqs"), action="store", type="character",
              dest="freqs", help="HPOs frequencies used to calculate HPO IC"),
  make_option(c("-t", "--threshold"), action="store", type="double", default = 1e-3,
              dest="threshold", help="Ratios threshold set. [Default: %default]"),
  make_option(c("-m", "--model"), action="store",type="character", default = "rdm_l_",
              dest="model", help="Random model (flag) used to filter. [Default: %default]"),
  make_option(c("-T", "--type"), action="store",type="character", default = "GO",
              dest="type", help="Signal type. [Default: %default]"),
  make_option(c("-a", "--append"), action="store_true",type="logical", default = FALSE,
              dest="append", help="Append results to output file. [Default: %default]"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store ratios calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
 # Load filtered signal
network <- read.table(file = opt$signal, header = F, sep = "\t", stringsAsFactors = F)
colnames(network) <- c("Model","HPO","Signal","Pval")
network <- network[network$Model == opt$model,]

# Load original ratios
if(!is.null(opt$ratios)){
  ratios <- read.table(file = opt$ratios, header = T, sep = "\t", stringsAsFactors = F)
  ratios <- ratios[which(ratios$Signal_type == opt$type & ratios$RDM_type == gsub("_$","",opt$model) & ratios$Pval_thr == opt$threshold),]  
}else{
  ratios <- NULL
}

# Load original DECIPHER HPO frequencies
hpos <- read.table(file = opt$freqs, header = F, sep = "\t", stringsAsFactors = F)
colnames(hpos) <- c("HPO","Freq")

# Load HPO ontology and add frequencies to obtain ICs
#data(hpo)
#ics <- get_term_info_content(hpo,unlist(lapply(seq_along(hpos$HPO),function(i){return(rep(hpos$HPO[i],hpos$Freq[i]))})))
ics <- -log(hpos$Freq/max(hpos$Freq))
names(ics) <- hpos$HPO

################################################################
##                         COMBINE DATA                       ##
################################################################
# Generate DF with <Model> <PvalThr> <HPO> <IC> <Ratio> <SignalAmount> <Signal>
info <- as.data.frame(do.call(rbind,lapply(unique(network$HPO),function(term){
	# Obtain signal
	signal <- network$Signal[which(network$HPO == term)]
	# Return info
	return(data.frame(Model = opt$model,
		              PvalThr = opt$threshold,
		              HPO = term,
		              IC = ics[term],
		              Ratio = ifelse(is.null(ratios),NA,ratios$Ratio[which(ratios$HPO == term)]),
		              SignalAmount = length(signal),
		              Signal = paste(signal,collapse = ";"),
		              stringsAsFactors = F))
})))


################################################################
##                          WRITE DATA                        ##
################################################################
write.table(info, file = opt$output, col.names = !opt$append, row.names = F, sep = "\t", quote = F, append = opt$append)