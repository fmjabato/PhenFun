#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-e", "--enrich"), action="store", type="character",
              dest="enrich", help="Enrichment table with necessary info"),
  make_option(c("-m", "--model"), action="store", type="character", default = "dec",
              dest="model", help="Model ID to be compared against rest of models. If it's a file, full file info will be used to compare data with HPO;IC;SignalAmount;Signal header. [Default: %default]"),
  make_option(c("-A", "--avoid"), action="store", type="character", default = "hpo",
              dest="avoid", help="List of models to be avoided separated by commas. [Default: %default]"),
  make_option(c("-r", "--randoms"), action="store", type="character",
              dest="randoms", help="File with Random models info"),
  make_option(c("-s", "--signal"), action="store", type="character",
              dest="signal", help="String with signal type to be used as flag in final table"),
  make_option(c("-i", "--indexes"), action="store",type="character", default = "1,2,3,7",
              dest="indexes", help="String with Model, Term, Signal and PVal column indexes, separated by commas [Default: %default]"),
  make_option(c("-t","--thresholds"), action="store",type="character", default = "1e-3,1e-4,1e-5",
              dest="thresholds", help="Thresholds to be applied over Pvalues, separated by commas. [Default: %default]"),
  make_option(c("-a","--append"), action="store_true",type="logical", default = FALSE,
              dest="append", help="Activate append mode to not overwrite output file"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store ratios calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
# Load enrichments
enrichment <- read.table(file = opt$enrich, header = T, sep = "\t", quote = "", stringsAsFactors = F)

# Prepare indexes, thresholds and avoid models
indx <- as.double(unlist(strsplit(opt$indexes,",")))
indx <- list(Model = indx[1], HPO = indx[2], Signal = indx[3], Pval = indx[4])

thrs  <- as.double(unlist(strsplit(opt$thresholds,",")))

# Load random models info
rdms <- read.table(file = opt$randoms, header = T, sep = "\t", quote = "", stringsAsFactors = F)
rdm_sizes <- unlist(lapply(rdms$Format,function(model){
	# Identify
	indx_rdm <- which(grepl(paste(model,"[0-9]",sep=""),enrichment[,indx$Model]))
	# Obtain different iddentifiers
	return(length(unique(enrichment[indx_rdm,indx$Model])))
}))

# Remove avoid models
avoid <- unlist(strsplit(opt$avoid,","))
avoid <- which(enrichment[,indx$Model] %in% avoid)
if(length(avoid) > 0){
	enrichment <- enrichment[-avoid,]	
}

if(file.exists(opt$model)){
	# Load file
	target_info <- read.table(file = opt$model, sep = "\t", quote = "", header = T, stringsAsFactors = F)
	# Change format to fit current enrichments
	target_info <- as.data.frame(do.call(rbind,lapply(seq(nrow(target_info)),function(i){
		signals <- unlist(strsplit(target_info$Signal[i],";"))
		return(data.frame(Model  = rep("dec",length(signals)),
			              HPO    = rep(target_info$HPO[i],length(signals)),
			              Signal = signals,
			              Pval   = rep(min(thrs),length(signals)),
			              stringsAsFactors = F))
	})))

	# Transform enrichments set to fit
	enrichment <- enrichment[,unlist(indx)]

	# Append
	colnames(target_info) <- colnames(enrichment)
	enrichment <- rbind(target_info,enrichment)

	# Update target model id and indexes
	opt$model <- "dec"
	indx <- list(Model = 1, HPO = 2, Signal = 3, Pval = 4)
}


# Obtain target and remove models to be avoided
indx_t <- which(enrichment[,indx$Model] == opt$model)


# Obtain list of models
models <- unique(rdms$Name)

# Obtain list of Terms
hpos <- unique(enrichment[,indx$HPO])


################################################################
##                        TRANSFORM DATA                      ##
################################################################
# Prepare speedup indexes
indx_thr   <- lapply(thrs,function(thr){return(which(enrichment[,indx$Pval] <= thr))})
indx_model <- lapply(rdms$Format,function(model){return(which(grepl(paste(model,"[0-9]",sep=""),enrichment[,indx$Model])))})
indx_hpos  <- lapply(hpos, function(hp){return(which(enrichment[,indx$HPO] == hp))})

#require(pbapply)
#pboptions(type="timer")

# Per each hpo, compare target model against each other
#dfratios <- as.data.frame(do.call(rbind,pblapply(hpos,function(hp){
dfratios <- as.data.frame(do.call(rbind,lapply(hpos,function(hp){
	# Find target
	indx_t_hpo <- intersect(indx_t,indx_hpos[[which(hpos == hp)]])

	# Obtain target signal
	signal_t <- unlist(lapply(thrs,function(thr){
		signal <- length(intersect(indx_t_hpo,indx_thr[[which(thrs == thr)]])) 
		return(ifelse(signal == 0,NA,signal))
	}))
	# Find NAs
	nas_t <- is.na(signal_t)

	# Compare agains rest of models
	info <- as.data.frame(do.call(rbind,lapply(models,function(model){
		# Find model
		indx_m_hpo <- intersect(indx_model[[which(models == model)]],indx_hpos[[which(hpos == hp)]])
		# Find signals
		signal_m <- unlist(lapply(thrs,function(thr){
			# Find
			indx_current <- intersect(indx_m_hpo,indx_thr[[which(thrs == thr)]])
			# Obtain subset
			aux <- enrichment[indx_current,]
			# Obtain signals per each RDM model
			signal <- sum(unlist(lapply(unique(aux[,indx$Model]),function(rdmodel){
				# Find signal
				return(length(which(aux[,indx$Model] == rdmodel)))
			})))
			signal <- signal / rdm_sizes[which(rdms$Name == model)]
			return(ifelse(signal == 0,NA,signal))
		}))
		# Find NAs
		nas_m <- is.na(signal_m)
		# Check special case
		if(all(c(nas_t,nas_m))) return(data.frame())
		# Obtain final info
		ratios <- as.data.frame(do.call(rbind,lapply(seq_along(signal_t),function(i){
			# Calc
			if(nas_t[i] | nas_m[i]){
				ratio <- NA
			}else{
				ratio <- log(signal_t[i] / signal_m[i])
			}
			# Return
			return(data.frame(Signal_type = opt$signal,
							  RDM_type    = model,
							  HPO         = hp,
							  Pval_thr    = thrs[i],
							  Target      = signal_t[i],
							  RDM         = signal_m[i],
							  Ratio       = ratio, stringsAsFactors = F))
		})))
		# Return info
		return(ratios)
	})))
	# Return info
	return(info)
})))


################################################################
##                          WRITE DATA                        ##
################################################################

# Write output
write.table(dfratios,file = opt$output, row.names = FALSE, quote = FALSE, sep = "\t", col.names = !opt$append, append = opt$append)