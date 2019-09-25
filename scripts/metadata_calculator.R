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
  make_option(c("-l", "--locis"), action="store", type="character", default = NULL,
              dest="locis", help="File with locis extended info (genes included)"),
  make_option(c("-g", "--go"), action="store", type="character", default = NULL,
              dest="go", help="File with GO enrichment (special cluster profiler format)"),
  make_option(c("-k", "--kegg"), action="store", type="character", default = NULL,
              dest="kegg", help="File with KEGG enrichment (cluster profiler format)"),
  make_option(c("-r", "--reac"), action="store", type="character", default = NULL,
              dest="reac", help="File with REACTOME enrichment (cluster profiler format)"),
  make_option(c("-v","--verbose"), action="store_true",type="logical", default = FALSE,
              dest="verbose", help="Activate verbose mode"),
  make_option(c("-o","--output"), action="store",type="character", default = "meta",
              dest="output", help="Destiny file base name. [Default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$locis) & is.null(opt$go) & is.null(opt$kegg)){
	stop("ERROR [metadata_calculator]: any input file given")
}

if(opt$verbose){
	suppressPackageStartupMessages(require(pbapply))
	pboptions(type="timer")
	apf <- pblapply
}else{
	apf <- lapply
}

options(stringsAsFactors = F)

####################################################################
##                               LOAD                             ##
####################################################################
# Load all data files and reset target files
if(!is.null(opt$locis)){
	locis <- read.table(file = opt$locis, sep = "\t", quote = "", header = T, stringsAsFactors = F)	
	write.table(paste(c("Model","HPO","Locis","LocisWithGenes","Genes","MeanLocisGenes","GenesList","LocisList","LocisWithGenesList"),collapse = "\t"),
			file = paste(paste(opt$output, "_loci", sep = ""), sep = "."), 
	        sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}
if(!is.null(opt$go)){
	go    <- read.table(file = opt$go, sep = "\t", quote = "", header = T, stringsAsFactors = F)
	go$Genes <- as.character(go$Genes)
	go$Fisher.Test <- as.double(gsub("<","",as.character(go$Fisher.Test)))
	write.table(paste(c("Model","HPO","Pval","BP","CC","MF","Genes","GenesPercentageOfTotal","GenesList","SignalList"),collapse = "\t"),
			file = paste(paste(opt$output, "_go", sep = ""), sep = "."), 
	        sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}
if(!is.null(opt$kegg)){
	kegg  <- read.table(file = opt$kegg, sep = "\t", quote = "", header = T, stringsAsFactors = F)
	write.table(paste(c("Model","HPO","Pval","Terms","Genes","GenesPercentageOfTotal","GenesList","SignalList"),collapse = "\t"),
			file = paste(paste(opt$output, "_kegg", sep = ""), sep = "."), 
	        sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}

if(!is.null(opt$reac)){
	reac  <- read.table(file = opt$reac, sep = "\t", quote = "", header = T, stringsAsFactors = F)
	write.table(paste(c("Model","HPO","Pval","Terms","Genes","GenesPercentageOfTotal","GenesList","SignalList"),collapse = "\t"),
			file = paste(paste(opt$output, "_reac", sep = ""), sep = "."), 
	        sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}


####################################################################
##                        CONFIGURE PROG                          ##
####################################################################
# Obtain list of locis
if(!is.null(opt$loci)){
	last_index <- 1
	max_freq   <- max(table(locis$Loci))
	lsize      <- nrow(locis)
	# Obtain LOCIS metadata per CNV's IDs
	locis_nodes <- as.data.frame(do.call(rbind,apf(unique(locis$Loci),function(node){
		# Find 
		if(last_index + max_freq + 1 > lsize){
			indx <- which(locis$Loci[last_index:lsize] == node)[1] + last_index - 1			
		}else{
			indx <- which(locis$Loci[last_index:(last_index+max_freq+1)] == node)[1] + last_index - 1
		}
		if(is.na(indx)){
			indx <- which(locis$Loci == node)[1]
		}
		last_index <<- indx
		# Obtain genes info
		num_genes <- 0
		genes <- NA
		if(!is.na(locis$Genes[indx]) | nchar(locis$Genes[indx]) > 0){
			genes     <- unlist(strsplit(locis$Genes[indx],":"))
			num_genes <- length(genes) 
		}
		# Create data frame
		info <- data.frame(NodeID  = node, 
			               NumGenes = num_genes,
			               Genes    = NA,
			               stringsAsFactors = F)
		# Check
		if(num_genes > 0){
			info$Genes <- list(genes)
		}

		# Return
		return(info)
	})))	
}


## Prepare necessary variables
# P-value thresholds
pval_thrs <- c(1e-3,1e-4,1e-5)

# Obtain HPOs and models
hpos   <- c()
models <- c()
if(!is.null(opt$loci)){
    hpos   <- c(hpos,locis$HPO)
    models <- c(models,locis$tag)
}
if(!is.null(opt$go)){
	if('hpo' %in% go$Term){
		hpos   <- c(hpos,go$Term[-which(go$Model == 'hpo')])
	}else{
		hpos   <- c(hpos,go$Term)
	}
	models <- c(models,go$Model)

	# Obtain subset
	go <- go[which(go$Fisher.Test <= pval_thrs[1]),]
	# go <- go[-which(go$Model == 'hpo'),]
}
if(!is.null(opt$kegg)){
    hpos   <- c(hpos,kegg$Term[-which(kegg$Model == 'hpo')])
	models <- c(models,kegg$Model)

	# Obtain subset
	kegg <- kegg[which(kegg$p.adjust <= pval_thrs[1]),]
}

if(!is.null(opt$reac)){
    hpos   <- c(hpos,reac$Term[-which(reac$Model == 'hpo')])
	models <- c(models,reac$Model)

	# Obtain subset
	kegg <- kegg[which(reac$p.adjust <= pval_thrs[1]),]
}


hpos   <- unique(hpos)
models <- unique(models)
if('hpo' %in% models){
	models <- models[-which(models == 'hpo')]
}


####################################################################
##                         CALC METADATA                          ##
####################################################################

# Per each Model, obtain metadata info
invisible(apf(models, function(model){
	############################ LOCIS
	if(!is.null(opt$loci) & !grepl("rdm_g_",model)){
		# Obtain target idnex
		indx_loci <- which(locis$tag == model)	
		tmp_locis <- locis[indx_loci,]
		# Obtain LOCIS metadata : {Number of locis per HPO, Number of locis with genes, 
		#	                       Number of different genes, Mean of genes per loci}
		meta_locis <- as.data.frame(do.call(rbind,lapply(hpos,function(hp){
			# Find HPO indexes
			indx <- which(tmp_locis$HPO == hp)
			# Obtain number of locis
			cnvs <- unique(tmp_locis$Loci[indx])
			# Obtain locis info
			e_locis <- which(locis_nodes$NodeID %in% cnvs & !is.na(locis_nodes$Genes))
			# Obtain genes list
			genes <- unique(as.vector(na.omit(unlist(locis_nodes$Genes[e_locis]))))
			# Obtain genes metadata (Genes per loci, Number of locis with genes)
			num_genes <- unlist(lapply(locis_nodes$Genes[e_locis],function(genesl){
				if(length(genesl) == 1 & is.na(genesl[1])) return(0)
				else return(length(unlist(genesl)))
			}))
			# Return info
			return(data.frame(Model          = model,
				              HPO            = hp,
				              Locis          = length(cnvs),
				              LocisWithGenes = length(e_locis),
				              Genes          = length(genes),
				              MeanLocisGenes = mean(num_genes),
				              GenesList      = paste(genes,collapse = ":"),
				              LocisList      = paste(cnvs,collapse = ":"),
				              LocisWGList    = paste(locis_nodes$NodeID[e_locis],collapse = ":"),
				              stringsAsFactors = F))
		})))

		# Write all info
		write.table(meta_locis, file = paste(paste(opt$output, "_loci", sep = ""), sep = "."), sep = "\t", append = T, quote = F, col.names = F, row.names = F)
	}


	############################ GO
	if(!is.null(opt$go)){
		# Obtain target indexes
		indx_go <- which(go$Model == model)
		tmp_go  <- go[indx_go,]


		# GO types indexes
		indx_go_bp <- which(tmp_go$Source == "BP")
		indx_go_cc <- which(tmp_go$Source == "CC")
		indx_go_mf <- which(tmp_go$Source == "MF")

		# Obtain GO metadata : {Number of terms enriched of each type, Num genes which enrich something,
		#                       % of genes which enrich something}
		meta_go <- as.data.frame(do.call(rbind, lapply(hpos,function(hp){
			# Find HPO indexes
			indx <- which(tmp_go$Term == hp)
			# Obtain genes list
			if(exists("meta_locis") & !grepl("rdm_g_",model)){
				#genes <- unlist(strsplit(meta_locis$GenesList[which(meta_locis$HPO == hp)],":"))
				genes <- meta_locis$Genes[which(meta_locis$HPO == hp)]
			}else{
				genes <- NULL
			}

			# Obtain data per each possible_pval
			info <- as.data.frame(do.call(rbind, lapply(pval_thrs,function(thr){
				# Find targets
				indx_thr <- intersect(which(tmp_go$Fisher.Test <= thr), indx)
				# Obtain number of BP, CC and MF terms
				num_bp <- intersect(indx_thr,indx_go_bp)
				num_cc <- intersect(indx_thr,indx_go_cc)
				num_mf <- intersect(indx_thr,indx_go_mf)
				# Concat signals
				signals <- paste(tmp_go$GO[c(num_bp,num_cc,num_mf)],collapse = ";")
				# Obtain genes which enrich something
				e_genes <- unique(unlist(strsplit(tmp_go$Genes[indx_thr],":")))
				# Percentage
				if(is.null(genes)){
					percent <- NA
				}else if(genes == 0){
					percent <- 0
				}else{
					percent <- length(e_genes) / genes
				}
				# Return info
				return(data.frame(Model = model,
					              HPO   = hp,
					              Pval  = thr,
					              BP    = length(num_bp),
					              CC    = length(num_cc),
					              MF    = length(num_mf),
					              Genes = length(e_genes),
					              GenesPercentageOfTotal = percent,
					              GenesList              = paste(e_genes,collapse = ":"),
					              SignalList             = signals,
					              stringsAsFactors = F))
			})))
			return(info)
		})))

		# Write all info
		write.table(meta_go, file = paste(paste(opt$output, "_go", sep = ""), sep = "."), sep = "\t", append = T, quote = F, col.names = F, row.names = F)
	}

	if(!is.null(opt$kegg)){
		# Obtain target indexes
		indx_kegg <- which(kegg$Model == model)
		tmp_kegg <- kegg[indx_kegg,]
		# Obtain KEGG metadata : {Number of terms enriched, Num genes which enrich something,
		#                         % of genes which enrich something}
		meta_kegg <- as.data.frame(do.call(rbind, lapply(hpos,function(hp){
			# Find HPO indexes
			indx <- which(tmp_kegg$Set == hp)
			# Obtain genes list
			if(exists("meta_locis") & !grepl("rdm_g_",model)){
                #genes <- unlist(strsplit(meta_locis$GenesList[which(meta_locis$HPO == hp)],":"))
                genes <- meta_locis$Genes[which(meta_locis$HPO == hp)]
            }else{
                genes <- NULL
            }

			# Obtain data per each possible_pval
			info <- as.data.frame(do.call(rbind, lapply(pval_thrs,function(thr){
				# Find targets
				indx_thr <- intersect(which(tmp_kegg$p.adjust <= thr), indx)
				# Concat signals
				signals <- paste(tmp_kegg$ID[indx_thr],collapse = ";")
				# Obtain genes which enrich something
				e_genes <- unique(unlist(strsplit(tmp_kegg$geneID[indx_thr],"/")))
				# Percentage
                if(is.null(genes)){
                        percent <- NA
                }else if(genes == 0){
                        percent <- 0
                }else{
                        percent <- length(e_genes) / genes
                }


				# Return info
				return(data.frame(Model = model,
					              HPO   = hp,
					              Pval  = thr,
					              Terms = length(indx_thr),
					              Genes = length(e_genes),
					              GenesPercentageOfTotal = percent,
					              GenesList              = paste(e_genes,collapse = ":"),
					              SignalList             = signals,
					              stringsAsFactors = F))
			})))
			return(info)
		})))	

		# Write all info
		write.table(meta_kegg, file = paste(paste(opt$output, "_kegg", sep = ""), sep = "."), sep = "\t", append = T, quote = F, col.names = F, row.names = F)
	}

	if(!is.null(opt$reac)){
		# Obtain target indexes
		indx_reac <- which(reac$Model == model)
		tmp_reac <- reac[indx_kegg,]
		# Obtain REACTOME metadata : {Number of terms enriched, Num genes which enrich something,
		#                         % of genes which enrich something}
		meta_reac <- as.data.frame(do.call(rbind, lapply(hpos,function(hp){
			# Find HPO indexes
			indx <- which(tmp_reac$Set == hp)
			# Obtain genes list
			if(exists("meta_locis") & !grepl("rdm_g_",model)){
                #genes <- unlist(strsplit(meta_locis$GenesList[which(meta_locis$HPO == hp)],":"))
                genes <- meta_locis$Genes[which(meta_locis$HPO == hp)]
            }else{
                genes <- NULL
            }
			# Obtain data per each possible_pval
			info <- as.data.frame(do.call(rbind, lapply(pval_thrs,function(thr){
				# Find targets
				indx_thr <- intersect(which(tmp_reac$p.adjust <= thr), indx)
				# Concat signals
				signals <- paste(tmp_reac$ID[indx_thr],collapse = ";")
				# Obtain genes which enrich something
				e_genes <- unique(unlist(strsplit(tmp_reac$geneID[indx_thr],"/")))
				# Percentage
                if(is.null(genes)){
                        percent <- NA
                }else if(genes == 0){
                        percent <- 0
                }else{
                        percent <- length(e_genes) / genes
                }

				# Return info
				return(data.frame(Model = model,
					              HPO   = hp,
					              Pval  = thr,
					              Terms = length(indx_thr),
					              Genes = length(e_genes),
					              GenesPercentageOfTotal = percent,
					              GenesList              = paste(e_genes,collapse = ":"),
					              SignalList             = signals,
					              stringsAsFactors = F))
			})))
			return(info)
		})))	

		# Write all info
		write.table(meta_reac, file = paste(paste(opt$output, "_reac", sep = ""), sep = "."), sep = "\t", append = T, quote = F, col.names = F, row.names = F)
	}
}))


