#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
suppressPackageStartupMessages(require(optparse))

dicts_ext <- "_dict"

option_list <- list(
  make_option(c("-g", "--go"), action="store", type="character", default = NULL,
              dest="go", help="GO final unified enrichment set"),
  make_option(c("-k", "--kegg"), action="store", type="character", default = NULL,
              dest="kegg", help="KEGG final unified enrichment set"),
  make_option(c("-r", "--reac"), action="store", type="character", default = NULL,
              dest="reac", help="Reactome final unified enrichment set"),
  make_option(c("-a", "--all"), action="store_true", default = FALSE,
              dest="all", help="Flag used to indicate that max frequency must be the number of genes with any FunSys term linked instead the maximum frequency observed"),
  make_option(c("-s", "--sep"), action="store", type="character", default = ";",
              dest="sep", help="Signal list separator. [Default: %default]"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help=paste("Output file basename to store dictionaries. Extension '<FunSysType>_",dicts_ext,"' will be added",sep=""))
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           FUNCTIONS                        ##
################################################################

#' @description function used to obtain unique list of enriched FunSys terms from a results dataset 
#' @param unified dataframe with results set
#' @param sep sperator character used in signal list
#' @return a list of unique terms contained into signal lists along the resutls set
unified2signalList <- function(unified,sep=opt$sep){
	# Obtain unique list of enriched terms
	signal_list <- unique(unlist(lapply(seq(nrow(unified)),function(i){
		# Take signal list
		sig <- unified$Signal[i]
		# Decompose
		sig <- unlist(strsplit(sig,sep))
		# Return
		return(sig)
	})))
}


calcICByGenes <- function(terms, term_genes, defaultIC = 0){
	# Obtain term ICs 
	term_genes$IC <- -log(term_genes$Genes / max(term_genes$Genes))
	# Substitute Inf ICs by default IC
	term_genes$IC[which(term_genes$IC == Inf)] <- defaultIC
	# Return target genes
	return(term_genes[which(term_genes$Term %in% terms),-which(names(term_genes) =='Genes')])
}






################################################################
##                           LOAD DATA                        ##
################################################################

# Prepare containers
go_terms   <- NULL
kegg_terms <- NULL
reac_terms <- NULL
	
# Load data sets
if(!is.null(opt$go) & file.info(opt$go)$size > 1){
	# Load necessary packages
	suppressPackageStartupMessages(require(org.Hs.eg.db))
	suppressPackageStartupMessages(require(ontologyIndex))
	# Load and transform data
	go_terms <- read.table(file = opt$go, header = T, quote = "", sep = "\t", stringsAsFactors = F)
	go_terms <- unified2signalList(go_terms)
}else if(!is.null(opt$go)){
	write.table(paste(c("Term","IC","Name"),collapse = "\t"),
      file = paste(opt$output,"_go",dicts_ext,sep=""), sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}
if(!is.null(opt$kegg) & file.info(opt$kegg)$size > 1){
	# Load necessary packages
	suppressPackageStartupMessages(require(clusterProfiler))
	# Load and transform data
	kegg_terms <- read.table(file = opt$kegg, header = T, quote = "", sep = "\t", stringsAsFactors = F)
	kegg_terms <- unified2signalList(kegg_terms)
}else if(!is.null(opt$kegg)){
	write.table(paste(c("Term","IC","Name"),collapse = "\t"),
      file = paste(opt$output,"_kegg",dicts_ext,sep=""), sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}
if(!is.null(opt$reac) & file.info(opt$reac)$size > 1){
	# Load necessary packages
	suppressPackageStartupMessages(require(ReactomePA))
	# Load and transform data
	reac_terms <- read.table(file = opt$reac, header = T, quote = "", sep = "\t", stringsAsFactors = F)
	reac_terms <- unified2signalList(reac_terms)
}else if(!is.null(opt$reac)){
	write.table(paste(c("Term","IC","Name"),collapse = "\t"),
      file = paste(opt$output,"_reactome",dicts_ext,sep=""), sep = "\t", quote = F, row.names = F, append = F, col.names = F)
}



################################################################
##                           CALC ICS                         ##
################################################################

if(!is.null(go_terms)){
	# Obtain term related genes
	# genes2go <- as.list(org.Hs.egGO[mappedkeys(org.Hs.egGO)])
	genes2go <- as.list(org.Hs.egGO2ALLEGS[mappedkeys(org.Hs.egGO2ALLEGS)])
	# go_genes <- table(unlist(lapply(genes2go,function(item){return(names(item))})))
	# go_genes <- data.frame(Term = names(go_genes), Genes = as.vector(go_genes), stringsAsFactors= F)
	go_genes <- data.frame(Term = names(genes2go), Genes = unlist(lapply(genes2go,function(item){length(item)})), stringsAsFactors = F)

	# Use all genes instead observed frequencies as maximum
	if(opt$all){
		# go_genes <- rbind(data.frame(Term = "All", Genes = length(genes2go), stringsAsFactors = F), go_genes)
		go_genes <- rbind(data.frame(Term = "All", Genes = length(unique(unlist(lapply(genes2go,function(item){as.vector(item)})))), stringsAsFactors = F), go_genes)
	}

	# Add result terms with no genes associated (If it happens, you must update your DB)
	no_genes <- which(!go_terms %in% go_genes$Term)
	if(length(no_genes) > 0){
		go_genes <- rbind(go_genes,data.frame(Term = go_terms[no_genes], Genes = rep(0,length(no_genes)), stringsAsFactors = F))
	}

	# Calculate ICs
	go_ics <- calcICByGenes(go_terms,go_genes)

	# Add term names
	data(go)
	go2name <- data.frame(GO = names(go$name), Name = as.vector(go$name), stringsAsFactors = F)
	go_ics$Name <- unlist(lapply(seq(nrow(go_ics)),function(i){
		# Look for name
		indx <- which(go2name$GO == go_ics$Term[i])
		# Check
		return(ifelse(length(indx)==0,go_ics$Term[i],go2name$Name[indx[1]]))
	}))

	# Write dictionary
	write.table(go_ics,file = paste(opt$output,"_go",dicts_ext,sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}



if(!is.null(kegg_terms)){
	# Obtain term related genes
	# kegg_info <- clusterProfiler:::prepare_KEGG(clusterProfiler:::organismMapper('hsa'), "KEGG", 'kegg')
	kegg_info <- clusterProfiler:::get_data_from_KEGG_db(clusterProfiler:::organismMapper('hsa'))
	kegg2genes <- data.frame(Term = names(kegg_info$PATHID2EXTID),
							 Genes = unlist(lapply(kegg_info$PATHID2EXTID,function(item){return(length(unlist(item)))})),
							 stringsAsFactors = F)

	# Use all genes instead observed frequencies as maximum
	if(opt$all){
		kegg2genes <- rbind(data.frame(Term = "All", Genes = length(unique(lapply(kegg_info$PATHID2EXTID,function(item){return(length(unlist(item)))}))), stringsAsFactors = F), kegg2genes)
	}

	# Add result terms with no genes associated (If it happens, you must update your DB)
	no_genes <- which(!kegg_terms %in% kegg2genes$Term)
	if(length(no_genes) > 0){
		kegg2genes <- rbind(kegg2genes,data.frame(Term = kegg_terms[no_genes], Genes = rep(0,length(no_genes)), stringsAsFactors = F))
	}


	# Calculate ICs
	kegg_ics <- calcICByGenes(kegg_terms,kegg2genes)

	# Add term names
	kegg2name <- data.frame(Term = names(kegg_info$PATHID2NAME),
							Name = unlist(kegg_info$PATHID2NAME),
							stringsAsFactors = F)

	kegg_ics$Name <- unlist(lapply(seq(nrow(kegg_ics)),function(i){
		# Look for name
		indx <- which(kegg2name$Term == kegg_ics$Term[i])
		# Check
		return(ifelse(length(indx)==0,kegg_ics$Term[i],kegg2name$Name[indx[1]]))
	}))

	# Write dictionary
	write.table(kegg_ics,file = paste(opt$output,"_kegg",dicts_ext,sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}



if(!is.null(reac_terms)){
	# Obtain term related genes
	reac_info <- ReactomePA:::get_Reactome_DATA("human")
	reac2genes <- data.frame(Term = names(reac_info$PATHID2EXTID),
							 Genes = unlist(lapply(reac_info$PATHID2EXTID,function(item){return(length(unlist(item)))})),
							 stringsAsFactors = F)
	
	# Use all genes instead observed frequencies as maximum
	if(opt$all){
		reac2genes <- rbind(data.frame(Term = "All", Genes = length(unique(lapply(reac_info$PATHID2EXTID,function(item){return(length(unlist(item)))}))), stringsAsFactors = F), reac2genes)
	}

	# Add result terms with no genes associated (If it happens, you must update your DB)
	no_genes <- which(!reac_terms %in% reac2genes$Term)
	if(length(no_genes) > 0){
		reac2genes <- rbind(reac2genes,data.frame(Term = reac_terms[no_genes], Genes = rep(0,length(no_genes)), stringsAsFactors = F))
	}

	# Calculate ICs
	reac_ics <- calcICByGenes(reac_terms,reac2genes)	

	# Add term names
	reac2name <- data.frame(Term = names(reac_info$PATHID2NAME),
							Name = unlist(reac_info$PATHID2NAME),
							stringsAsFactors = F)

	reac_ics$Name <- unlist(lapply(seq(nrow(reac_ics)),function(i){
		# Look for name
		indx <- which(reac2name$Term == reac_ics$Term[i])
		# Check
		return(ifelse(length(indx)==0,reac_ics$Term[i],reac2name$Name[indx[1]]))
	}))

	# Write dictionary
	write.table(reac_ics,file = paste(opt$output,"_reactome",dicts_ext,sep=""),col.names=T, row.names=F,quote=F,sep="\t")
}

