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
  make_option(c("-n", "--network"), action="store", type="character",
              dest="flinks", help="Network file"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Destiny file"),
  make_option(c("-e","--enrich"), action="store", type = "character", default = "go",
              dest="enrich", help="Enrichment type. Values allowed: 'go', 'kegg', 'reactome'. [Default = %default]"),
  make_option(c("-p","--pval"), action="store", type = "double", default = 0.05,
              dest="pval", help="P-val cutoff"),
  make_option(c("-a","--padjust"), action="store", type = "character", default = "BH",
              dest="padjust", help="Pvalue adjust method. Allowed: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' [Default = %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              dest="verbose",help="Activate verbose mode"),
  make_option(c("-f", "--force"), action = "store_true", default = FALSE, help = "Force generate new enrichment instead already stored without filter"),
  make_option(c("-t","--tag"), action="store", type = "character", default = NULL,
              dest="tag", help="Tag to be added as a column. Optional")
)

opt <- parse_args(OptionParser(option_list=option_list))

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                       error=function(e) # works when using R CMD
                         normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
bname <- dirname(full.fpath)

# Import functions
if(is.null(bname)){
  source_files <- list.files(path = paste("enrich_by_onto",.Platform$file.sep,sep=""), pattern = "\\.R$", full.names = T)
}else{
  source_files <- list.files(path = paste(bname,"enrich_by_onto",sep=.Platform$file.sep), pattern = "\\.R$",full.names = T)
}
invisible(lapply(source_files, source))
to_remove <- c("to_remove","source_files","script.name","script.path","initial.options")
rm(list = ls()[which(ls() %in% to_remove)])

# Obtain execution path
exec_path <- getwd()

# Configure pb package
if(opt$verbose){
  pboptions(type="timer")
}


# Check inputs
if(opt$enrich != "kegg" & opt$enrich != "go" & opt$enrich != "reactome"){
  stop("Only GO, KEGG and REACTOME enrichments are available now.")
}


##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################
# Check modes
mode_GO       <- F
mode_REACTOME <- F
mode_KEGG     <- F

if(opt$enrich == "kegg"){
  # Activate mode
  mode_KEGG <- T
  
  # Load necessary packages
  require(clusterProfiler)
  
  # Necessary file
  kegg_file <- file.path(exec_path,"KEGG_Enrich.RData")
  
  # Check force mode
  if(opt$force & file.exists(kegg_file)){
   file.remove(kegg_file)
  }
  
  # Check if already exists files to improve resource consumption
  if(file.exists(kegg_file)){
    load(kegg_file)
  }
}else if(opt$enrich == "go"){
  # Activate mode
  mode_GO <- T
  
  # Load necessary packages
  require(clusterProfiler)
  require(org.Hs.eg.db)
  
  # Necessary file
  go_file <- file.path(exec_path,"GO_Enrich.RData")
  
  # Check force mode
  if(opt$force & file.exists(go_file)){
   file.remove(go_file)
  }
  # Check if already exists files to improve resource consumption
  if(file.exists(go_file)){
    load(go_file)
  }
}else if(opt$enrich == "reactome"){
  mode_REACTOME <- T

  # Load necessary packages
  require(ReactomePA)
  
  # Necessary file
  reac_file <- file.path(exec_path, "REACT_Enrich.RData")

  # Check force mode
  if(opt$force & file.exists(reac_file)){
   file.remove(reac_file)
  }
  # Check if already exists files to improve resource consumption
  if(file.exists(reac_file)){
    load(reac_file)
  }
}

##############################################################################
##                          LOAD NETWORK DATA                               ##
##############################################################################
# Load network table
if(!exists("enrichment")){
	network <- read.table(file   = opt$flinks,
	                      header = F,
	                      sep    = "\t",
	                      stringsAsFactors = F)

	names(network) <- c("ID","Gene")

	# Collapse by Key
	genes_sep <- ":"
	tnetwork <- as.data.frame(do.call(rbind,lapply(unique(network$ID),function(id){
	  # Find tuples
	  indx <- which(network$ID == id)
	  # Collapse tuples
	  genes <- paste(network$Gene[indx], collapse = genes_sep)
	  # Return info
	  return(list(ID    = id,
	              Genes = genes))
	})))
	# Unlist columns
	for(i in seq(ncol(tnetwork))){
	  tnetwork[,i] <- unlist(tnetwork[,i])
	  rm(i)
	}

	##############################################################################
	##                                 ENRICH                                   ##
	##############################################################################
	# Enrich by key
	if(mode_KEGG){
	  # Enrich
	  enrichment <- annot_sets_KEGG(genes_sets    = tnetwork$Genes,
	                                KEGG_DATA     = NULL,
	                                pvalueCutoff  = opt$pval,
	                                pAdjustMethod = opt$padjust,
	                                verbose       = opt$verbose,
	                                set_names     = tnetwork$ID,
	                                split         = genes_sep) 
	}

	if(mode_GO){
	  # Enrich
	  enrichment <- annot_sets_GO(genes_sets    = tnetwork$Genes,
	                              OrgDb         = org.Hs.eg.db,
	                              GO_DATA       = NULL,
	                              pvalueCutoff  = opt$pval,
	                              pAdjustMethod = opt$padjust,
	                              verbose       = opt$verbose,
	                              split         = genes_sep,
	                              set_names     = tnetwork$ID,
	                              ont           = "BP")
	}
	if(mode_REACTOME){
	  enrichment <-annot_sets_REACTOME(genes_sets    = tnetwork$Genes,
	  	                               Reactome_DATA = NULL,
	  	                               organism      = "human",
	  	                               pAdjustMethod = opt$padjust,
	  	                               pvalueCutoff  = opt$pval,
	  	                               readable      = TRUE, 
                                       verbose       = opt$verbose, 
                                       split         = genes_sep, 
                                       set_names     = tnetwork$ID)
	}
} # End IF exists enrichment

##############################################################################
##                                  STORE                                   ##
##############################################################################
# Write enrichment without filtering
if(mode_GO){
 if(!file.exists(go_file)){
  save(enrichment, file = go_file)
 }
}
if(mode_KEGG){
 if(!file.exists(kegg_file)){
  save(enrichment, file = kegg_file)
 }
}
if(mode_REACTOME){
 if(!file.exists(reac_file)){
  save(enrichment, file = reac_file)
 }
}

if(!is.null(opt$tag)){
  enrichment <- cbind(Model = rep(opt$tag,nrow(enrichment)), enrichment)
}

# Write table with header
#write.table(x         = enrichment[which(enrichment$pvalue <= opt$pval),],
write.table(x         = enrichment[which(enrichment$p.adjust <= opt$pval),],
            file      = opt$output,
            quote     = F,
            sep       = "\t",
            col.names = T,
            row.names = F)
