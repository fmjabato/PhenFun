#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(org.Hs.eg.db))

option_list <- list(
  make_option(c("-e", "--enrich"), action="store", type="character",
              dest="enrich", help="GO Enrichment file. 'dec' flag will be used"),
  make_option(c("-g", "--genesSep"), action="store", type="character",
              dest="genesSep", help="Genes IDs separator used on genes column"),
  make_option(c("-i", "--index"), action="store", type="numeric",
              dest="index", help="Index of GenesList column"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file. Same input file data but with gene names translated")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                      TRANSFORM DATA                        ##
################################################################
# Load data
dataset <- read.table(file = opt$enrich, sep = "\t", quote = "", header = T, stringsAsFactors = F)

# Prepare genes dictionary
genes_all <- mapIds(org.Hs.eg.db, names(as.list(org.Hs.egGO[mappedkeys(org.Hs.egGO)])),'SYMBOL','ENTREZID')
genes_all <- data.frame(ENTREZ = names(genes_all), SYMBOL = as.vector(genes_all), stringsAsFactors = F)

# Per each column, translate
j <- opt$index
invisible(lapply(seq(nrow(dataset)),function(i){
  # Obtain genes
  genes <- unlist(strsplit(dataset[i,j],opt$genesSep))
  # Translate
  genes <- genes_all$ENTREZ[which(genes_all$SYMBOL %in% genes)]
  # Update
  dataset[i,j] <<- paste(genes,collapse = opt$genesSep)
}))

# Write output
write.table(dataset, file = opt$output, sep = "\t", col.names = T, row.names = F, quote = F)