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
  make_option(c("-p", "--path"), action="store", type="character", default = ".",
              dest="path", help="Path where target files are stored (subdirectories included)[Default: %default]"),
  make_option(c("-f", "--file"), action="store", type="character",
              dest="file", help="Path where target files are stored (subdirectories included)"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Destiny file"),
  make_option(c("-P","--pval"), action="store",type="double",default = 1,
              dest="pval", help="Pvalue to be applied. [Default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))


##############################################################################
##                                LOAD FILE                                 ##
##############################################################################

# Column names
topgo <- c("GO.ID",
           "Term",
           "Annotated",
           "Significant",
           "Expected",
           "Fisher.Test",
           "Significant.Genes",
           "ID",
           "sub_ontology",
           "tag")

clprof <- c("Model",
            "Term",
            "Source",
            "GO",
            "GO_Name",
            "Annotated",
            "Significant",
            "Expected",
            "Fisher.Test",
            "Genes")

# Create relationship between format columns sorting topgo
formats <- data.frame(ClProfiler = clprof,
                      TopGO      = topgo[c(10, 8, 9, seq(7))],
                      stringsAsFactors = F)

# Find target files
target_files <- list.files(path = opt$path, pattern = opt$file, recursive = TRUE, full.names = T)

# Check
if(length(target_files) == 0){
  stop("ERROR [topGO2clusterProfiler]: any file have been found to be translated")
}

# Reset destiny file
write.table(paste(formats$ClProfiler, collapse = "\t"),
            file = opt$output, 
            sep = "\t", quote = F, row.names = F, append = F, col.names = F)

# Per each file apply ETL process
invisible(lapply(target_files, function(f){
  # Load file
  enrichment <- read.table(file = f, sep = "\t", header = T, quote = "", stringsAsFactors = F)
  enrichment$Fisher.Test <- as.double(gsub("<","",as.character(enrichment$Fisher.Test)))
  # Check if enrichment format is correct
  if(!all(colnames(enrichment) %in% topgo)){
  stop("ERROR [topGO2clusterProfiler: enrichment file has not a correct header. Missing elements")
  }
  
  # Match source column names positions
  new_order <- match(formats$TopGO, colnames(enrichment))

  # Sort
  new_format_enrichment <- enrichment[,new_order]

  # Check
  if(length(colnames(enrichment)) > nrow(formats)){
    warning("Warning [topGO2clusterProfiler] Enrichment file given has more columns than used in formated new file. Some information can be lost..")
  }

  # Write into output file
  write.table(new_format_enrichment[which(new_format_enrichment$Fisher.Test <= opt$pval),], file = opt$output, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
}))
