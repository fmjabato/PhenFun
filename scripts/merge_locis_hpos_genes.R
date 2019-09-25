#! /usr/bin/env Rscript

#' @description Merge specific format files
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse


###########################################################################
##                          CONFIGURE PROGRAM                            ##
###########################################################################
options(stringsAsFactors = F)

# Load necessary packages
suppressPackageStartupMessages(require(optparse))        # Parse script inputs

# Prepare input commands
option_list <- list(
  make_option(c("-n", "--hpoloci"), action="store", type="character",
              help="File with network HPO - Loci"),
  make_option(c("-l","--locigenes"), action="store", type = "character",
              help="File with network Loci - Genes"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Destiny file"),
  make_option(c("-r","--reconstruct"), action = "store_true", type = "logical", default = F,
              dest = "reconstruct", help = "Flag which activate reconstruction mode. Which substitute HPO-Gene output by HPO-Loci-Genes output (separated by ':')"),
  make_option(c("-t","--tag"), action = "store", type = "character", dest = "tag",default = NULL, help = "Tag to be added")
)

opt <- parse_args(OptionParser(option_list=option_list))

###########################################################################
##                          LOAD & TRANSFORM                             ##
###########################################################################

# Load both files
locis_hpo   <- read.table(file = opt$hpoloci, header = F, stringsAsFactors = F, sep = "\t")
locis_genes <- read.table(file = opt$locigenes, header = F, stringsAsFactors = F, sep = "\t")

if(!opt$reconstruct){
   # Per each loci, generate link HPO - Gene
   hpo_genes <- as.data.frame(do.call(rbind,lapply(unique(locis_hpo[,1]), function(loci){
     # Find locis index into both sets
     hpo_indx <- which(locis_hpo[,1] == loci)
     gen_indx <- which(locis_genes[,1] == loci)
     # Concatenate HPO-Gene
     info <- data.frame(HPO  = unlist(lapply(hpo_indx, function(i){return(rep(locis_hpo[i,2],length(gen_indx)))})),
                        Gene = rep(locis_genes[gen_indx,2],length(hpo_indx)),
                        stringsAsFactors = F)
     # return
     return(info)
   })))

   # Sort by HPO
   hpo_genes <- hpo_genes[order(hpo_genes[,1],hpo_genes[,2]),]

   # Write file into output file
   write.table(hpo_genes, file = opt$output, sep = "\t", quote = F, col.names = F, row.names = F)
}else{ # Reconstruct mode
   # Create container
   hpo_loci_genes <- data.frame(HPO = locis_hpo[,2], Loci = locis_hpo[,1], Genes = character(nrow(locis_hpo)), stringsAsFactors = F)
   if(!is.null(opt$tag)){
	hpo_loci_genes <- cbind(tag = rep(opt$tag,nrow(hpo_loci_genes)), hpo_loci_genes)
   }
   # Per each Loci, add genes to DF
   invisible(lapply(unique(hpo_loci_genes$Loci),function(loci){
      # Find loci into each set
      indx_hpo <- which(hpo_loci_genes$Loci == loci)
      indx_gen <- which(locis_genes[,1] == loci)
      # Collapse genes
      genes <- paste(locis_genes[indx_gen,2], collapse = ":")
      # Add genes to each loci row
      hpo_loci_genes$Genes[indx_hpo] <<- genes
   }))
   # Write extended file
   write.table(hpo_loci_genes, file = opt$output, sep = "\t", quote = F, col.names = T, row.names = F)
}
