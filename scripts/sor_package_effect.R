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
              dest="enrich", help="Meta-Enrichment table with necessary info"),
  make_option(c("-n", "--net"), action="store", type="character",
              dest="net", help="Networks file with model, phenotypes and SORs linked"),
  make_option(c("-t", "--threshold"), action="store", type="double",
              dest="threshold", help="Threshold to select from meta-enrichments data"),
  make_option(c("-p", "--pack"), action="store", type="numeric", default = 2,
              dest="pack", help="Value used to determine if a SORs is showing package effect. enrich_genes(SOR) >= VALUE"),
  make_option(c("-v", "--verbose"), action="store_true", type="logical", default = FALSE,
              dest="verbose", help="Activate verbose mode"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store ratios calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))

################################################################
##                           LOAD DATA                        ##
################################################################

# Load network
net <- read.table(file = opt$net, sep = "\t", quote = "", header = T, stringsAsFactors = F)

# Load enrichments
meta_enrich <- read.table(file = opt$enrich, sep = "\t", quote = "", header = T, stringsAsFactors = F)

# Select only meta_enrich wanted
meta_enrich <- meta_enrich[which(meta_enrich$Pval == opt$threshold),]

# Check
if(nrow(meta_enrich) <= 0){
  stop(paste("Meta-Enrichment set is empty after apply threshold:",opt$threshold))
}



################################################################
##                   CALCULATE PACKAGE EFFECT                 ##
################################################################
if(opt$verbose){
  require(pbapply)
  pboptions(type="timer")
  appFun <- pblapply
}else{
  appFun <- lapply
}

models <- unique(net$tag)
if('hpo' %in% models){
  models <- models[-which(models == 'hpo')]
}


# WRITE OUTPUT HEADER
# write.table(paste(c("Model","Pheno","SOR","Pheno_E_Genes","EnrichGenes","Package_Effect"),collapse = "\t"),
write.table(paste(c("Model","Pheno","SOR","Pheno_E_Genes","EnrichGenes","Package_Effect","GenesList"),collapse = "\t"),
      file = opt$output, sep = "\t", quote = F, row.names = F, append = F, col.names = F)



# Calculata packages
invisible(appFun(models,function(model){
  # Obtain model subsets
  sub_net <- net[which(net$tag == model),]
  sub_enr <- meta_enrich[which(meta_enrich$Model == model),]

  sor_genes <- lapply(unique(sub_net$Loci),function(sor){
    genes <- unlist(strsplit(sub_net$Genes[which(sub_net$Loci == sor)[1]],":"))
  })
  names(sor_genes) <- unique(sub_net$Loci)

  pheno_genes <- lapply(unique(sub_net$HPO),function(pheno){
    genes <- unlist(strsplit(sub_enr$GenesList[which(sub_enr$HPO == pheno)],":"))
  })
  names(pheno_genes) <- unique(sub_net$HPO)


  info <- as.data.frame(do.call(rbind,lapply(unique(sub_net$HPO),function(pheno){
    # Obtain phenotype linked info
    sors <- unique(sub_net$Loci[which(sub_net$HPO == pheno)])
    # Per each SOR obtain number of genes given to final set
    pheno_info <- as.data.frame(do.call(rbind,lapply(sors,function(sor){
      egenes <- length(intersect(sor_genes[[sor]], pheno_genes[[pheno]]))
      glist <- intersect(sor_genes[[sor]], pheno_genes[[pheno]])
      if(length(glist) <= 0){
        glist <- ""
      }else{
        glist <- paste(glist,collapse=":")
      }
      return(data.frame(Model          = model,
                        Pheno          = pheno,
                        SOR            = sor,
                        Pheno_E_Genes  = length(pheno_genes[[pheno]]),
                        Enrich_Genes   = egenes,
                        Package_Effect = egenes >= opt$pack,
                        EGenes_List    = glist,
                        stringsAsFactors = FALSE))

    })))
    # Return info
    return(pheno_info)
  })))


  # Write info
  write.table(info, file = opt$output, sep = "\t", append = T, quote = F, col.names = F, row.names = F)
}))
