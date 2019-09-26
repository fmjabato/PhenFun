#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

################################################################
##                           CONFIGURE                        ##
################################################################
suppressPackageStartupMessages(require(optparse))

option_list <- list(
  make_option(c("-t", "--tripartite"), action="store", type="character",
              dest="tripartite", help="Tripartite network"),
  make_option(c("-n", "--extendedNet"), action="store", type="character",
              dest="extendedNet", help="Extended network (HPO,Loci,Genes)"),
  make_option(c("-e", "--enrich"), action="store", type="character",
              dest="enrich", help="GO Enrichment file. 'dec' flag will be used"),
  make_option(c("-g", "--genesSep"), action="store", type="character",
              dest="genesSep", help="Genes IDs separator used on genes column"),
  make_option(c("-i", "--index"), action="store", type="character",
              dest="index", help="Index of Model, Phenotype, FunSys and GenesList columns separated by commas"),
  make_option(c("-f", "--fenrich"), action="store", type="character",
              dest="fenrich", help="Filtered Enrichment file"),
  make_option(c("-P", "--pats"), action="store", type="character", default = NULL,
              dest="pats", help="Optional file to store <Pheno, FunSys, Pat> info"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file to store ratios calculated")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
# Prepara enrichment indexes
indx <- unlist(strsplit(opt$index,","))
names(indx) <- c("Model","Pheno","FunSys","Genes")



# Load Tripartite network
tripartite <- read.table(file = opt$tripartite, sep = "\t", quote = "", stringsAsFactors = F, header = F)

# Load extended bipartite network
extended_net <- read.table(file = opt$extendedNet, sep = "\t", quote = "", stringsAsFactors = F, header = T)
extended_net <- extended_net[which(extended_net$tag == 'dec'),]

# Load System enrichments
enrich <- read.table(file = opt$enrich, sep = "\t", quote = "", stringsAsFactors = F, header = T)
# Prepare colnames
aux <- colnames(enrich)
aux[as.numeric(indx['Model'])] <- "Model"
aux[as.numeric(indx['Pheno'])] <- "Term"
aux[as.numeric(indx['FunSys'])] <- "FunSys"
aux[as.numeric(indx['Genes'])] <- "Genes"
colnames(enrich) <- aux 
enrich <- enrich[which(enrich$Model == 'dec'),]

# Load filtered enrichments
if(file.info(opt$fenrich)$size > 1){
  fenrich <- read.table(file = opt$fenrich, sep = "\t", quote = "", stringsAsFactors = F, header = T)  
}else{
  write.table(paste(c("Patient","Phenos","Genos","GWithGenes","Genes","SigPhenos","SigGenes","PercSigGenes","LinkablePhenos","Systems","SystemsUnique","PhenosPerSys","ProfMeanIC","HPOsWithoutIC","PhenosList","PhenosSigList",
      "GenesList","GenesSigList","SystemsList"),collapse = "\t"),
      file = opt$output, sep = "\t", quote = F, row.names = F, append = F, col.names = F)
  if(!is.null(opt$pats)){
    write.table(paste(c("Pheno","FunSys","Patient","EnrichGenes","PatientSigGenes","SharedGenes","GenesCoverage","PatientGenesList"),collapse = "\t"),
      file = opt$pats, sep = "\t", quote = F, row.names = F, append = F, col.names = F)
  }
  stop("There are not enrichments to handle. Generating empty output files")
}

#fenrich <- fenrich[which(fenrich$Model == opt$rdm & fenrich$PvalThr == 1e-3),]

################################################################
##                        TRANSFORM DATA                      ##
################################################################

# Prepare HPO ICs
hpo_ics <- fenrich[,c("HPO","IC")]

# Transform tripartite to real tripartite dataframe
tripartite <- as.data.frame(do.call(rbind,lapply(unique(tripartite[!grepl("HP:",tripartite[,1]),1]),function(patient){
  # Find related HPOs
  hpos <- tripartite[which(tripartite[,2] == patient),1]
  # Find related Locis
  gvrs <- tripartite[which(tripartite[,1] == patient),2]
  # Prepare container
  info <- data.frame(Patient        = patient,
                     Phenos         = length(hpos),
                     Genos          = length(gvrs),
                     Genes          = -1,
                     LocisWithGenes = -1,
                     HPOs           = "",
                     Locis          = "",
                     GenesList      = "", stringsAsFactors = F)
  # Add Phenotypes and Genotypes
  info$HPOs = list(hpos)
  info$Locis = list(gvrs)
  # Return
  return(info)
})))


# Obtain loci's genes
locigenes <- as.data.frame(do.call(rbind,lapply(unique(extended_net$Loci),function(loci){
  # Obtain linked genes
  genes <- strsplit(extended_net$Genes[which(extended_net$Loci == loci)[1]],":")
  # Prepare container
  info <- data.frame(Loci = loci, Genes = length(genes[[1]]), GenesList = "", stringsAsFactors = F)
  info$GenesList <- genes
  # Return
  return(info)
})))


# Prepare extended tripartite network
invisible(lapply(seq(nrow(tripartite)),function(i){
  # Obtain locis genes
  locis <- locigenes[which(locigenes$Loci %in% tripartite$Locis[[i]]),]
  # Find locis with genes
  withGenes <- which(locis$Genes > 0)
  # Obtain genes
  genes <- sort(unique(unlist(lapply(locis$GenesList,unlist))))
  # Update info
  tripartite$Genes[i] <<- length(genes)
  tripartite$LocisWithGenes[i] <<- length(withGenes)
  tripartite$GenesList[i] <<- list(genes)
}))


# Obtain genes used to enrich systems
hpo_sys_genes <- as.data.frame(do.call(rbind,lapply(unique(fenrich$HPO),function(hp){
  aux <- enrich[which(enrich$Term == hp),]

  info <- as.data.frame(do.call(rbind,lapply(unlist(strsplit(fenrich$Signal[which(fenrich$HPO == hp)],";")),function(go){
    # Obtain genes
    genes <- strsplit(aux$Genes[which(aux$FunSys == go)],opt$genesSep)
    # Prepare container
    info <- data.frame(HPO       = hp,
                       System    = go,
                       Genes     = length(genes[[1]]),
                       GenesList = "",
                       stringsAsFactors = F)
    # ADd genes
    info$GenesList <- genes
    # Return
    return(info)
  })))
  return(info)
})))


# Prepare HPOs enrichment info
fenrich$GenesList <- rep("",nrow(fenrich))
invisible(lapply(seq(nrow(fenrich)),function(i){
  # Genes list
  genes <- unique(unlist(lapply(unlist(strsplit(fenrich$Signal[i],";")),function(sys){
    return(hpo_sys_genes$GenesList[[which(hpo_sys_genes$HPO == fenrich$HPO[i] & hpo_sys_genes$System == sys)]])
  })))
  # Update
  fenrich$GenesList[i] <<- list(genes)
}))


# Merge to obtain final data
patients_info <- as.data.frame(do.call(rbind,lapply(seq(nrow(tripartite)),function(i){
  # Obtain significative phenotypes
  hpos <- tripartite$HPOs[[i]][which(tripartite$HPOs[[i]] %in% fenrich$HPO)]
  # Obtain significative genes
  sig_genes <- tripartite$GenesList[[i]][which(tripartite$GenesList[[i]] %in% unlist(lapply(hpos,function(hp){return(fenrich$GenesList[[which(fenrich$HPO == hp)]])})))]
  # Per each phenotype, find linkable enriched systems <Pheno, FSystem>
  phenos_signals <- as.data.frame(do.call(rbind,lapply(hpos,function(hp){
    # Obtain systems enriched
    sys <- hpo_sys_genes[which(hpo_sys_genes$HPO == hp),]
    # Find linkable systems
    indx <- which(unlist(lapply(seq(nrow(sys)),function(j){return(any(sys$GenesList[[j]] %in% sig_genes))})))
    #Check
    if(length(indx) > 0){
      return(sys[indx,1:2])
    }else{
      return(data.frame())
    }
  })))

  # Obtain Phenotype profile IC
  meanic <- unlist(lapply(hpos,function(hp){
    # Find
    indx <- which(hpo_ics$HPO == hp)
    if(length(indx) == 0){
      return(NA)
    }
    return(hpo_ics$IC[which(hpo_ics$HPO == hp)])
  }))
  if(any(is.na(meanic))){
    nas <- which(is.na(meanic))
    meanic <- meanic[-nas]
    nas <- length(nas)
  }else{
    nas <- 0
  }
  if(length(meanic) > 0){
    meanic <- mean(meanic)
  }else{
    meanic <- 0
  }


  # Generate dataframe
  info <- data.frame(Patient        = tripartite$Patient[i],
                     Phenos         = tripartite$Phenos[i],
                     Genos          = tripartite$Genos[i],
                     GWithGenes     = tripartite$LocisWithGenes[i],
                     Genes          = tripartite$Genes[i],
                     SigPhenos      = length(hpos),
                     SigGenes       = length(sig_genes),
                     PercSigGenes   = length(sig_genes)/tripartite$Genes[i],
                     LinkablePhenos = length(unique(phenos_signals$HPO)),
                     Systems        = length(phenos_signals$System),
                     SystemsUnique  = length(unique(phenos_signals$System)),
                     PhenosPerSys   = mean(table(phenos_signals$System)),
                     ProfMeanIC     = meanic,
                     HPOsWithoutIC  = nas,
                     PhenosList     = paste(tripartite$HPOs[[i]], collapse = ";"),
                     PhenosSigList  = paste(hpos, collapse = ";"),
                     GenesList      = paste(tripartite$GenesList[[i]], collapse = ";"),
                     GenesSigList   = paste(sig_genes, collapse = ";"),
                     SystemsList    = paste(unique(phenos_signals$System),collapse=";"),
                     stringsAsFactors = F)
  # return
  return(info)
})))

patients_info <- patients_info[order(patients_info$Patient),]


################################################################
##                          WRITE DATA                        ##
################################################################
write.table(patients_info, file = opt$output, col.names = T, row.names = F, sep = "\t", quote = F)

################################################################
##                         SPECIAL CASE                       ##
################################################################
if(!is.null(opt$pats)){
  # Transform data
  triplets <- as.data.frame(do.call(rbind,lapply(seq(nrow(patients_info)),function(i){
    # Check special case
    if(patients_info$SigPhenos[i] == 0 | patients_info$Systems[i] == 0){
      return(data.frame())
    }
    # Prepare patient info
    genes_pat <- unlist(strsplit(patients_info$GenesSigList[i],";"))
    # Per each patient significative phenotype
    inf <- as.data.frame(do.call(rbind,lapply(unlist(strsplit(patients_info$PhenosSigList[i],";")),function(pheno){
      # Find related FunSys
      sys <- hpo_sys_genes[which(hpo_sys_genes$HPO == pheno),]
      # Obtain enrichment genes
      aux_inf <- as.data.frame(do.call(rbind,lapply(seq(nrow(sys)),function(j){
        # Calc coverage
        genes_rel <- unlist(sys$GenesList[j])
        shared    <- intersect(genes_pat,genes_rel)
        if(length(shared) == 0){
          cov <- 0
        }else{
          cov <- length(shared) / length(genes_rel)
        }
        # Obtain info
        return(data.frame(Pheno            = pheno,
                          FunSys           = sys$System[j],
                          Patient          = patients_info$Patient[i],
                          EnrichGenes      = length(genes_rel),
                          PatientSigGenes  = length(genes_pat),
                          SharedGenes      = length(shared),
                          GenesCoverage    = cov,
                          PatientGenesList = ifelse(length(shared)>0,paste(shared, collapse = ";"),NA),
                          # Test1 = paste(genes_rel, collapse = ";"),
                          # Test2 = paste(genes_pat, collapse = ";"),
                          stringsAsFactors = F))
      })))
      return(aux_inf)
    })))
    return(inf)
  })))

  triplets <- triplets[which(triplets$GenesCoverage > 0),]
  if(nrow(triplets) > 0)
    triplets <- triplets[order(triplets$Pheno,triplets$FunSys),]

  # Write info
  write.table(triplets, file = opt$pats, col.names = T, row.names = F, sep = "\t", quote = F)  
}