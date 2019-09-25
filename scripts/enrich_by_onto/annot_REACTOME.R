##' @description This file contains several functions related to REACTOME annotations with a set 
##' of genes or a pull of sets. This functions are based on ReactomePA 
##' package but implementing functionalitie to upgrade recourse consumption
##' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
##' @importFrom ReactomePA package
##' @seealso \link{https://bioconductor.org/packages/release/bioc/html/ReactomePA.html}
##' @method annot_REACTOME
##' @method annto_sets_REACTOME


##' Optimization of "ReactomePA" package REACTOME annotation method
##' avoiding computational expensive repetived tasks
##' @param genes a vector which includes a set of genes
##' @param Reactome_DATA database already loaded in ReactomePA specific structure
##' @param organism source
##' @param pAdjustMethod p-value adjust method {holm, hochberg, hommel, bonferroni, BH, BY, fdr, none}. Default: BH
##' @param pvalueCutoff p-value threshold. Default: 0.05
##' @param qvalueCutoff q-value threshold. Default: 0.2
##' @param readable whether mapping gene ID to gene Name. Default: TRUE
##' @return an enrichResult instance with REACTOME enrichment
##'  
##' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
##' @importFrom ReactomePA package
##' @seealso \link{https://bioconductor.org/packages/release/bioc/html/ReactomePA.html}
annot_REACTOME <- function(genes, Reactome_DATA = NULL, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE){
  # CHECK INPUTS (not implemented yet)
  
  # Load necessary packages
  require(ReactomePA)
  
  # Obtain wanted function
  enrichPathway <- ReactomePA::enrichPathway
  
  # Check if Reactome_DATA environment has been given
  if(!is.null(Reactome_DATA)){ # Use special KEGG Annotation function
    # Find not necessary task into code
    line_to_remove <- grep("Reactome_DATA *<-",body(enrichPathway))
    # Check
    if(length(line_to_remove) == 0){ # Warning, task not found
      warning("annot_REACTOME: Can not find GO annot task to be removed. Regular version will be used.")
    }else{ # Remove task from code
      body(enrichPathway)[[line_to_remove]] <- substitute(Reactome_DATA <- parent.frame()$Reactome_DATA)      
    }
  }# ELSE: Use regular REACTOME Annotation function
  
  # Enrich
  enrich <- enrichPathway(gene          = genes,
                          organism      = organism,
                          pvalueCutoff  = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          qvalueCutoff  = qvalueCutoff,
                          readable      = readable)
  
  # Return results
  return(enrich)
}


##' @description Optimization of "ReactomePA" package REACTOME annotation method
##' avoiding computational expensive repetived tasks and managing several genes sets
##' @param genes_sets a list with genes sets to be enriched using REACTOME
##' @param Reactome_DATA Reactome database already loaded in ReactomePA specific structure
##' @param organism source of genes
##' @param pAdjustMethod p-value adjust method {holm, hochberg, hommel, bonferroni, BH, BY, fdr, none}. Default: BH
##' @param pvalueCutoff p-value threshold. Default: 0.05
##' @param qvalueCutoff q-value threshold. Default: 0.2
##' @param readable whether mapping gene ID to gene Name. Default: TRUE
##' @param verbose activate verbose mode. Default: False
##' @param split a string used to split genes stored into genes_sets. If is NULL, split process is avoided. Default: NULL
##' @param set_names names used instead set IDs (index)
##' @return a dataframe with all enrichments generated (without filtering)
##'  
##' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
##' @importFrom ReactomePA package
##' @seealso \link{https://bioconductor.org/packages/release/bioc/html/ReactomePA.html}
##' @seealso annot_REACTOME
annot_sets_REACTOME <- function(genes_sets, Reactome_DATA = NULL, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE, 
                                verbose = F, split = NULL, set_names = NULL){
  # CHECK INPUTS (not implemented yet)
  
  # Load necessary packages
  require(ReactomePA)
  if(verbose){
    require(pbapply)
  }
  
  # Check list of genes
  if(!is.null(split)){
    genes_sets <- lapply(genes_sets, function(set){unlist(strsplit(set,split))})
  }
  
  # Check if Reactome_DATA environment has been given
  if(is.null(Reactome_DATA)){ # Load REACTOME set to avoid repetitive load
    Reactome_DATA <- ReactomePA:::get_Reactome_DATA(organism)
  }# ELSE: Use given REACTOME set
  
  
  # Prepare necessary functions
  options(stringsAsFactors = F)
  
  # Enrich
  if(verbose){
    enrichment <- as.data.frame(do.call(rbind,pblapply(seq_along(genes_sets), function(i){
      # Check
      if(length(genes_sets[[i]]) == 0){
        return(NULL)
      }
      # Enrich
      enrich <- annot_REACTOME(genes         = genes_sets[[i]],
                               Reactome_DATA = Reactome_DATA,
                               organism      = organism, 
                               pAdjustMethod = pAdjustMethod,
                               pvalueCutoff  = pvalueCutoff,
                               qvalueCutoff  = qvalueCutoff,
                               readable      = readable)
      # Check and prepare result
      if(is.null(enrich)){
        return(NULL)
      }else if(nrow(enrich@result) == 0){
        return(NULL)
      }else{
        result <- cbind(Set = rep(i,nrow(enrich@result)), enrich@result)
        return(result)
      }
    })))  
  }else{
    enrichment <- as.data.frame(do.call(rbind,lapply(seq_along(genes_sets), function(i){
      # Check
      if(length(genes_sets[[i]]) == 0){
        return(NULL)
      }
      # Enrich
      enrich <- annot_REACTOME(genes         = genes_sets[[i]],
                               Reactome_DATA = Reactome_DATA,
                               organism      = organism, 
                               pAdjustMethod = pAdjustMethod,
                               pvalueCutoff  = pvalueCutoff,
                               qvalueCutoff  = qvalueCutoff,
                               readable      = readable)
      # Check and prepare result
      if(is.null(enrich)){
        return(NULL)
      }else if(nrow(enrich@result) == 0){
        return(NULL)
      }else{
        result <- cbind(Set = rep(i,nrow(enrich@result)), enrich@result)
        return(result)
      }
    })))  
  }
  
  # Check names
  if(!is.null(set_names)){
    invisible(lapply(seq_along(genes_sets), function(i){
      # Find
      indexes <- which(enrichment$Set == i)
      # Upload
      if(length(indexes) > 0){
        enrichment$Set[indexes] <<- set_names[i]
      }
    }))
  }
  
  # Return generated enrichment table
  return(enrichment)
}
