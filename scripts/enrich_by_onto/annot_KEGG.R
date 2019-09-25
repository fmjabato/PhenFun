#' @author Fernando Moreno Jabato
#' @description Functions to handle KEGG enrichments based on a set of genes
#' @method annot_KEGG : used to enrich using KEGG based on a set of genes
#' @method annot_sets_KEGG : used to enrich a pull of sets of genes using KEGG  


annot_KEGG <- function(genes, KEGG_DATA = NULL, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                       keyType = "kegg"){
#' @description Optimization of "clusterProfiler" package KEGG annotation method
#' avoiding computational expensive repetived tasks
#' @param genes a vector (or set of vectors) which includes a set of genes
#' @param KEGG_DATA KEGG database already loaded in clusterProfiler specific structure
#' @param organism identifier. Default: hsa
#' @param pAdjustMethod p-value adjust method {holm, hochberg, hommel, bonferroni, BH, BY, fdr, none}. Default: BH
#' @param pvalueCutoff p-value threshold. Default: 0.05
#' @param qvalueCutoff q-value threshold. Default: 0.2
#' @param keyType input format {kegg, ncbi-geneid, ncbi-proteinid, uniprot}. Default: kegg
#' @return an enrichResult instance with KEGG enrichment or NULL if no enrichments i s possible.
#' @throws exceptions if any error occurs.
#'  
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom clusterProfiler package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}

  # CHECK INPUTS
  if(is.null(genes)){ # Genes
    stop("[annot_KEGG] Given genes set is NULL")
  }else if(is.list(genes)){
    genes <- unlist(genes)
    if(!is.vector(genes)){
      stop("[annot_KEGG] Genes is not, or is not convertible to, a vector")
    }
  }else if(!is.vector(genes)){
    stop("[annot_KEGG] Genes is not, or is not convertible to, a vector")
  }
  
  if(!is.null(KEGG_DATA)){ # KEGG_DATA
    if(!is.environment(KEGG_DATA)){
      stop("[annot_KEGG] KEGG_DATA set is not an environment")
    }
  }
  
  if(is.null(organism)){ # Organism
    stop("[annot_KEGG] Given Organism code is NULL")
  }else if(!is.character(organism)){
    stop("[annot_KEGG] Given organism is not a character string")
  }
  
  if(is.null(pAdjustMethod)){ # pAdjustMethod
    stop("[annot_KEGG] Given pAdjustMethod code is NULL")
  }else if(!is.character(pAdjustMethod)){
    stop("[annot_KEGG] Given pAdjustMethod is not a character string")
  }
  
  if(is.null(pvalueCutoff)){ # pvalueCutoff
    stop("[annot_KEGG] Given pvalueCutoff is NULL")
  }else if(!is.numeric(pvalueCutoff)){
    stop("[annot_KEGG] Given pvalueCutoff is not a numeric value")
  }else if(pvalueCutoff < 0){
    warning("[annot_KEGG] Given pvalueCutoff is a negative number. Zero value will be used")
    pvalueCutoff <- 0
  }
  
  if(is.null(qvalueCutoff)){ # qvalueCutoff
    stop("[annot_KEGG] Given qvalueCutoff is NULL")
  }else if(!is.numeric(qvalueCutoff)){
    stop("[annot_KEGG] Given qvalueCutoff is not a numeric value")
  }else if(qvalueCutoff < 0){
    warning("[annot_KEGG] Given qvalueCutoff is a negative number. Zero value will be used")
    qvalueCutoff <- 0
  }
  
  if(is.null(keyType)){ # keyType
    stop("[annot_KEGG] Given KeyType code is NULL.")
  }else if(!is.character(keyType)){
    stop("[annot_KEGG] Given KeyType is not a character string")
  }
  
  
  # Load necessary packages
  require(clusterProfiler)
  
  # Obtain wanted function
  enrichKEGG <- clusterProfiler::enrichKEGG
  
  # Check if KEGG_DATA environment has been given
  if(!is.null(KEGG_DATA)){ # Use special KEGG Annotation function
    # Find not necessary task into code
    line_to_remove <- grep("KEGG_DATA *<-",body(enrichKEGG))
    # Check
    if(length(line_to_remove) == 0){ # Warning, task not found
      warning("annot_KEGG: Can not find GO annot task to be removed. Regular version will be used.")
    }else{ # Remove task from code
      body(enrichKEGG)[[line_to_remove]] <- substitute(KEGG_DATA <- parent.frame()$KEGG_DATA)      
    }
  }# ELSE: Use regular KEGG Annotation function
  
  # Enrich
  enrich <-  enrichKEGG(gene          = genes,
                        organism      = organism,
                        keyType       = keyType,
                        pvalueCutoff  = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        qvalueCutoff  = qvalueCutoff)
  # Return results
  return(enrich)
}



annot_sets_KEGG <- function(genes_sets, KEGG_DATA = NULL, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                            keyType = "kegg", verbose = F, split = NULL, set_names = NULL){
#' @description Optimization of "clusterProfiler" package KEGG annotation method
#' avoiding computational expensive repetived tasks and managing several genes sets
#' @param genes_sets a list with genes sets to be enriched using GO
#' @param KEGG_DATA KEGG database already loaded in clusterProfiler specific structure
#' @param organism identifier. Default: hsa
#' @param pAdjustMethod p-value adjust method {holm, hochberg, hommel, bonferroni, BH, BY, fdr, none}. Default: BH
#' @param pvalueCutoff p-value threshold. Default: 0.05
#' @param qvalueCutoff q-value threshold. Default: 0.2
#' @param keyType input format {kegg, ncbi-geneid, ncbi-proteinid, uniprot}. Default: kegg
#' @param verbose activate verbose mode. Default: False
#' @param split a string used to split genes stored into genes_sets. If is NULL, split process is avoided. Default: NULL
#' @param set_names names used instead set IDs (index)
#' @return a dataframe with all enrichments generated (without filtering) or NULL if any set generates an enrichment
#'  
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom clusterProfiler package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}
#' @seealso annot_KEGG

  # CHECK INPUTS (only check NULL for inputs that will be checked at annot_KEGG)
  if(is.null(genes_sets)){ # genes_sets
    stop("[annot_sets_KEGG] Given genes sets is NULL")
  }
  
  # Load necessary packages
  require(clusterProfiler)
  if(verbose){
    require(pbapply)
  }
  
  # Check list of genes
  if(!is.null(split)){
    genes_sets <- lapply(genes_sets, function(set){unlist(strsplit(set,split))})
  }
  
  # Check if KEGG_DATA environment has been given
  if(is.null(KEGG_DATA)){ # Load KEGG set to avoid repetitive load
    # KEGG_DATA <- clusterProfiler:::prepare_KEGG(clusterProfiler:::organismMapper(organism), "KEGG", keyType)
    KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db(clusterProfiler:::organismMapper(organism))
  }# ELSE: Use given KEGG set
  
  
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
      enrich <- annot_KEGG(genes         = genes_sets[[i]],
                           KEGG_DATA     = KEGG_DATA,
                           organism      = organism,
                           pAdjustMethod = pAdjustMethod,
                           pvalueCutoff  = pvalueCutoff,
                           qvalueCutoff  = qvalueCutoff,
                           keyType       = keyType)
      
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
      enrich <- annot_KEGG(genes         = genes_sets[[i]],
                           KEGG_DATA     = KEGG_DATA,
                           organism      = organism,
                           pAdjustMethod = pAdjustMethod,
                           pvalueCutoff  = pvalueCutoff,
                           qvalueCutoff  = qvalueCutoff,
                           keyType       = keyType)
      
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
  
  # Unlist
  
  
  # Check names
  if(!is.null(set_names)){
    invisible(lapply(seq_along(genes_sets), function(i){
      # Find
      indexes <- which(enrichment$Set == i)
      # Upload
      if(length(indexes) > 0){
        enrichment$Set[indexes] <<- set_names[[i]]
      }
    }))
  }
  
  # Return generated enrichment table
  return(enrichment)
}
