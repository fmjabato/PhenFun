#' @description This file contains several functions related to GO annotations with a set 
#' of genes or a pull of sets. This functions are based on clusterProfiler 
#' package but implementing functionalitie to upgrade recourse consumption
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom clusterProfiler package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}



annot_GO <- function(genes, GO_DATA = NULL, OrgDb, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = T, keyType = "ENTREZID"){
#' @description Optimization of "clusterProfiler" package GO annotation method
#' avoiding computational expensive repetived tasks
#' @param genes a vector which includes a set of genes
#' @param GO_DATA GO database already loaded in clusterProfiler specific structure
#' @param OrgDB organisms genes database
#' @param ont GO subontology to be used {MF,BP,CC}. Default: BP
#' @param pAdjustMethod p-value adjust method {holm, hochberg, hommel, bonferroni, BH, BY, fdr, none}. Default: BH
#' @param pvalueCutoff p-value threshold. Default: 0.05
#' @param qvalueCutoff q-value threshold. Default: 0.2
#' @param readable whether mapping gene ID to gene Name. Default: TRUE
#' @param keyType format of input key genes. Default: ENTREZID
#' @return an enrichResult instance with GO enrichment
#'  
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom clusterProfiler package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}

  # Load necessary packages
  require(clusterProfiler)
  
  # Obtain wanted function
  enrichGO <- clusterProfiler::enrichGO
  
  # Check if GO_DATA environment has been given
  if(!is.null(GO_DATA)){ # Use special GO Annotation function
    # Find not necessary task into code
    line_to_remove <- grep("GO_DATA *<-",body(enrichGO))
    # Check
    if(length(line_to_remove) == 0){ # Warning, task not found
      warning("annot_GO: Can not find GO annot task to be removed. Regular version will be used.")
    }else{ # Remove task from code
      body(enrichGO)[[line_to_remove]] <- substitute(GO_DATA <- parent.frame()$GO_DATA)      
    }
  }# ELSE: Use regular GO Annotation function

  # Enrich
  enrich <- enrichGO(gene          = genes,
                     OrgDb         = OrgDb,
                     ont           = ont,
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff  = pvalueCutoff,
                     qvalueCutoff  = qvalueCutoff,
                     readable      = readable,
                     keyType       = keyType)
 
  # Return results
  return(enrich)
}



annot_sets_GO <- function(genes_sets, GO_DATA = NULL, OrgDb, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = T, 
                         keyType = "ENTREZID", verbose = F, split = NULL, set_names = NULL){
#' @description Optimization of "clusterProfiler" package GO annotation method
#' avoiding computational expensive repetived tasks and managing several genes sets
#' @param genes_sets a list with genes sets to be enriched using GO
#' @param GO_DATA GO database already loaded in clusterProfiler specific structure
#' @param OrgDB organisms genes database
#' @param ont GO subontology to be used {MF,BP,CC}. Default: BP
#' @param pAdjustMethod p-value adjust method {holm, hochberg, hommel, bonferroni, BH, BY, fdr, none}. Default: BH
#' @param pvalueCutoff p-value threshold. Default: 0.05
#' @param qvalueCutoff q-value threshold. Default: 0.2
#' @param readable whether mapping gene ID to gene Name. Default: TRUE
#' @param keyType format of input key genes. Default: ENTREZID
#' @param verbose activate verbose mode. Default: False
#' @param split a string used to split genes stored into genes_sets. If is NULL, split process is avoided. Default: NULL
#' @param set_names names used instead set IDs (index)
#' @return a dataframe with all enrichments generated (without filtering)
#'  
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom clusterProfiler package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}
#' @seealso annot_GO
  
  # Load necessary packages
  require(clusterProfiler)
  if(verbose){
    require(pbapply)
  }
  
  # Check list of genes
  if(!is.null(split)){
    genes_sets <- lapply(genes_sets, function(set){unlist(strsplit(set,split))})
  }
  
  # Check if GO_DATA environment has been given
  if(is.null(GO_DATA)){ # Load GO set to avoid repetitive load
    GO_DATA <- clusterProfiler:::get_GO_data(OrgDb, ont, keyType)
  }# ELSE: Use given GO set
  
  
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
      enrich <- annot_GO(genes         = genes_sets[[i]],
                         GO_DATA       = GO_DATA,
                         OrgDb         = OrgDb,
                         ont           = ont,
                         pAdjustMethod = pAdjustMethod,
                         pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff,
                         readable      = readable,
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
      enrich <- annot_GO(genes         = genes_sets[[i]],
                         GO_DATA       = GO_DATA,
                         OrgDb         = OrgDb,
                         ont           = ont,
                         pAdjustMethod = pAdjustMethod,
                         pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff,
                         readable      = readable,
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
