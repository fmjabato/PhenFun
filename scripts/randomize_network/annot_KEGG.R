annot_KEGG <- function(genes, KEGG_DATA = NULL, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "kegg"){
#' Optimization of "clusterProfiler" package KEGG annotation method
#' avoiding computational expensive repetived tasks
#' @param genes a vector (or set of vectors) which includes a set of genes
#' @param KEGG_DATA
#' @param organism
#' @param pAdjustMethod
#' @param pvalueCutoff
#' @param qvalueCutoff
#' @param keyType
#' @return dataframe with enrichment result
#'  
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom clusterProfiler package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html}
  
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
      body(enrichKEGG)[[line_to_remove]] <- substitute("")      
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
