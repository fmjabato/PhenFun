annot_GO <- function(genes, GO_DATA = NULL, OrgDb, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = T, keyType = "ENTREZID"){
#' Optimization of "clusterProfiler" package GO annotation method
#' avoiding computational expensive repetived tasks
#' @param genes a vector (or set of vectors) which includes a set of genes
#' @param GO_DATA
#' @param OrgDB
#' @param ont
#' @param pAdjustMethod
#' @param pvalueCutoff
#' @param qvalueCutoff
#' @param readable
#' @param keyType
#' @return dataframe with enrichment result
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
      body(enrichGO)[[line_to_remove]] <- substitute("")      
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
