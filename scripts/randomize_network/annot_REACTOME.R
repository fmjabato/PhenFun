annot_REACTOME <- function(genes, Reactome_DATA = NULL, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE){
#' Optimization of "ReactomePA" package REACTOME annotation method
#' avoiding computational expensive repetived tasks
#' @param genes a vector (or set of vectors) which includes a set of genes
#' @param Reactome_DATA
#' @param organism
#' @param pAdjustMethod
#' @param pvalueCutoff
#' @param qvalueCutoff
#' @param readable
#' @return dataframe with enrichment result
#'  
#' @author Fernando Moreno Jabato <fmjabato(at)gmail(dot)com>
#' @importFrom ReactomePA package
#' @seealso \link{https://bioconductor.org/packages/release/bioc/html/ReactomePA.html}
  
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
      body(enrichPathway)[[line_to_remove]] <- substitute("")      
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
