#' @author Fernando Moreno Jabato <jabato@uma.com>
#' @description function to randomize links into a network with an specific format
#' @RVersion 3.4.2
#' @method randomize_network
#' @method randomize_nodes



randomize_network <- function(network, rdm_num_rels = FALSE){
#' Method to randomize a weighted network in a specific format.
#' After randomization process, a network with the same elements,
#' number of links and format but with different configuration is
#' returned.
#' If a weighted network is formed by {X,Y} nodes and {W} weight 
#' values, number of links of each Xi and Yi will be changed and 
#' weights will be assigned randomly
#' Specific formats allowed are:
#'  > HPO-HPO network
#'   - Num_columns = 2 or 3 (weighted or not)
#'   - Columns = <HPO_1>; <HPO_2>; <Value>
#'  > HPO-Region network
#'   - Num_columns = 6
#'   - Columns = <Chr>; <Start>; <End>; <HPO>; <Value>; <Node>
#' @param network to be randomized. Must be a dataframe
#' @param rdm_num_rels flag to indicate if number of relationships of each node must be randomized too
#' @return a network with same properties but different configuration (randomized)
#'   or a character with an error message if any error occurs
#' @require pbapply package if verbose mode is activated
  
  
  # Config values (A = HPO_HPO ; B = HPO_Region)
  ncol_formatA <- 3
  ncol_formatB <- 6
  ncol_specialFormat <- 2
  
  needed_formatA <- c("HPO_1","HPO_2","Value")
  needed_formatB <- c("Chr","Start","End","HPO","Value","Node")
  needed_formatS <- c("HPO_1","HPO_2")
  
  # Check inputs
  if(is.null(network)){
    return("ERROR [randomize_network]: Given network is NULL pointer") 
  }else if(!is.data.frame(network)){
    return("ERROR [randomize_network]: Given network is not a data frame")
  }else if(ncol(network) != ncol_formatA & ncol(network) != ncol_formatB & ncol(network) != ncol_specialFormat){
    return("ERROR [randomize_network]: Given network has not allowed dimensions")
  }else if(!all(colnames(network) %in% needed_formatA) & !all(colnames(network) %in% needed_formatB) & !all(colnames(network) %in% needed_formatS)){
    return("ERROR [randomize_network]: Given network has not allowed columns format")
  }

  # Identify 
  
  # Identify format
  if(ncol(network) == ncol_formatA | ncol(network) == ncol_specialFormat){ # HPO-HPO
    if(ncol(network) == ncol_specialFormat){
      weighted <- F
    }else{
      weighted <- T
    }
    withoutErrs <- FALSE
    while(!withoutErrs){
      # Obtain number of relationships per term
      hpo_rels <- sort(table(c(network$HPO_1,network$HPO_2)),decreasing = T)
      if(weighted){
        # Obtain vector of relationships weights
        hpo_weig <- network$Value
      }
      
      # Randomize number of relationships of each HPO term
      if(rdm_num_rels){
        names(hpo_rels) <- sample(names(hpo_rels))
        if(weighted){
          # Randomize weights
          hpo_weig <- sample(hpo_weig) 
        }
      }

      # Prepare hpos
#      hpos <- names(hpo_rels)
#      names(hpo_rels) <- seq_along(hpos)
      
      # Generate new data frame
      options(stringsAsFactors = F)
      if(weighted){
        rdm_net <- data.frame(HPO_1 = character(nrow(network)),
                              HPO_2 = character(nrow(network)),
                              Value = numeric(nrow(network)),
                              stringsAsFactors = F) 
      }else{
        rdm_net <- data.frame(HPO_1 = character(nrow(network)),
                              HPO_2 = character(nrow(network)),
                              stringsAsFactors = F)
      }
      already_matched <- as.list(rep("",length(hpo_rels)))
      
      
      # Randomize
      result = tryCatch({
        trash <- lapply(seq(nrow(network)),function(i){
          # Check
          if(length(hpo_rels) <= 1){
            stop()
          }else if(length(hpo_rels) == 2){
            rdm_nodes <- c(1,2)
          }else{
            # Take random samples
            rdm_nodes <- c(1,sample(seq_along(hpo_rels)[-1],1))  
          }
          
          # Update already matched nodes
          already_matched[[rdm_nodes[1]]] <<- c(already_matched[[rdm_nodes[1]]],names(hpo_rels)[rdm_nodes[2]]) 
          already_matched[[rdm_nodes[2]]] <<- c(already_matched[[rdm_nodes[2]]],names(hpo_rels)[rdm_nodes[1]])
          # Generate new edge
          entry <- list(HPO_1 = names(hpo_rels)[rdm_nodes[1]],
                        HPO_2 = names(hpo_rels)[rdm_nodes[2]])
          if(weighted){
            entry <- c(entry,Value = hpo_weig[i])
          }
          
          # Store edge
          rdm_net$HPO_1[i] <<- entry$HPO_1
          rdm_net$HPO_2[i] <<- entry$HPO_2
          if(weighted){
            rdm_net$Value[i] <<- entry$Value
          }

          # Update current node frequencies
          hpo_rels[rdm_nodes[1]] <<- hpo_rels[rdm_nodes[1]] - 1
          hpo_rels[rdm_nodes[2]] <<- hpo_rels[rdm_nodes[2]] - 1
          # Check if any must be removed
          to_remove <- c()
          if(hpo_rels[[rdm_nodes[1]]] <= 0){
            to_remove <- rdm_nodes[1]
          }
          if(hpo_rels[[rdm_nodes[2]]] <= 0){
            to_remove <- c(to_remove,rdm_nodes[2])
          }
          # Remove if it's necessary
          if(length(to_remove) > 0){
            hpo_rels <<- hpo_rels[-to_remove]
            already_matched <<- already_matched[-to_remove]
          }
        })
        # Everything finished OK
        withoutErrs <- TRUE
      }, 
      warning = function(w,verbose = F){if(verbose){warning("A warning has been launched. Retrying model")}},
      error = function(e,verbose = F){if(verbose){warning(paste("An error has been launched(",e,"). Retrying model"))}})
    }
    
  }else{ # HPO-Region
    ############################################################### WARNING!!
    message("Error, this functionality needs revision. Implementation is not correct. Returning NULL")
    return(NULL)
    ###############################################################
    
    
    
    # Obtain number of relationships per HPO term and per Node
    hpo_rels  <- table(network$HPO)
    node_rels <- table(network$Node) 
    # Obtain links weights
    rel_weights <- network$Value
    # Prepare data frame with NODEs info
    nodes_info <- as.data.frame(do.call(rbind,lapply(names(node_rels), function(node){
      # Find node into network
      indx <- which(network$Node == node)[1]
      # Take extra data and return
      return(list(Chr   = network$Chr[indx],
                  Start = network$Start[indx],
                  End   = network$End[indx],
                  Node  = node))
    })))
    # Randomize number of realtionships
    if(rdm_num_rels){
      # Randomize HPO number of relationships
      names(hpo_rels) <- sample(names(hpo_rels))
      # Randomize Node number of relationships
      names(node_rels) <- sample(names(node_rels)) 
    }
    # Randomize weights
    rel_weights <- sample(rel_weights)
    # Randomize relatioships    
    rdm_net <- as.data.frame(do.call(rbind,lapply(sample(rel_weights),function(i){
      # Obtain random samples
      rdm_hpo <- sample(seq_along(hpo_rels),1)
      rdm_nod <- sample(seq_along(node_rels),1)
      # Find necessary info
      node_indx <- which(nodes_info$Node == names(node_rels)[rdm_nod])
      # Generate new entry
      entry <- list(Chr   = nodes_info$Chr[node_indx],
                    Start = nodes_info$Start[node_indx],
                    End   = nodes_info$End[node_indx],
                    HPO   = names(hpo_rels)[rdm_hpo],
                    Value = i,
                    Node  = nodes_info$Node[node_indx])
      # Update relationships
      hpo_rels[rdm_hpo]  <<- hpo_rels[rdm_hpo] - 1
      node_rels[rdm_nod] <<- node_rels[rdm_nod] - 1
      # Remove if it's necessary
      if(node_rels[rdm_nod] == 0){
        node_rels <<- node_rels[-rdm_nod]
      }
      if(hpo_rels[rdm_hpo] == 0){
        
        hpo_rels <<- hpo_rels[-rdm_hpo]
      }
      # Return new entry
      return(entry)
    })))
  }
  

  # Return random network
  return(rdm_net)
}



randomize_nodes <- function(network){

  # Configure allowed formats
  ncol_formatA          <- 2
  ncol_formatA_weighted <- 3
  ncol_formatB          <- 6
  
  needed_formatW <- c("HPO_1","HPO_2","Value")
  needed_formatB <- c("Chr","Start","End","HPO","Value","Node")
  needed_formatA <- c("HPO_1","HPO_2")
  needed_formatC <- c("Loci","HPO")

  # Check inputs
  if(is.null(network)){
    return("ERROR [randomize_network]: Given network is NULL pointer") 
  }else if(!is.data.frame(network)){
    return("ERROR [randomize_network]: Given network is not a data frame")
  }else if(ncol(network) != ncol_formatA & ncol(network) != ncol_formatB & ncol(network) != ncol_formatA_weighted){
    return("ERROR [randomize_network]: Given network has not allowed dimensions")
  }else if(!all(colnames(network) %in% needed_formatA) & !all(colnames(network) %in% needed_formatB) & !all(colnames(network) %in% needed_formatW) & !all(colnames(network) %in% needed_formatC)){
    return("ERROR [randomize_network]: Given network has not allowed columns format")
  }

  # Identify format
  if(ncol(network) == ncol_formatA & all(colnames(network) %in% needed_formatC)){
    # Obtain nodes
    nodes <- unique(network$Loci)
    # Substitute network nodes using indexes IDs
    rdm_net <- network
    invisible(lapply(seq_along(nodes),function(i){
      rdm_net$Loci[which(rdm_net$Loci == nodes[i])] <<- i
    }))
    # Randomize nodes
    nodes <- nodes[sample(seq(length(nodes)))]
    # Substitute nodes
    invisible(lapply(seq_along(nodes), function(i){
      rdm_net$Loci[which(rdm_net$Loci == i)] <<- nodes[i]
    }))
  }else if(ncol(network) == ncol_formatA | ncol(network) == ncol_formatA_weighted){
    # Obtain nodes
    nodes <- unique(c(network$HPO_1,network$HPO_2))
    # Substitute network nodes using indexing IDs
    rdm_net <- network
    invisible(lapply(seq_along(nodes),function(i){
      rdm_net$HPO_1[which(rdm_net$HPO_1 == nodes[i])] <<- i
      rdm_net$HPO_2[which(rdm_net$HPO_2 == nodes[i])] <<- i
    }))
    # Randomize nodes
    nodes <- nodes[sample(seq(length(nodes),length(nodes)))]
    # Substitute new nodes
    invisible(lapply(seq_along(nodes),function(i){
      rdm_net$HPO_1[which(rdm_net$HPO_1 == i)] <<- nodes[i]
      rdm_net$HPO_2[which(rdm_net$HPO_2 == i)] <<- nodes[i]
    }))
  }else if(ncol(network) == ncol_formatB){
    # Obtain nodes
    nodes <- unique(network(HPO))
    # Substitute network nodes using indexing IDs
    rdm_net <- network
    invisible(lapply(seq_along(nodes),function(i){
      rdm_net$HPO[which(rdm_net$HPO == nodes[i])] <<- i
    }))
    # Randomize nodes
    nodes <- nodes[sample(seq(length(nodes),length(nodes)))]
    # Substitute new nodes
    invisible(lapply(seq_along(nodes),function(i){
      rdm_net$HPO[which(rdm_net$HPO == i)] <<- nodes[i]
    }))
  }

  # Return generated network
  return(rdm_net)
}
