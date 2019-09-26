#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse, ontologyIndex

################################################################
##                           CONFIGURE                        ##
################################################################
require(optparse)
require(ontologyIndex)

option_list <- list(
  make_option(c("-s", "--source"), action="store", type="character",
              dest="source", help="Source network (not enriched)"),
  make_option(c("-n", "--net"), action="store", type="character",
              dest="net", help="Used network"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file")
)

opt <- parse_args(OptionParser(option_list=option_list))


################################################################
##                           LOAD DATA                        ##
################################################################
# Load source network
source <- read.table(file = opt$source, sep = "\t", quote = "", header = FALSE, stringsAsFactors = FALSE)
source_indx <- which(grepl("HP",source[,1]))

# Load used network
net <- read.table(file = opt$net, sep = "\t", quote = "", header = FALSE, stringsAsFactors = FALSE)
net_indx <- which(grepl("HP",net[,1]))



################################################################
##                         PROCESS DATA                       ##
################################################################
# Prepare metadata container
meta_subject <- c()
meta_value   <- c()

# Obtain SOURCE info
pheno_source <- unique(source[source_indx,1])
sor_source   <- unique(source[-source_indx,2])
pats_source  <- unique(source[source_indx,2])

# Obtain NET info
pheno_net <- unique(net[net_indx,1])
sor_net   <- unique(net[-net_indx,2])
pats_net  <- unique(net[net_indx,2])

# Make base comparissons
pheno_shared <- intersect(pheno_source,pheno_net)
sor_shared   <- intersect(sor_source,sor_net)
pats_shared  <- intersect(pats_source,pats_net)

# Add base metadata
meta_subject <- c(meta_subject,"Source net Phenotypes","Used net Phenotypes","Shared Phenotypes","Source net Patients","Used net Patients","Shared Patients","Source net SORs","Used net SORs","Shared SORs")
meta_value   <- c(meta_value,length(pheno_source),length(pheno_net),length(pheno_shared),length(pats_source),length(pats_net),length(pats_shared),length(sor_source),length(sor_net),length(sor_shared))




################################################################
##                      REMOVE PARENTALS                      ##
################################################################
# Define useful function
remove_parentals <- function(terms,ontology){
  # Per each term: check parentals
  to_remove <- lapply(terms,function(term){
    parentals <- head(get_ancestors(ontology,term),-1)
    # Check if are into terms
    indx <- which(terms %in% parentals)
    # Return
    return(indx)
  })
  # Obtain terms to remove
  to_remove <- unique(unlist(to_remove))
  # Remove and return
  if(length(to_remove > 0)){
    filtered_terms <- terms[-to_remove]
    return(filtered_terms)    
  }

  return(terms)
}

# Load ontology
data(hpo)

# Remove parentals
pheno_source_filtered <- remove_parentals(pheno_source,hpo)
pheno_net_filtered    <- remove_parentals(pheno_net,hpo)

# Calculate intersects
pheno_shared_filtered     <- intersect(pheno_source_filtered,pheno_net_filtered)
pheno_not_shared_filtered <- setdiff(pheno_source_filtered,pheno_net_filtered)


# Add meta_info
meta_subject <- c(meta_subject,"Source net Phenotypes after remove parentals", "Used net Phenotypes after remove parentals")
meta_value   <- c(meta_value,length(pheno_source_filtered), length(pheno_net_filtered))

# # Check used child terms not included into source net 
# pheno_recovered <- setdiff(pheno_net_filtered,pheno_source)
# pheno_lost      <- setdiff(pheno_source_filtered,pheno_net_filtered)
# meta_subject <- c(meta_subject,"Recovered childs into used not included into source", "Lost childs into used included into source")
# meta_value   <- c(meta_value,length(pheno_recovered),length(pheno_lost))

pheno_exclusive <- unique(unlist(lapply(pheno_not_shared_filtered,function(term){
  childs <- get_descendants(hpo,term)
  if(any(childs %in% pheno_net)){
    return(term)
  }
  return(c())
})))

# Update info
meta_subject <- c(meta_subject,"Used net Phenotypes filtered present into used net filtered","Used net Phenotypes filtered with childs present into used net filtered","Used net exclusive Phenotypes filtered")
meta_value <- c(meta_value,length(pheno_shared_filtered),length(pheno_net_filtered)-length(pheno_shared_filtered)-length(pheno_exclusive),length(pheno_exclusive))




################################################################
##                            WRITE                           ##
################################################################
# Create output dataframe
out <- data.frame(Concept = meta_subject, Value = meta_value, stringsAsFactors = FALSE)
# Export
write.table(out, file = opt$output,quote = FALSE, row.names = FALSE)