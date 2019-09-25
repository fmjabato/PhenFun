#! /usr/bin/env Rscript

#' @description Perform different randomizations over a network
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse

##############################################################################
##                           CONFIGURE PROGRAM                              ##
##############################################################################

# Load necessary packages
suppressPackageStartupMessages(require(optparse))        # Parse script inputs

# Prepare input commands
option_list <- list(
  make_option(c("-n", "--network"), action="store", type="character",
              dest="flinks", help="Network file"),
  make_option(c("-r","--RDMNumRels"), action="store", type = "character", default = "links",
              dest="rdm_num_rels", help="Randomization type. Values allowed: 'links', 'nodes', 'all'. [Default = %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              dest="verbose",help="Activate verbose mode"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Destiny file")
)

opt <- parse_args(OptionParser(option_list=option_list))

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                       error=function(e) # works when using R CMD
                         normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
bname <- dirname(full.fpath)

# Import functions
if(is.null(bname)){
  source_files <- list.files(path = paste("randomize_network",.Platform$file.sep,sep=""), pattern = "\\.R$", full.names = T)
}else{
  source_files <- list.files(path = paste(bname,"randomize_network",sep=.Platform$file.sep), pattern = "\\.R$",full.names = T)
}
invisible(lapply(source_files, source))
to_remove <- c("to_remove","bname","source_files","script.name","script.path","initial.options")
rm(list = ls()[which(ls() %in% to_remove)])

# Configure pb package
if(opt$verbose){
  pboptions(type="timer")
}

##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################

# Verbose point
if(opt$verbose){
  message("Loading original network")
}

# Load source network
net <- load_links_file_format(opt$flinks)

if(is.character(net)){
  stop(net)
}

##############################################################################
##                             RANDOMIZE                                    ##
##############################################################################
if(opt$verbose){
  message("Randomizing network")
}

# Generate random network
if(opt$rdm_num_rels == "links"){
  rdm_net <- randomize_network(network = net, rdm_num_rels = FALSE)
}else if(opt$rdm_num_rels == "all"){
  rdm_net <- randomize_network(network = net, rdm_num_rels = TRUE)
}else if(opt$rdm_num_rels == "nodes"){
  rdm_net <- randomize_nodes(net)
}else{
  stop("Main error: randomization types given is not allowed")
}

if(opt$verbose){
  message(paste("Writing netowrk into:",opt$output))
}
# Write networks
write.table(rdm_net, append = F, sep = "\t",row.names = F,col.names = F,quote = F, file = opt$output)

 
