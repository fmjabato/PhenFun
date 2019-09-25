#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse, topGO

library(optparse)
library(topGO)
######################################################################################
## FUNCTIONS
######################################################################################
get_collapsed_data <- function(data){ 

        # Collapse by Key
        genes_sep <- ":"
        tnetwork <- as.data.frame(do.call(rbind,lapply(unique(data$ID),function(id){
          # Find tuples
          indx <- which(data$ID == id)
          # Collapse tuples
          genes <- paste(data$Gene[indx], collapse = genes_sep)
          # Return info
          return(list(ID    = id,
                      Genes = genes))
        })))
        # Unlist columns
        for(i in seq(ncol(tnetwork))){
          tnetwork[,i] <- unlist(tnetwork[,i])
          rm(i)
        }
	tnetwork$Genes <- lapply(tnetwork$Genes, function(set){unlist(strsplit(set, genes_sep))})

	#tnetwork$Genes <- lapply(tnetwork$Genes, function(set){strsplit(set, genes_sep)})
	return(tnetwork)
}

########################################################################
## OPTPARSE
########################################################################

option_list <- list(
	make_option(c("-d", "--data"), type="character",
		help="Path to input file"),
	make_option(c("-o","--output"), type="character",
		help="Output report"),
	make_option(c("-s","--sub_ontology"), type = "character",
		help="GO sub-ontology to perform analysis. Can be MF, CC or BP"),
	make_option(c("-p","--pval"), type = "double", default = 0.05,
		help="GO sub-ontology to perform analysis. Can be MF, CC or BP"),
	make_option(c("-t","--tag"), type = "character",
		help="Tag to identify whole analysis")
)

opt <- parse_args(OptionParser(option_list=option_list))


######################################################################################
## MAIN
#####################################################################################
data <- read.table(opt$data, sep = "\t")
names(data) <- c("ID","Gene")
tag2gene <- get_collapsed_data(data)
#https://rdrr.io/bioc/topGO/man/annFUN.html
xx <- annFUN.org(opt$sub_ontology, mapping = "org.Hs.eg.db")
allGenes <- unique(unlist(xx))
all_results_table <- data.frame()
for(i in 1:nrow(tag2gene)){ #analyse each  gene vector
	record <- tag2gene[i, ]
	id <- record$ID
	geneNames <- unlist(record$Genes)
	# creating contigency vector for gene identifiers
	geneList <- factor(as.integer(allGenes %in% geneNames))
	names(geneList) <- allGenes
	if(length(levels(geneList)) == 2){ # launch analysis
		TopGOobject <- new("topGOdata", ontology = opt$sub_ontology, allGenes = geneList, annot = annFUN.org, mapping = "org.Hs.eg.db")
		classic <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher_Test")
		resultsclassic <- getSigGroups(TopGOobject, classic)
		tops <- length(score(resultsclassic)[score(resultsclassic) < opt$pval]) # score is a slot of resultsclassic
		if(tops > 0){
			results_table <- GenTable(TopGOobject, classic = resultsclassic, orderBy = "classic", ranksOf = "classic", topNodes = tops)
			# adding go data to gene table
			genes_in_goes <- genesInTerm(TopGOobject, results_table$GO.ID)
			colnames(results_table)[6] <- "Fisher Test" # change the column name 'classic'
			intersected <- lapply(genes_in_goes, function(x) intersect(x, geneNames))
			genes_of_goes <- unlist(lapply(intersected, function(x) paste(x, collapse = ":")))
			results_table <- cbind(results_table, genes_of_goes)
			colnames(results_table)[7] <- "Significant Genes"
			results_table$ID <- rep(id, tops)
			results_table$sub_ontology <- rep(opt$sub_ontology, tops)
			if(!is.null(opt$tag)){
				results_table$tag <- rep(opt$tag, tops)
			}
			if(nrow(all_results_table) == 0){
				all_results_table <- results_table
			}else{
				all_results_table <- rbind(all_results_table, results_table) 
			}
		}
	}
}
write.table(all_results_table, file = opt$output, sep = "\t", quote = FALSE, row.names = F, col.names = T)
