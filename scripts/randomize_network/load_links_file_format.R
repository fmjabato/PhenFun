

load_links_file_format <- function(file,sep="\t",header=F){
#' @author Fernando Moreno Jabato <jabato@uma.com>
#' @RVersion 3.4.2
#' @description Method to load a table of relationships loaded in a specific format into a given file.
#' After load, info is parsed and returned in the correct format.
#' Allowed formats are: HPO_HPO_Format and HPO_Region_Format which have the following
#' formats:
#'   > HPO_HPO_Format:
#'     - Header: NO
#'     - Sep: "\t"
#'     - Line info format: <HPO> <HPO> <Intensity>
#'     - Line info types: {c,c,d}
#'     - Note: third column (intesity) is optional
#'   > HPO_Region_Format:
#'     - Header: NO
#'     - Sep: "\t"
#'     - Line info format: <Chr> <Start> <End> <HPO> <Intensity> <NodeID>
#'     - Line info types: {c,i,i,c,d,c}
#'   > Loci_HPO_format:
#'     - Header: NO
#'     - Sep: "\t"
#'     - Line info format: <NodeID> <HPO>
#'     - Line info types: {c,c}
#' @param file with info to be loaded and parsed
#' @param sep columns separator used into file
#' @param header indicates if file has, or not, header
#' @return a dataframe with formated info or a string with an error description if any error occurs
  
  # Config values (A = HPO_HPO ; B = HPO_Region)
  ncol_formatA <- 3
  ncol_formatB <- 6
  ncol_specialFormat <- 2
  
  # Check inputs
  if(is.null(file)){ # File
    return("ERROR [load_links_file_format]: Given file is NULL pointer")
  }else if(!is.character(file)){
    return("ERROR [load_links_file_format]: Given file is not a character value")
  }else if(!file.exists(file)){
    return("ERROR [load_links_file_format]: Given file does not exists")
  }else if(file.access(file,mode=4) != 0){ # Mode 4 = Read permission
    return("ERROR [load_links_file_format]: Given file can not be read")
  }else if(file.size(file) == 0){
    return("ERROR [load_links_file_format]: Given file is empty file")
  }
  
  if(is.null(sep)){
    return("ERROR [load_links_file_format]: Given separator is NULL pointer")
  }else if(!is.character(sep)){
    return("ERROR [load_links_file_format]: Given seperator is not a character value")
  }else if(nchar(sep) == 0){
    return("ERROR [load_links_file_format]: Given separator is an empty string")
  }
  
  if(is.null(header)){
    return("ERROR [load_links_file_format]: Given Header flag is NULL pointer")
  }else if(!is.logical(header)){
    return("ERROR [load_links_file_format]: Given Header flag is not a logical value")
  }
  
  # Load file info
  info <- read.table(file,sep = sep, header = header, stringsAsFactors = F)
  
  # Check 
  if(length(dim(info)) != 2){
    return("ERROR [load_links_file_format]: Read info is not a table")
  }else if(any(dim(info)==0)){
    return("ERROR [load_links_file_format]: Table is empty or is vector")
  }
  
  # Check format
  if(ncol(info) == ncol_formatA | ncol(info) == ncol_specialFormat){   # FORMAT: HPO-HPO (weighted or not) or Loci-HPO
    # Check if it's a Loci-HPO file
    if(ncol(info) == ncol_specialFormat){
      if(all(grepl("^HP:[0-9]{7}$", info[,2])) & !all(grepl("^HP:[0-9]{7}$",info[,1]))){
        info <- as.data.frame(info,stringsAsFactors = F)
        colnames(info) <- c("Loci","HPO")
        info[,1] <- as.character(info[,1])
        info[,2] <- as.character(info[,2])
        return(info)
      }else if(all(grepl("^HP:[0-9]{7}$", info[,1]))){
        info <- as.data.frame(info,stringsAsFactors = F)
#        colnames(info) <- c("HPO","Gene")
        colnames(info) <- c("Loci","HPO")
        info[,1] <- as.character(info[,1])
        info[,2] <- as.character(info[,2])
        return(info)
      }else return("ERROR: Test")
    }
    # Check each column
    if(!all(grepl("^HP:[0-9]{7}$",info[,1]))){       # Column 1: HP code
      return("ERROR [load_links_file_format]: Column 1 has not correct HPO-HPO format")
    }else if(!all(grepl("^HP:[0-9]{7}$",info[,2]))){ # Column 2: HP code
      return("ERROR [load_links_file_format]: Column 2 has not correct HPO-HPO format")
    }else if(ncol(info) == 3){
      weighted <- T
      if(any(is.na(as.numeric(info[,3])) | !is.double(as.numeric(info[,3])))){ # Column 3: Relationship weight
        return("ERROR [load_links_file_format]: Column 3 has not correct HPO-HPO format")
      }
    }else{ # Everything OK
      weighted <- F
    }
    # Parse format
    info <- as.data.frame(info)
    info[,1] <- as.character(info[,1])
    info[,2] <- as.character(info[,2])
    if(weighted){
      info[,3] <- as.double(info[,3])
      colnames(info) <- c("HPO_1","HPO_2","Value") 
    }else{
      colnames(info) <- c("HPO_1","HPO_2")
    }
    # return
    return(info)
  }else if(ncol(info) == ncol_formatB){ # FORMAT: HPO-REGION
    # Check each column
    if(!all(grepl("(^[1-9]$)|(^1[0-9]$)|(^2[0-2]$)|(^[X,Y]$)",info[,1]))){ # Column 1: Chromosome
      return("ERROR [load_links_file_format]: Column 1 has not correct HPO-Region format")
    }else if(!all(grepl("^[0-9]*$",info[,2]))){              # Column 2: Start coord
      return("ERROR [load_links_file_format]: Column 2 has not correct HPO-Region format")
    }else if(!all(grepl("^[0-9]*$",info[,3]))){              # Column 3: End coord
      return("ERROR [load_links_file_format]: Column 3 has not correct HPO-Region format")
    }else if(!all(grepl("^HP:[0-9]{7}$",info[,4]))){         # Column 4: HP code
      return("ERROR [load_links_file_format]: Column 4 has not correct HPO-Region format")
    }else if(any(is.na(as.numeric(info[,5])) | !is.double(as.numeric(info[,5])))){ # Column 5: Relationship weight
      return("ERROR [load_links_file_format]: Column 5 has not correct HPO-Region format")
    }else if(!all(grepl("^(([1-9])|(1[0-9])|(2[0-2])|([X,Y]))\\.[0-9]{1,3}\\.[A-Z]\\.[0-9]{1,3}$",info[,6]))){
      return("ERROR [load_links_file_format]: Column 6 has not correct HPO-Region format")
    }else{ # Everything OK
      # Parse info
      info <- as.data.frame(info)
      info[,1] <- as.character(info[,1])
      info[,2] <- as.integer(info[,2])
      info[,3] <- as.integer(info[,3])
      info[,4] <- as.character(info[,4])
      info[,5] <- as.double(info[,5])
      info[,6] <- as.character(info[,6])
      colnames(info) <- c("Chr","Start","End","HPO","Value","Node")
      # Return
      return(info)
    }
  }else{                                # FORMAT: NOT ALLOWED
    return("ERROR [load_links_file_format]: Info format is not allowed")
  }
}
