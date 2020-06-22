snp_contig_location <- function(flag, pos, adjusted_bp_location, alignment_length){
  
    if (is.numeric(adjusted_bp_location)){

    } else {
        return("-")
    }

    if (flag == 0 || flag == 256){
      return((pos + adjusted_bp_location - 1))
    } else if (flag == 16 || flag == 272){
      return((pos + alignment_length - adjusted_bp_location))
    } else {
      return("-")
    }
}

compliment_name <- function(name, flag){
      #if the alignment is a reverse, add _comp to the end of its identification
      if (flag == 16 || flag == 272){
          return(paste0(name,"_comp"))
      } else {
          return(name)
      }
}


match_snp <- function(allele){
  allele <- toupper(allele)
  if (allele == "A"){
    return("T")
  } else if (allele == "T"){
    return("A")
  } else if (allele == "C"){
    return("G")
  } else if (allele == "G"){
    return("C")
  } else if (allele == "N"){
    return("N")
  } else { 
    return(paste0("Need valid nucleotide (ATGC) or N\n, not valid = ",allele))
    
  }
  
  
}


allele_comp_check <- function(in_allele, flag){
  # if alignment is a reverse, flip the alleles to the complimentary nucleotides """
  if (flag == 0 || flag == 256){
      return(in_allele)
  } else if (flag == 16 || flag == 272){
    in_allele <- unlist(strsplit(in_allele,",",fixed=TRUE),use.names=FALSE)
    if (length(in_allele) == 1){
       return(match_snp(in_allele))
    } else {
       out_alleles <- NULL
       for (i in 1:length(in_allele)){
          out_alleles[i] <- (match_snp(in_allele[i]))
       }
       out_alleles <- paste(out_alleles,collapse=",")
       return(out_alleles)
    }
  }
}
