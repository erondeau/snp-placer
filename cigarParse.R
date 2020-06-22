
#BiocManager::install(c("GenomicAlignments"))


cigar_cutter <- function(cigar){
    library(GenomicAlignments)
    cigar_match <- explodeCigarOps(cigar)
    cigar_length <- explodeCigarOpLengths(cigar)
    cigar_expand <- strrep(unlist(cigar_match,use.names = FALSE),unlist(cigar_length,use.names = FALSE))
    cigar_collapse <- paste(cigar_expand,collapse = "")
    cigar_tmp <- unlist(strsplit(cigar_collapse,split=""),use.names = FALSE)
    return(cigar_tmp)
}


fringe_snp_check <- function(bp_of_snp, cigar_dat){
  # look at the first and last tuples, if snp falls outside aligned region, return True"""
  sequence_len = 0
  cigar_log <- cigar_dat == "S" 
  if (min(which(cigar_log == FALSE)) > bp_of_snp){
    return(TRUE)
  } else if (max(which(cigar_log == FALSE)) < bp_of_snp){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


adjust_bp <- function(bp_of_snp, cigar_dat, flag){
  change_to_bp <- 0
  bp_scan <- 0
  if (flag == 16 || flag == 272){
      cigar_dat_2 <- rev(cigar_dat)
  } else {
      cigar_dat_2 <- cigar_dat
  }
  for (j in 1:length(cigar_dat_2)){
    
    if(cigar_dat_2[j] == 'M'){
      # count the matches towards the scan, no change to location
      bp_scan <- bp_scan + 1	
    } else if (cigar_dat_2[j] == 'D'){
      #add one to location, no change to scan count
      change_to_bp <- change_to_bp + 1
    } else if (cigar_dat_2[j] == 'I'){
      #minus one from location, move bp scan count up
      change_to_bp <- change_to_bp - 1
      bp_scan <- bp_scan + 1
      #exception: if the inserted bp housed the SNP
      if (bp_scan == bp_of_snp){
        return("snp_outside_aligned_region")
      }
    } else if (cigar_dat_2[j] == 'S' && (j != (length(cigar_dat_2) - 1))){
      # find the soft clipping strings, subtract from bp location"""
      change_to_bp <- change_to_bp - 1
      bp_scan <- bp_scan + 1 
    } else if (cigar_dat_2[j] == 'S'){
      bp_scan <- bp_scan + 1
    }
    if (bp_scan >= bp_of_snp){
      # if the scan has passed the snp, return the result as no more changes are needed
      return(as.numeric((bp_of_snp + change_to_bp)))
    }
  }
}




cigar_string_change <- function(bp_of_snp, cigar_string, flag){
    # take in the original string, and the snp location, adjust location based on
    # cigar data, returns a new bp integer that can be used relative to the start
    #	of the sequence's alignment to place the bp of the snp
    # NOTE: both the input and output string are NOT zero indexed """
    
    cigar_dat <- cigar_cutter(cigar_string)
    #first, identify the snps with no indels or font trimming
    if (all(cigar_dat == 'M') && length(cigar_dat) > bp_of_snp){
        return(bp_of_snp)
    } else {
        new_bp <- adjust_bp(bp_of_snp, cigar_dat, flag)
    }
    if (fringe_snp_check(bp_of_snp, cigar_dat) == TRUE){
        new_bp <- 'snp_outside_aligned_region'
    }
    return(new_bp)
}


alignment_length <- function(cigar_string){
    # take a list of cigar data tuples count total length of alignment: 
    #	'M' is match to the reference so these are counted
    #	'D' is deletion from  the reference, so these are counted in the length
    #	
    #	'I' is more sequence on the short read not on the referece, thefore
    #	don't add to length of sequence covered on the reference
    #	
    #	'S' is trimmed so don't add to length of sequence covered on the reference
    
    align_length <- 0
    for (pair in cigar_cutter(cigar_string)){
        if (pair == 'M' || pair == 'D')
        align_length = align_length + 1
    }
    return(align_length)
  
}



