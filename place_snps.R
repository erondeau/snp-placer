library(GenomicAlignments)
library(dplyr)

# rm(list=ls())

current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)


calculate_new_bp_data  <- function(sam_dataframe){
  
  sam_dataframe$adjusted_bp_SNP_location <- NA
  for (i in 1:length(rownames(sam_dataframe))){
      sam_dataframe$adjusted_bp_SNP_location[i] <- cigar_string_change(bp_of_snp = sam_dataframe$bp[i]
                                                                    , cigar_string = sam_dataframe$Cigar[i]
                                                                    , flag = sam_dataframe$Flag[i])
  }
  sam_dataframe$adjusted_bp_SNP_location <- as.numeric(sam_dataframe$adjusted_bp_SNP_location)
  
  sam_dataframe$alignment_length <- NA
  for (i in 1:length(rownames(sam_dataframe))){
    sam_dataframe$alignment_length[i] <- alignment_length(cigar_string = sam_dataframe$Cigar[i])
  }
  
  return(sam_dataframe)
}

snp_placement_dataframe <- function(sam_dataframe){
  # apply snp_contig_location across a dataframe
  sam_dataframe$contig_location <- NA
  for (i in 1:length(rownames(sam_dataframe))){
    sam_dataframe$contig_location[i] <- snp_contig_location(flag = sam_dataframe$Flag[i]
                                                            , pos = sam_dataframe$Pos[i]
                                                            , adjusted_bp_location = sam_dataframe$adjusted_bp_SNP_location[i]
                                                            , alignment_length = sam_dataframe$alignment_length[i])
  }
  return(sam_dataframe)
}


output_to_vcf <- function(output_df){
    #take the information in the dataframe and turn it into .vcf file format with the following header: 
    output_df$adj_name <- NA
    output_df$full_adj_name <- NA
    for (i in 1:length(rownames(output_df))){
      output_df$adj_name[i] <- paste0(output_df$Qname[i],"_",output_df$Polymorphism[i],"_",output_df$bp[i])
      output_df$full_adj_name[i] <- compliment_name(name=output_df$adj_name[i],flag=output_df$Flag[i])
    }
    
    vcf_out <- output_df[, c("Rname","contig_location","full_adj_name","Polymorphism","MapQ","Flag")]
    vcf_out$FILTER <- "."
    vcf_out$INFO <- "."
    vcf_out$QUAL <- "."
    vcf_out$REF <- NA
    vcf_out$ALT <- NA
    
    for (i in 1:length(rownames(output_df))){
        vcf_out$REF_CHECK[i] <- (unlist(strsplit(vcf_out$Polymorphism[i],split="/"),use.names = FALSE))[1]
        vcf_out$ALT_CHECK[i] <- paste0((unlist(strsplit(vcf_out$Polymorphism[i],split="/"),use.names = FALSE))[-1],collapse = ",")
    
        vcf_out$REF[i] <- allele_comp_check(in_allele = vcf_out$REF_CHECK[i], flag = vcf_out$Flag[i])
        vcf_out$ALT[i] <- allele_comp_check(in_allele = vcf_out$ALT_CHECK[i], flag = vcf_out$Flag[i])
        
    }
    vcf_out <- vcf_out[, c("Rname","contig_location","full_adj_name","REF","ALT","MapQ","FILTER","INFO")]
    colnames(vcf_out) <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
    return(vcf_out)
    
}


source("samParse.R")
source("cigarParse.R")

sam_header <- c('Qname','Flag','Rname','Pos','MapQ','Cigar','Rnext',
              'Pnext', 'TLEN', 'SEQ', 'QUAL','tag','type','value','value2')
sam_dat <- read.delim(file="one_location_alignments.sam",header = FALSE,stringsAsFactors = FALSE)

colnames(sam_dat) <- sam_header

snp_input_dat <- read.delim(file="SNPs.txt",header = TRUE,stringsAsFactors = FALSE)

sam_data_on_contigs <- filter(sam_dat,Qname %in% snp_input_dat$SNP)

all_polymorphism_data <- merge(sam_data_on_contigs,snp_input_dat,by.x="Qname",by.y="SNP")

all_polymorphism_data <- calculate_new_bp_data(sam_dataframe = all_polymorphism_data)

all_polymorphism_data <- snp_placement_dataframe(all_polymorphism_data)

polymorphism_vcf <- output_to_vcf(all_polymorphism_data)

write.table(polymorphism_vcf,file="placed_snps.vcf",sep='\t',quote = FALSE, row.names = FALSE, col.names = TRUE)


