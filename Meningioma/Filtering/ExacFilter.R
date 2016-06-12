## Noah Greenwald and PK Agarwalla
## Checks MAF file against EXAC database of known germline mutations, with threshold AF in order to be called as germline
## Assuming maximum number of alleles (120,000), a variant seen in 12 samples will have an AF of 12/120,000 = .0001

run.exac <- function(maf, af = 0){
  
  total.exac <- read.delim('C:/Users/Noah/OneDrive/Work/Coding/R/dbs/germline_37k_split_mult_alleles.AF.INFO',
                                      stringsAsFactors = FALSE)
    
  idx <- which(total.exac$POS %in% maf$Start_position)
  exac <- total.exac[idx, ]
  
  # Goal is to annotate the maf with the germline data
  # Loop through entire maf
  for (i in 1:nrow(maf)){
    
   ## Prints out every 1,000 cycles
    if (i %% 1000 == TRUE){
      cat('Checking maf row', i, '\n')
    }
    
    ## Get all matches at current position
    match.idx <- which((exac$CHROM == maf$Chromosome[i]) & (exac$POS == maf$Start_position[i]) & exac$AF > af)
  
  ## Deal with SNPs first   
  if (maf$Variant_Type[i] == "SNP"){
    
    ## Keeps only those exac hits which have single nucleotide changes
    temp.idx <- which(nchar(exac$REF[match.idx]) & nchar(exac$ALT[match.idx]) == 1)
    snp.idx <- match.idx[temp.idx]
    
    # Make sure SNP at given position is same base pair change
    tmp.idx <- which(maf$Tumor_Seq_Allele2[i] == exac$ALT[snp.idx])
    snp.idx <- snp.idx[tmp.idx]
    
    if (length(snp.idx) > 0){
      maf$POS[i] <- exac$POS[snp.idx]
      maf$REF[i] <- exac$REF[snp.idx]
      maf$ALT[i] <- exac$ALT[snp.idx]
      maf$AF[i] <- exac$AF[snp.idx]
      maf[i,'germline']=TRUE
    }else{
      maf$POS[i] <- NA
      maf$REF[i] <- NA
      maf$ALT[i] <- NA
      maf$AF[i] <- 0
      maf[i,'germline']=FALSE
    }
    
    ## Now look at indels
  } else {
    
    ## Keep only those exac hits which are indels
    temp.idx <- which(nchar(exac$REF[match.idx]) | nchar(exac$ALT[match.idx]) != 1)
    indel.idx <- match.idx[temp.idx]
    
    ## Flag as indel, don't check specifics
    if (length(indel.idx > 1)){
    indel.idx <- indel.idx[1]
    }
    ## TODO: indel criteria??
    
    if (length(indel.idx) > 0){
      maf$POS[i] <- exac$POS[indel.idx]
      maf$REF[i] <- exac$REF[indel.idx]
      maf$ALT[i] <- exac$ALT[indel.idx]
      maf$AF[i] <- exac$AF[indel.idx]
      maf[i,'germline']=TRUE
    }else{
      maf$POS[i] <- NA
      maf$REF[i] <- NA
      maf$ALT[i] <- NA
      maf$AF[i] <- 0
      maf[i,'germline']=FALSE
    }
    
  }
  }

return(maf)

}