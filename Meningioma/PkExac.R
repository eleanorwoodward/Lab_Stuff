pk.exac <- function(maf){
  exac2 <- read.delim('C:/Users/Noah/OneDrive/Work/R/germline_37k_split_mult_alleles.AF.INFO',
                      stringsAsFactors = FALSE)
  starts <- which(exac2$POS %in% maf$Start_position)
  ends <- which(exac2$POS %in% maf$End_position)
  
  idx <- unique(c(starts, ends))
  
  exac <- exac2[idx, ]
  # Goal is to annotate the maf with the germline data
  
  for (i in 1:nrow(maf)){
    cat('Checking maf row', i, '\n')
    match.idx <- which(exac$CHROM == maf$Chromosome[i] & (exac$POS ==
                                                            maf$Start_position[i] | exac$POS == maf$End_position[i]))
    if(length(match.idx)>1){
      highest_AF_idx <- which.max(exac[match.idx,'AF'])
      match.idx <- match.idx[highest_AF_idx]
    }
    if (length(match.idx) > 0){
      maf$POS[i] <- exac$POS[match.idx]
      maf$REF[i] <- exac$REF[match.idx]
      maf$ALT[i] <- exac$ALT[match.idx]
      maf$AF[i] <- exac$AF[match.idx]
      maf[i,'germline']=TRUE
    }else{
      maf$POS[i] <- NA
      maf$REF[i] <- NA
      maf$ALT[i] <- NA
      maf$AF[i] <- NA
      maf[i,'germline']=FALSE
    }
    
    #maf$germline[i] <- as.logical(is.na(maf$POS[i]))
    
    #if (maf[i, 'i_COSMIC_n_overlapping_mutations'] > 0 | (maf[i,
    #'germline'] == FALSE & (maf[i, 'dbSNP_RS'] == "" |
                              #is.na(maf[i,'dbSNP_RS'])))){
                                # maf[i, 'KEEP'] = TRUE
                                #}else{maf[i, 'KEEP'] = FALSE}
                              #}
}
return(maf)
}