## ESP Filter

## Noah Greenwald and PK Agarwalla
## Checks MAF file against ESP database of known germline mutations
run.esp <- function(maf, af = 0){

# print("Reading in ESP Database. How about a snack break?")
# total.esp <-  read.table("C:/Users/Noah/OneDrive/Work/R/dbs/ESPdb.txt", 
#                          stringsAsFactors=FALSE, fill = TRUE)
# 
# ## Reads in column names from ESP db, as they are not formatted correctly for read.table
# col.names <- scan("C:/Users/Noah/OneDrive/Work/R/dbs/ESPdb.txt", what = "character", n = 150)
# col.names <- col.names[109:141]
# colnames(total.esp) <- col.names[1:31]
# colnames(total.esp)[1] <- "NCBI.Base"
# 
# ## Remove unecessary rows and columns to streamline processing
# # Optional: save a copy of master database
# # original.esp <- total.esp
# total.esp <- total.esp[-1, ]
# total.esp <- total.esp[, c(1,4,8, 12)]
# total.esp <- total.esp[as.integer(total.esp$AvgSampleReadDepth) > 30, ]
# 
# # Format database for comparison. Splits position into chromosme and base, 
# # splits change from c>g format to c, g in separate columns for easy comparison
# # and extracts allelic fraction for overall aggregate population
# pos <- strsplit(total.esp$NCBI.Base, ":")
# change <- strsplit(total.esp$Alleles, ">")
# fraction <- strsplit(total.esp$'MAFinPercent(EA/AA/All)', "/")
# 
# 
# ## Other option is to use gsub to both parse and split into new vector
# ## chr = gsub("([0-9]+):([0-9]+)","\\1",base)
# 
# print("Using the power of wind, air, earth and water to create a better formatted data frame")
# esp <- data.frame( chrom =sapply( pos, "[", 1), 
#                base = sapply( pos, "[", 2), 
#                ref  = sapply( change, "[", 1), 
#                alt = sapply(change, "[", 2), 
#                af  = sapply(fraction,  "[", 3), stringsAsFactors = FALSE
# )
# 
# 
# print("Removing deadweight. Harder. Better. Faster. Stronger.")
# dups <- duplicated(esp)
# esp <- esp[!dups, ]
# write.csv(esp, "C:/Users/Noah/OneDrive/Work/R/dbs/ESP_Lite.csv")
# ## Free memory
# rm(total.esp, pos, change, fraction)
## Above pre-processing only needs to run once.

esp <- read.csv("C:/Users/Noah/Documents/Big Files/dbs/ESP_Lite.csv", stringsAsFactors = FALSE)

idx <- which(esp$base %in% maf$Start_position)
esp <- esp[idx, ]

# Now that ESP is formatted, annotate MAF
# Loop through entire maf  
for (i in 1:nrow(maf)){
  
  ## Prints out every 1,000 cycles
  if (i %% 1000 == 0){
    cat('Checking maf row', i, '\n')
  }
  
  ## Get all matches at current position
  match.idx <- which((esp$chrom == maf$Chromosome[i]) & (esp$base == maf$Start_position[i]) & (esp$af > af))
  
  ## Deal with SNPs first   
  if (maf$Variant_Type[i] == "SNP"){
    
    ## Keeps only those esp hits which have single nucleotide changes
    temp.idx <- which(nchar(esp$ref[match.idx]) & nchar(esp$alt[match.idx]) == 1)
    snp.idx <- match.idx[temp.idx]
    
    ## Checks to make sure SNP at given site represents same basepair change
    tmp.idx <- which(maf$Tumor_Seq_Allele2[i] == esp$alt[snp.idx])
    snp.idx <- snp.idx[tmp.idx]
    
    if (length(snp.idx) > 0){
      maf$esp_POS[i] <- esp$base[snp.idx]
      maf$esp_REF[i] <- esp$ref[snp.idx]
      maf$esp_ALT[i] <- esp$alt[snp.idx]
      maf$esp_AF[i] <- esp$af[snp.idx]
      maf$esp_germline[i] <- TRUE
    }else{
      maf$esp_POS[i] <- NA
      maf$esp_REF[i] <- NA
      maf$esp_ALT[i] <- NA
      maf$esp_AF[i] <- 0
      maf$esp_germline[i] <- FALSE
    }
    
    ## Now look at indels
  } else {
    ## Keep only those exac hits which are indels
    temp.idx <- which(nchar(esp$ref[match.idx]) | nchar(esp$alt[match.idx]) != 1)
    indel.idx <- match.idx[temp.idx]
    
    ## Flag as indel, don't check specifics
    if (length(indel.idx > 1)){
      indel.idx <- indel.idx[1]
    }

    if (length(indel.idx) > 0){
      maf$esp_POS[i] <- esp$base[indel.idx]
      maf$esp_REF[i] <- esp$ref[indel.idx]
      maf$esp_ALT[i] <- esp$alt[indel.idx]
      maf$esp_AF[i] <- esp$af[indel.idx]
      maf$esp_germline[i] <- TRUE
    }else{
      maf$esp_POS[i] <- NA
      maf$esp_REF[i] <- NA
      maf$esp_ALT[i] <- NA
      maf$esp_AF[i] <- 0
      maf$esp_germline[i] <-FALSE
    }
    
  }
}

return(maf)

}