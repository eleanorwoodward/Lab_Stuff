## Runs on LSF farm
#run.esp <- function(maf){
  print("Reading in ESP Database. How about a snack break?")
  total.esp <-  read.table("ESPdb.txt", 
                           stringsAsFactors=FALSE, fill = TRUE)
  
  ## Reads in column names from ESP db, as they are not formatted correctly for read.table
  col.names <- scan("ESPdb.txt", what = "character", n = 150)
  col.names <- col.names[109:141]
  colnames(total.esp) <- col.names[1:31]
  colnames(total.esp)[1] <- "NCBI.Base"
  
  ## Remove unecessary rows and columns to streamline processing
  # Optional: save a copy of master database
  original.esp <- total.esp
  total.esp <- total.esp[-1, ]
  total.esp <- total.esp[, c(1,4,8, 12)]
  total.esp <- total.esp[as.integer(total.esp$AvgSampleReadDepth) > 30, ]
  
  # Format database for comparison. Splits position into chromosme and base, and splits change from c>g format to c, g in separate columns for easy comparison
  pos <- strsplit(total.esp$NCBI.Base, ":")
  change <- strsplit(total.esp$Alleles, ">")
  fraction <- strsplit(total.esp$'MAFinPercent(EA/AA/All)', "/")
  for (i in 1:length(pos)){
    current <- pos[[i]]
    mutation <- change[[i]]
    af <- fraction[[i]]
    total.esp$chrom[i] <- current[1]
    total.esp$base[i] <- current [2]
    total.esp$ref[i] <- mutation[1]
    total.esp$alt[i] <- mutation[2]
    total.esp$af[i] <- af[3]
    if (i %% 100 == 0) print(c("I hate whoever made this data base. Currently doing their job for them on position", i))
  }
  
write.csv(total.esp, "annotatedESPData.csv")
  
#   maf <- val.snp.no.filter
#   idx <- which(total.esp$base %in% maf$Start_position)
#   print("Removing deadweight. Harder. Better. Faster. Stronger.")
#   idx <- which(total.esp$base %in% maf$Start_position)
  