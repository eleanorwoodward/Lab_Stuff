d## Takes a MAF file, runs the exac filter on it, and spits out a tsv file. 

source("C:/Users/Noah/OneDrive/Work/R/Scripts/ExactFilter.R")

validation.snp.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/ValidationSNPs.xls", stringsAsFactors=FALSE, comment.char = "#")

discovery.snps.filtered <- run.exac(maf = snps.discovery)

write.table(validation.snp.filtered, file = "C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/filteredsnp.txt",row.names=F,
          sep = "\t")

