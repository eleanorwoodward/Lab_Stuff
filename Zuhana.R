## Quickie For Zuhana. Takes two samples, runs ExAc on them, and returns with filtered. 

folder <- ("C:/Users/Noah/OneDrive/Work/temp")
files <- list.files(folder)

indels.1783 <- read.delim(paste(folder, files[1], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)
indels.46 <- read.delim(paste(folder, files[3], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)
snps.46 <- read.delim(paste(folder, files[4], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)
snps.1783 <- read.delim(paste(folder, files[2], sep = "/"), comment.char = '#', stringsAsFactors = FALSE)

exac.afs <- c(.000005, .00005, .0005, .005, .01, .05)

indel.variants <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Start_Codon_Del", "Stop_Codon_Del")
snp.variants <- c("De_novo_Start_OutOfFrame", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Start_Codon_SNP")

indels.1783.fl <- FilterMaf(indels.1783, indel.variants, "Variant_Classification")
indels.46.fl <- FilterMaf(indels.46, indel.variants, "Variant_Classification")
snps.1783.fl <- FilterMaf(snps.1783, snp.variants, "Variant_Classification")
snps.46.fl <- FilterMaf(snps.46, snp.variants, "Variant_Classification")

barplot(count, names.arg = exac.afs, main = "Indels in Sample 46 removed")

indels.46.filtered <- run.exac(indels.46.fl, ".0001")
indels.1783.filtered <- run.exac(indels.1783.fl, ".0001")
snps.46.filtered <- run.exac(snps.46.fl, ".0001")
snps.1783.filtered <- run.exac(snps.1783.fl, ".0001")


write.table(indels.1783.filtered, file = "C:/Users/Noah/OneDrive/Work/temp/sample1783_indels.txt",row.names=F,
            sep = "\t")

write.table(indels.46.filtered, file = "C:/Users/Noah/OneDrive/Work/temp/sample46_indels.txt",row.names=F,
            sep = "\t")

write.table(snps.1783.filtered, file = "C:/Users/Noah/OneDrive/Work/temp/sample1783_snps.txt",row.names=F,
            sep = "\t")

write.table(snps.46.filtered, file = "C:/Users/Noah/OneDrive/Work/temp/sample46_snps.txt",row.names=F,
            sep = "\t")
