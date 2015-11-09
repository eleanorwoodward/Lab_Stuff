## Noah Greenwald
## Wrangles data into correct format for subsequent analysis

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

indel.variants <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Start_Codon_Del", "Stop_Codon_Del")
snp.variants <- c("De_novo_Start_OutOfFrame", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Start_Codon_SNP")


## Generate differently filtered MAFs for analysis
val.snp.original <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Val829.annotated", 
                                stringsAsFactors=FALSE, comment.char = "#")

val.snp.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/ContEstHigh.annotated",
                           stringsAsFactors = FALSE)

val.snp.no.filter$Tumor_Sample_Barcode <- sapply(val.snp.no.filter$Tumor_Sample_Barcode, PairSetFormat, 2, USE.NAMES = FALSE)

val.snp.no.filter <- run.pon(val.snp.no.filter, -2.5)

val.snp.all.muts <- val.snp.no.filter[val.snp.no.filter$pon_germline == FALSE, ]

val.snp <- FilterMaf(val.snp.no.filter, snp.variants, "Variant_Classification")
val.snp.1.5 <- run.pon(val.snp, -1.5)
val.snp <- val.snp[val.snp$pon_germline == FALSE, ]
val.snp.1.5 <- val.snp.1.5[val.snp.1.5$pon_germline == FALSE, ]



## Filter indels
val.indel.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/ValIndels.annotated", 
                                  stringsAsFactors=FALSE, comment.char = "#")

val.indel.no.filter <- run.exac(val.indel.no.filter, .001)
val.indel.no.filter <- run.esp(val.indel.no.filter)

val.indel.no.filter$Tumor_Sample_Barcode <- sapply(val.indel.no.filter$Tumor_Sample_Barcode, PairSetFormat, 3, USE.NAMES = FALSE)

val.indel.all.muts <- val.indel.no.filter[val.indel.no.filter$esp_germline == FALSE, ]
val.indel.all.muts <- val.indel.no.filter[val.indel.all.muts$germline == FALSE, ]

val.indel <- FilterMaf(val.indel.no.filter, indel.variants, "Variant_Classification")

val.indel <- val.indel[val.indel$germline == FALSE, ]
val.indel <- val.indel[val.indel$esp_germline == FALSE, ]



## Discovery analysis 
disc.snp.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/DiscoverySnps.txt", 
                                 stringsAsFactors=FALSE, comment.char = "#")

disc.indel.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/DiscoveryIndels.txt", 
                                 stringsAsFactors=FALSE, comment.char = "#")

disc.indel <- FilterMaf(disc.indel.no.filter, indel.variants, "Variant_Classification")

disc.snp <- FilterMaf(disc.snp.no.filter, snp.variants, "Variant_Classification")
disc.snp.silent <- FilterMaf(disc.snp.no.filter, c(snp.variants, "Silent"), "Variant_Classification")

recurrences <- c("MEN0030-TumorB", "MEN0042-TumorB", "MEN0042-TumorC", "MEN0042-TumorC", "MEN0045-TumorB", "MEN0045-TumorC", "MEN0045-TumorD", "MEN0045-TumorE",
   "MEN0048-TumorB", "MEN0048-TumorC", "MEN0048-TumorD", "MEN0093-TumorB", "MEN0093-TumorC", "MEN0093-TumorD", "MEN0093-TumorE", "MEN0097-TumorA", 
   "MEN0097-TumorB", "MEN0097-TumorC", "MEN0101-TumorB", "MEN0118-TumorB", "MEN0119-TumorB", "MEN0120-TumorB")

disc.snp.primary <- FilterMaf(disc.snp, recurrences, "Tumor_Sample_Barcode", FALSE)
disc.snp.primary.silent <- FilterMaf(disc.snp.silent, recurrences, "Tumor_Sample_Barcode", FALSE)
disc.indel.primary <- FilterMaf(disc.indel, recurrences, "Tumor_Sample_Barcode", FALSE)

## Combine SNPs and Indels into one maf, keeping relevant columns
silent.snps.disc <- MiniMaf(disc.snp.primary.silent, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
silent.snps.disc[, 6] <- 0
silent.indels.disc <- MiniMaf(disc.indel.primary, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
silent.indels.disc[, 6] <- 1
snindels.disc.silent <- rbind(silent.snps.disc, silent.indels.disc)
names(snindels.disc.silent)[6] <- "Indel"
names(snindels.disc.silent)[5] <- "tumor_f"


short.snps.disc <- MiniMaf(disc.snp.primary, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.snps.disc[, 6] <- 0
short.indels.disc <- MiniMaf(disc.indel.primary, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.indels.disc[, 6] <- 1
snindels.disc <- rbind(short.snps.disc, short.indels.disc)
names(snindels.disc)[6] <- "Indel"
names(snindels.disc)[5] <- "tumor_f"


## Read in Strelka calls if needed
strelka.snp <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/StrelkaSNP.tsv", 
                          stringsAsFactors=FALSE, comment.char = "#")

strelka.snp <- FilterMaf(strelka.snp, c("De_novo_Start_OutOfFrame", "Missense_Mutation", "Nonsense_Mutation", 
                                        "Nonstop_Mutation", "Splice_Site"), "Variant_Classification")

## CNV analysis

sample.list <- c("M109-tumor", "M114-tumor", "M121-tumor", "M133-tumor", "M144-tumor", "M145-tumor", "M147-tumor",  
                 "M148-tumor",
                 "M149-tumor", "M156-tumor", 
                 "M157-tumor", "M159-tumor", "M162-tumor", "M167-tumor", "M168-tumor", "M17-tumor", "M173-tumor", 
                 "M18-tumor",  "M187-tumor", "M188-tumor", "M189-tumor", "M190-tumor", "M191-tumor", "M192-tumor", 
                 "M194-tumor", "M195-tumor", "M196-tumor", "M197-tumor", "M198-tumor", "M2-tumor",  "M5-tumor",
                 "M201-tumor", 
                 "M203-tumor", "M204-tumor", "M205-tumor", "M206-tumor", "M207-tumor", "M208-tumor", "M209-tumor", 
                 "M213-tumor", "M215-tumor", "M216-tumor", "M217-tumor", "M218-tumor", "M219-tumor", "M222-tumor", 
                 "M223-tumor", "M226-tumor", "M227-tumor",  "M23-tumor", "M233-tumor", "M24-tumor",  "M246-tumor", 
                 "M250-tumor", "M26-tumor",  "M264-tumor", "M265-tumor", "M266-tumor", "M269-tumor", "M27-tumor",  
                 "M270-tumor", "M272-tumor", "M43-tumor", "M44-tumor", "M45-tumor",  "M62-tumor",  "M63-tumor",  
                 "M67-tumor", "M71-tumor",  "M73-tumor", "M83-tumor", "M85-tumor", "M92-tumor", "M95-tumor")




cnv.backup <- read.delim("C:/Users/Noah/OneDrive/Work/ccgd/geneCNV.txt", stringsAsFactors = FALSE, header = TRUE)
cnv.filtered <- cnv.backup[cnv.backup$GeneCall != "NormalCopy", ]
cnv.filtered <- cnv.filtered[cnv.filtered$GeneCall != "NormalCopy+", ]
cnv.filtered <- cnv.filtered[cnv.filtered$GeneCall != "NormalCopy-", ]
cnv.backup <- cnv.filtered

likely.somatic <- read.delim("C:/Users/Noah/OneDrive/Work/ccgd/likelysomatic.txt", 
                             stringsAsFactors = FALSE, header = TRUE)

format <- function(cnv) {
    for (i in 1:nrow(cnv)){
    cur <- cnv$Sample[[i]]
    temp <- strsplit(cur, "-")[[1]][1]
    temp <- paste("M", temp, "-tumor", sep = "")
    cnv$Sample[i] <- temp
    }
    cnv
}

mutFormat <- function(cnv) {
  for (i in 1:nrow(cnv)){
    cur <- cnv$tumor_sample_name[[i]]
    temp <- strsplit(cur, "-")[[1]][1]
    temp <- paste("M", temp, "-tumor", sep = "")
    cnv$tumor_sample_name[i] <- temp
  }
  cnv
}

cnv.format <- format(cnv.filtered)
idx <- which(cnv.format$Sample %in% sample.list)
final.cnv <- cnv.format[idx, ]

somatic.format <- mutFormat(likely.somatic)
idx <- which(somatic.format$tumor_sample_name %in% sample.list)
final.somatic <- somatic.format[idx, ]


## interval file
intervals <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Firehose/TargetList.tsv", 
                        stringsAsFactors = FALSE, header = FALSE)

## convert to bed
intervals <- intervals[-(1:88), ]

chr <-intervals[, 1]
start.pos <- sapply(as.numeric(intervals[, 2]), '-', 1)
end.pos <- sapply(as.numeric(intervals[,3]), '-', 1)
targets <- seq_along(chr)
text <- rep("target_", length(targets))
targets <- paste(text, targets, sep = "")
bait.file <- cbind(chr, start.pos, end.pos, targets)
write.table(bait.file, "custombait.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

