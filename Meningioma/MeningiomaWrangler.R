## Noah Greenwald
## Wrangles data into correct format for subsequent analysis

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")
## for all analysis
indel.variants <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Start_Codon_Del", "Stop_Codon_Del")
snp.variants <- c("De_novo_Start_OutOfFrame", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Start_Codon_SNP")

## For MutationsIndels
discovery.snps.folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutect/Discovery")
discovery.indel.folder <- ("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Indelocator/Discovery")


discovery.indel <- NULL
for (i in 1:length(list.files(discovery.indel.folder))){
    temp <- read.delim(paste(discovery.indel.folder, list.files(discovery.indel.folder)[i], sep = "/"),
                     stringsAsFactors = F, comment.char = "#")
## renames to pair    
    if (nrow(temp) > 0){
        temp[, 16] <-  strsplit(list.files(discovery.indel.folder)[i], ".indel")[[1]][1]
        discovery.indel <- rbind(discovery.indel, temp)
    }
}

discovery.coding.indels <- FilterMaf(discovery.indel, indel.variants, "Variant_Classification")

discovery.snps <- read.delim(paste(discovery.snps.folder, list.files(discovery.snps.folder)[61], sep = "/"),
                             stringsAsFactors = F, comment.char = "#")
discovery.snps[, 16] <- "MEN0120-P1"
discovery.snps <- discovery.snps[, c(1:90, 266,283,290)]
for (i in 1:length(list.files(discovery.snps.folder))){
    temp <- read.delim(paste(discovery.snps.folder, list.files(discovery.snps.folder)[i], sep = "/"),
                      stringsAsFactors = F, comment.char = "#")
    if (ncol(temp) == 272){
        temp[, 16] <-  strsplit(list.files(discovery.snps.folder)[i], ".snp")[[1]][1]
        temp[, 273:290] <- NA
        colnames(temp)[c(283, 290)] <- c("oxoGCut", "i_ffpe_cut")
        discovery.snps <- rbind(discovery.snps, temp[, c(1:90, 266,283,290)])
    }else if (ncol(temp) == 290){
        temp[, 16] <-  strsplit(list.files(discovery.snps.folder)[i], ".ffpe")[[1]][1]
        discovery.snps <- rbind(discovery.snps, temp[, c(1:90, 266,283,290)])
    }else{
        print("hello")
    }
}    

discovery.coding.snps <- FilterMaf(discovery.snps, snp.variants, "Variant_Classification")


mini.disc.snps <- discovery.coding.snps[, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position", 
                                            "COSMIC_total_alterations_in_gene")]

## snowman
# mini.disc.indels <- discovery.coding.indels[, c("Hugo_Symbol", "i_read_depth", "i_allelic_depth", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position",
#                                                 "COSMIC_total_alterations_in_gene")]
# mini.disc.indels[, 2] <- mini.disc.indels[, 3] / mini.disc.indels[, 2]
# mini.disc.indels <- mini.disc.indels[, -3]
# colnames(mini.disc.indels)[2] <- "i_tumor_f"

mini.disc.indels <- discovery.coding.indels[, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position",
                                                "COSMIC_total_alterations_in_gene")]

disc.snindels <- rbind(mini.disc.indels, mini.disc.snps)
orphans <- FilterMaf(disc.snindels, c("MEN0093G-P2", "MEN0109-P", "MEN0110-P"), "Tumor_Sample_Barcode")
disc.snindels <- FilterMaf(disc.snindels, hg.list[!is.na(hg.list)], "Tumor_Sample_Barcode")


ph.snps.folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutect/LG")
ph.indel.folder <- ("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Indelocator/LG")

ph.indel <- NULL
for (i in 1:length(list.files(ph.indel.folder))){
    temp <- read.delim(paste(ph.indel.folder, list.files(ph.indel.folder)[i], sep = "/"),
                       stringsAsFactors = F, comment.char = "#")
    if (nrow(temp) > 0){
        temp[, 16] <-  strsplit(list.files(ph.indel.folder)[i], ".indel")[[1]][1]
        ph.indel <- rbind(ph.indel, temp)
    }
}

ph.coding.indels <- FilterMaf(ph.indel, indel.variants, "Variant_Classification")

ph.snps <- read.delim(paste(ph.snps.folder, list.files(ph.snps.folder)[63], sep = "/"),
                             stringsAsFactors = F, comment.char = "#")
ph.snps[, 16] <- "MENex006-pair"
ph.snps <- ph.snps[, c(1:90, 266,283,290)]
for (i in 1:length(list.files(ph.snps.folder))){
    temp <- read.delim(paste(ph.snps.folder, list.files(ph.snps.folder)[i], sep = "/"),
                       stringsAsFactors = F, comment.char = "#")
    if (ncol(temp) == 272){
        temp[, 16] <-  strsplit(list.files(ph.snps.folder)[i], ".snp")[[1]][1]
        temp[, 273:290] <- NA
        colnames(temp)[c(283, 290)] <- c("oxoGCut", "i_ffpe_cut")
        ph.snps <- rbind(ph.snps, temp[, c(1:90, 266,283,290)])
    }else if (ncol(temp) == 290){
        temp[, 16] <-  strsplit(list.files(ph.snps.folder)[i], ".ffpe")[[1]][1]
        ph.snps <- rbind(ph.snps, temp[, c(1:90, 266,283,290)])
    }else{
        print("hello")
    }
}

ph.coding.snps <- FilterMaf(ph.snps, snp.variants, "Variant_Classification")


## create all low grade mutations
mini.ph.snps <- ph.coding.snps[, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position")]
mini.ph.indels <- ph.coding.indels[, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position")]

# for snowman
#  mini.ph.indels <- ph.coding.indels[, c("Hugo_Symbol", "i_read_depth", "i_allelic_depth", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position")]
# mini.ph.indels[, 2] <- mini.ph.indels[, 3] / mini.ph.indels[, 2]
# mini.ph.indels <- mini.ph.indels[, -3]
# colnames(mini.ph.indels)[2] <- "i_tumor_f"

ph.snindels <- rbind(mini.ph.indels, mini.ph.snps)
total.snindels <- rbind(disc.snindels[, -6], ph.snindels)
total.snindels <- rbind(total.snindels, orphans[, -6])

total.coding.snps <- rbind(mini.ph.snps, mini.disc.snps[, -6])

## Load rearrangement data sets

## read in high grade
discovery.rearrangements.folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/Discovery/rearrangements/")
discovery.rearrangements <- NULL
for (i in 1:length(list.files(discovery.rearrangements.folder))){
    temp <- read.csv(paste(discovery.rearrangements.folder, list.files(discovery.rearrangements.folder)[i], sep = "/"),
                     stringsAsFactors = F)
    temp[, 28] <-  strsplit(list.files(discovery.rearrangements.folder)[i], ".csv")[[1]]
    colnames(temp)[28] <- "Sample"
    colnames(temp)[29] <- "vcf.info"
    discovery.rearrangements <- rbind(discovery.rearrangements, temp)    
}

## Keep rearrangements only with likely somatic score

discovery.rearrangements <- discovery.rearrangements[discovery.rearrangements$somatic_lod > 4, ]
discovery.rearrangements[, c(39, 40)] <- NA



## read in low grade
ph.rearrangements.folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/LG/rearrangements/")
ph.rearrangements <- NULL
for (i in 1:length(list.files(ph.rearrangements.folder))){
    temp <- read.csv(paste(ph.rearrangements.folder, list.files(ph.rearrangements.folder)[i], sep = "/"),
                     stringsAsFactors = F)
    temp[, 28] <-  strsplit(list.files(ph.rearrangements.folder)[i], ".csv")[[1]]
    colnames(temp)[28] <- "Sample"
    colnames(temp)[29] <- "vcf.info"
    ph.rearrangements <- rbind(ph.rearrangements, temp)    
}

## Keep rearrangements only with likely somatic score

ph.rearrangements <- ph.rearrangements[ph.rearrangements$somatic_lod > 4, ]
ph.rearrangements[, c(39, 40)] <- NA

ph.rearrangements[ph.rearrangements$gene1 == "NF2" | ph.rearrangements$gene2 == "NF2", 27:36]


## Import Copy Number Data, order same as master table
gistic.calls <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/broad_values_by_arm (3).txt", stringsAsFactors = F)
gistic.calls <- cbind(gistic.calls[, -(2:12)], gistic.calls[, 2:12])
arms <- gistic.calls[, 1]
gistic.calls <- t(gistic.calls[, -1])
colnames(gistic.calls) <- arms


## Bare minimum validation processing

validation.snps.folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutect/Validation")
validation.snps <- NULL
for (i in 1:length(list.files(validation.snps.folder))){
    temp <- read.delim(paste(validation.snps.folder, list.files(validation.snps.folder)[i], sep = "/"),
                     stringsAsFactors = F, comment.char = "#")
    print(i)
    validation.snps <- rbind(validation.snps, temp)    
}


validation.coding.snps <- FilterMaf(validation.snps, snp.variants, "Variant_Classification")

validation.coding.snps <- run.pon(validation.coding.snps, -1.5)
validation.coding.snps <- run.exac(validation.coding.snps, .0001)
validation.coding.snps <- run.esp(validation.coding.snps, .0001)

validation.filtered.coding.snps <- validation.coding.snps[validation.coding.snps$pon_germline == F &
                                                                       validation.coding.snps$germline == F & validation.coding.snps$esp_germline == F, ]

validation.indels.folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Indelocator/Validation")
validation.indels <- NULL
for (i in 1:length(list.files(validation.indels.folder))){
    temp <- read.delim(paste(validation.indels.folder, list.files(validation.indels.folder)[i], sep = "/"),
                       stringsAsFactors = F, comment.char = "#")
    validation.indels <- rbind(validation.indels, temp)    

}


validation.coding.indels <- FilterMaf(validation.indels, indel.variants, "Variant_Classification")

validation.coding.indels <- run.pon(validation.coding.indels, -1.5)
validation.coding.indels <- run.exac(validation.coding.indels, .0001)
validation.coding.indels <- run.esp(validation.coding.indels, .0001)
validation.filtered.coding.indels <- validation.coding.indels[validation.coding.indels$germline == F &
                                                                           validation.coding.indels$pon_germline == F & validation.coding.indels$esp_germline == F, ]

## combine validation mutations
mini.val.snps <- validation.filtered.coding.snps[, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position")]
mini.val.indels <- validation.filtered.coding.indels[, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position")]
ccgd.snindels <- rbind(mini.val.snps, mini.val.indels)


for (i in 1:nrow(ccgd.snindels)){
    cur <- ccgd.snindels$Tumor_Sample_Barcode[[i]]
    temp <- as.numeric(strsplit(cur, "-")[[1]][2])
    temp <- formatC(temp, width = 3, flag = 0)
    temp <- paste("MG-", temp, "-tumor", sep = "")
    ccgd.snindels$Tumor_Sample_Barcode[[i]] <- temp
}

ccgd.snindels <- FilterMaf(ccgd.snindels, unique(master.table[master.table$Cohort %in% c("ccgd.lg", "ccgd.hg", "ccgd.tbd"), ]$Pair.Name), "Tumor_Sample_Barcode")
val.snindels <- FilterMaf(ccgd.snindels, master.table[master.table$Cohort %in% c("ccgd.hg", "ccgd.tbd"), ]$Pair.Name,"Tumor_Sample_Barcode" )



## original stuff

## Generate differently filtered MAFs for analysis
val.snp.original <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Val829.annotated", 
                                stringsAsFactors=FALSE, comment.char = "#")

val.snp.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/ContEstHigh.annotated",
                           stringsAsFactors = FALSE)

val.snp.no.filter$Tumor_Sample_Barcode <- sapply(val.snp.no.filter$Tumor_Sample_Barcode, PairSetFormat, 2, USE.NAMES = FALSE)
temp1 <- FilterMaf(val.snp.no.filter, snp.variants, "Variant_Classification")

val.snp.no.filter <- run.pon(val.snp.no.filter, -2.5)
val.snp.no.filter <- run.exac(val.snp.no.filter, .0001)

val.snp.all.muts <- val.snp.no.filter[val.snp.no.filter$pon_germline == FALSE, ]

val.snp <- FilterMaf(val.snp.no.filter, snp.variants, "Variant_Classification")
val.snp.1.5 <- run.pon(val.snp, -1.5)
val.snp <- val.snp[val.snp$pon_germline == FALSE, ]
val.snp.1.5 <- val.snp.1.5[val.snp.1.5$pon_germline == FALSE, ]



## Filter indels
val.indel.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/ValIndels.annotated", 
                                  stringsAsFactors=FALSE, comment.char = "#")
temp2 <- FilterMaf(val.indel.no.filter, indel.variants, "Variant_Classification")

val.indel.no.filter <- run.exac(val.indel.no.filter, .0001)

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

