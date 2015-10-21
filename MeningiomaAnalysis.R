## Summary Mutation Analysis

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

##  Take formatted data from Meningioma Wrangler output
snps <- val.snp
indels <- val.indel
snps.nf <- val.snp.all.muts
indels.nf <- val.indel.all.muts

## Combine SNPs and Indels into one maf, keeping relevant columns
short.snps <- MiniMaf(snps, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.snps[, 6] <- 0
short.indels <- MiniMaf(indels, c("Hugo_Symbol", "tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.indels[, 6] <- 1
names(short.snps)[5] <- "tumor_f"
snindels <- rbind(short.snps, short.indels)
names(snindels)[6] <- "Indel"

short.snps.nf <- MiniMaf(snps.nf, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.snps.nf[, 6] <- 0
short.indels.nf <- MiniMaf(indels.nf, c("Hugo_Symbol", "tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.indels.nf[, 6] <- 1
names(short.snps.nf)[5] <- "tumor_f"
snindels.nf <- rbind(short.snps.nf, short.indels.nf)
names(snindels.nf)[6] <- "Indel"

short.snps.15 <- MiniMaf(val.snp.1.5, c("Hugo_Symbol", "i_tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.snps.15[, 6] <- 0
short.indels <- MiniMaf(indels, c("Hugo_Symbol", "tumor_f", "Tumor_Sample_Barcode", "Variant_Classification", "Start_position"))
short.indels[, 6] <- 1
names(short.snps.15)[5] <- "tumor_f"
snindels.15 <- rbind(short.snps.15, short.indels)
names(snindels.15)[6] <- "Indel"
 

## Remove mitochondrial associated genes
mt.bool <- substr(snindels$Hugo_Symbol, 1, 3) == "MT-"
snindels <- snindels[!mt.bool, ]

mt.bool <- substr(snindels.15$Hugo_Symbol, 1, 3) == "MT-"
snindels.15 <- snindels.15[!mt.bool, ]

mt.bool <- substr(snindels.nf$Hugo_Symbol, 1, 3) == "MT-"
snindels.nf <- snindels.nf[!mt.bool, ]


## Divide samples into paired and unpaired. 148, 224, 255, 267, 274, 276 were not included in validation set. 
paired.list <- c("M2-tumor", "M17-tumor", "M26-tumor", "M44-tumor", "M45-tumor", "M133-tumor", "M148-tumor", "M159-tumor", "M203-tumor", "M224-tumor", 
                 "M226-tumor", "M255-tumor", "M267-tumor", "M274-tumor", "M276-tumor")

paired <- FilterMaf(snindels, paired.list, "Tumor_Sample_Barcode")
nonpaired <- FilterMaf(snindels, paired.list, "Tumor_Sample_Barcode", FALSE)

paired.table <- table(paired$Tumor_Sample_Barcode)
nonpaired.table <- table(nonpaired$Tumor_Sample_Barcode)
paired.tbl <- data.frame(paired.table)
nonpaired.tbl <- data.frame(nonpaired.table)

boxplot(paired.tbl$Freq, nonpaired.tbl$Freq, names = c("Paired Samples", "Unpaired Samples"),
        main = "Comparison of mutations per sample for paired and unpaired")

x <- c(nrow(paired) / length(unique(paired$Tumor_Sample_Barcode)), nrow(nonpaired) / length(unique(nonpaired$Tumor_Sample_Barcode)))
barplot(x, names.arg = c("Paired Samples", "Unpaired Samples"), main = "Average number of mutations per sample")

## Number of paired samples vs unpaired samples with mutations
length(unique(paired$Tumor_Sample_Barcode))
length(unique(nonpaired$Tumor_Sample_Barcode))

## Average allelic fraction of paired vs unpaired samples
mean(paired$i_tumor_f)
mean(nonpaired$i_tumor_f)

## Number of mutated genes in paired vs unpaired
length(unique(paired$Hugo_Symbol))
length(unique(nonpaired$Hugo_Symbol))

## Fig 1. Barplot of mutations per sample
PlotMaf(snindels, "Tumor_Sample_Barcode", 2, "Fig 1. SNP + Indel Mutations Per Sample")
abline(h = nrow(snindels) / length(unique(snindels$Tumor_Sample_Barcode)))

length(unique(snindels$Tumor_Sample_Barcode))
sum(snindels$Tumor_Sample_Barcode == "M156-tumor")
sum(snindels$Tumor_Sample_Barcode == "M204-tumor")
sum(snindels$Tumor_Sample_Barcode == "M269-tumor")
sum(snindels$Tumor_Sample_Barcode == "M231-tumor")

## Mutations per gene
snin.unique <- PerSampleMaf(snindels, "Hugo_Symbol")
snindels5 <- ReccurentMaf(snin.unique, "Hugo_Symbol", 4)
snindels10 <- ReccurentMaf(snin.unique, "Hugo_Symbol", 9)
PlotMaf(snindels5, "Hugo_Symbol", 2, "SNP + Indel Mutations in Genes with at least 5 hits")

# Shows difference in mutation calls for differnt filter levels
FilterCutoffMaf(snindels, snindels.15, 4, "PoN -2.5 Cutoff", "PoN -1.5 Cutoff", "Effect of PoN filtering on highly mutated genes")

# Calculate ratio of silent to coding
coding.variants <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Nonstop_Mutation", 
                     "De_Novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")
noncoding.variants <- "Silent"

RatioMaf(snindels.nf, "Variant_Classification", coding.variants, noncoding.variants)

silent.maf <- FilterMaf(snindels.nf, c(coding.variants, noncoding.variants), "Variant_Classification")


small <- FilterMaf(silent.maf, snindels5$Hugo_Symbol, "Hugo_Symbol")
silent <- FilterMaf(small, noncoding.variants, "Variant_Classification")
loud <- FilterMaf(small, coding.variants, "Variant_Classification")
table1 <- table(silent$Hugo_Symbol)
table2 <- table(loud$Hugo_Symbol)
table3 <- EqualizeTable(table2, table1)
barplot(table3, beside = TRUE, main = "Silent vs NonSilent Mutations, multiple per sample", 
        legend.text = c("Coding Mutations", "Silent Mutations"), las = 2)

## Find all mutations with given ratio of silent to coding
pass <- c()
for(i in 1:length(colnames(table3))){
  silent <- table3[2*i]
  coding <- table3[2*i - 1]
  if (silent == 0 | coding / silent >= 5 & coding > 7){
    pass <- c(pass, colnames(table3)[i])
  }
}


## Comparison of frameshift/nonsese vs others
working.maf <- ReccurentMaf(snindels, "Hugo_Symbol", 5)
working.maf <- PerSampleMaf(working.maf, "Hugo_Symbol")
table2 <- table(working.maf$Hugo_Symbol)
table1 <- table(FilterMaf(working.maf, c("Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "Frame_Shift_Del", "De_novo_Start_OutOfFrame"), 
            "Variant_Classification")$Hugo_Symbol)
table3 <- EqualizeTable(table1,table2)
barplot(table3, main = "Frameshift/Nonsense/Splice/DeNovoStart vs Total Coding", 
        legend.text = c("Damaging Mutations", "All coding"), las = 2, beside = TRUE, args.legend = list(x = "topleft"))


## Hotspot finder
snindels11 <- ReccurentMaf(snindels, "Hugo_Symbol", 4)
gene.list <- sort(unique(snindels11$Hugo_Symbol))
mtrx <- matrix(0, 40, length(gene.list))
for(i in 1:length(gene.list)){
  temp <- FilterMaf(snindels11, gene.list[i], "Hugo_Symbol")
  x <- table(temp$Start_position)
  x <- sort(unname(x), decreasing = TRUE)
  if (length(x) > 40){
    x <- x[1:40]
  }
  
  if (length(x) < 40){
    temp <- 40 - length(x)
    x <- c(x, rep(0, temp))
  }
  mtrx[,i] <- x
}

barplot(mtrx, main = "Mutations grouped by base start position", names.arg = gene.list, las = 2)


## Fig 2. Average allelic fraction of mutations in a given sample
snindels.156 <- snindels[snindels$Tumor_Sample_Barcode == "M156-tumor", ]
genelist.156 <- unique(snindels.156$Hugo_Symbol)
avgs.156 <- AverageMaf(snindels.156, genelist.156, "Hugo_Symbol","i_tumor_f")
barplot(avgs.156, names.arg = genelist.156, main = "Fig 2. Average Allelic Fraction of SNPs in M156", las = 2)


## Fig 3. Average allelic fraction of genes mutated at least twice
genlist <- sort(c("CRIPAK", "FRG1", "FRG1B", "NF2", "PRB1", "STK19", "APOBR", "FAM155A"))
avgs <- AverageMaf(snindels5, genlist, "Hugo_Symbol", "tumor_f")
barplot(avgs, names.arg = genlist, main = "Fig 3. Allelic Fraction of SNPs + Indels", las = 2)

## Fig 4. Mutations in NF2 mutated samples
nf2.muts <- FilterMaf(snindels, "NF2", "Hugo_Symbol")
nf2.muts.list <- nf2.muts$Tumor_Sample_Barcode
nf2.muts <- FilterMaf(snindels, nf2.muts.list, "Tumor_Sample_Barcode")
PlotMaf(nf2.muts, "Tumor_Sample_Barcode", 2, "Fig 4. SNPs + Indels in NF2 Mutated Samples")
avg <- nrow(nf2.muts) / length(unique(nf2.muts$Tumor_Sample_Barcode))
avg
abline(h = avg)

## Fig 5. Mutations in nonf2 samples
nonf2.muts <- FilterMaf(snindels, nf2.muts.list, "Tumor_Sample_Barcode", FALSE)
PlotMaf(nonf2.muts, "Tumor_Sample_Barcode", 2, "Fig 5. SNPs + Indels in Non-NF2 Mutated Samples")
avg <- nrow(nonf2.muts) / length(unique(nonf2.muts$Tumor_Sample_Barcode))
abline(h = avg)
avg



## Fig 6a,b. Genes mutated in NF2 vs NonF2 samples
PlotMaf(paired, "Hugo_Symbol", 2, "SNPS + Indels in Paired Samples")
abline(h = nrow(paired) / length(unique(paired$Hugo_Symbol)))

PlotMaf(nonpaired, "Hugo_Symbol", 2, "SNPS + Indels in NonPaired Samples")
abline(h = nrow(nonpaired) / length(unique(nonpaired$Hugo_Symbol)))

x <- FilterMaf(snindels, "NF2", "Hugo_Symbol")
y <- x[, c("Start_position", "Variant_Classification", "Tumor_Sample_Barcode")]
y[order (y[,1]), ]


## CCGD version
somatic <- final.somatic
## Fig 1. Barplot of mutations per sample
PlotMaf(somatic, "tumor_sample_name", 2, "Fig 1. SNP + Indel Mutations Per Sample")
abline(h = nrow(somatic) / length(unique(somatic$tumor_sample_name)))

length(unique(somatic$Tumor_Sample_Barcode))
sum(somatic$tumor_sample_name == "M156-tumor")
sum(somatic$tumor_sample_name == "M204-tumor")
sum(somatic$tumor_sample_name == "M269-tumor")
sum(somatic$tumor_sample_name == "M231-tumor")

## Mutations per gene
somatic6 <- ReccurentMaf(somatic, "BestEffect_Hugo_Symbol", 5)
PlotMaf(somatic6, "BestEffect_Hugo_Symbol", 2, "SNP + Indel Mutations in Genes with at least 2 hits")

somatic10 <- ReccurentMaf(somatic, "BestEffect_Hugo_Symbol", 9)
unique(somatic10$BestEffect_Hugo_Symbol)


## Remove other problem genes from data set
bad.boys <- c("CDC27")
snps <- FilterMaf(snps, bad.boys, "Hugo_Symbol", FALSE)


## Fig 4. Mutations in NF2 mutated samples
nf2.muts <- FilterMaf(somatic, "NF2", "BestEffect_Hugo_Symbol")
nf2.muts.list <- nf2.muts$tumor_sample_name
nf2.muts <- FilterMaf(somatic, nf2.muts.list, "tumor_sample_name")
PlotMaf(nf2.muts, "tumor_sample_name", 2, "Fig 4. SNPs + Indels in NF2 Mutated Samples")
avg <- nrow(nf2.muts) / length(unique(nf2.muts$tumor_sample_name))
avg
abline(h = avg)

## Fig 5. Mutations in nonf2 samples
nonf2.muts <- FilterMaf(somatic, nf2.muts.list, "tumor_sample_name", FALSE)
PlotMaf(nonf2.muts, "tumor_sample_name", 2, "Fig 5. SNPs + Indels in Non-NF2 Mutated Samples")
avg <- nrow(nonf2.muts) / length(unique(nonf2.muts$tumor_sample_name))
abline(h = avg)
avg

## Divide samples into paired and unpaired. 148, 224, 255, 267, 274, 276 were not included in validation set. 
paired.list <- c("M2-tumor", "M17-tumor", "M26-tumor", "M44-tumor", "M45-tumor", "M133-tumor", "M148-tumor", "M159-tumor", "M203-tumor", "M224-tumor", 
                 "M226-tumor", "M255-tumor", "M267-tumor", "M274-tumor", "M276-tumor")

paired <- FilterMaf(somatic, paired.list, "tumor_sample_name")
nonpaired <- FilterMaf(somatic, paired.list, "tumor_sample_name", FALSE)

## Number of paired samples vs unpaired samples with mutations
length(unique(paired$tumor_sample_name))
length(unique(nonpaired$tumor_sample_name))

## Average allelic fraction of paired vs unpaired samples
mean(paired$tumor_f)
mean(nonpaired$tumor_f)

## Number of mutated genes in paired vs unpaired
length(unique(paired$Hugo_Symbol))
length(unique(nonpaired$Hugo_Symbol))

## Fig 6a,b. Genes mutated in NF2 vs NonF2 samples
PlotMaf(paired, "Hugo_Symbol", 2, "SNPS + Indels in Paired Samples")
abline(h = nrow(paired) / length(unique(paired$Hugo_Symbol)))

PlotMaf(nonpaired, "Hugo_Symbol", 2, "SNPS + Indels in NonPaired Samples")
abline(h = nrow(nonpaired) / length(unique(nonpaired$Hugo_Symbol)))







## Copy Number Analysis
cnv <- final.cnv

PlotMaf(cnv, "Sample", 2, "Fig 1. Copy Number Variants per Sample")
abline(h = nrow(cnv) / length(unique(cnv$Sample)))

cnv5 <- ReccurentMaf(cnv, "Gene", 4)
PlotMaf(cnv5, "Gene", 2, "CNVs Per Gene with at least 3 alterations")
abline(h = nrow(cnv5) / length(unique(cnv5$Gene)))

cnv.gain <- FilterMaf(cnv, "gain", "GeneCall")
cnv.loss <- FilterMaf(cnv, "loss", "GeneCall")

PlotMaf(ReccurentMaf(cnv.gain, "Gene"), "Gene", 2, "CN Gain per Gene, at least 2")

PlotMaf(ReccurentMaf(cnv.loss, "Gene", 2), "Gene", 2, "CN Loss per Gene, at least 3")



APOBR.snv <- FilterMaf(val.indel, "APOBR", "Hugo_Symbol")
Nf2.cnv.list <- NF2.cnv$Sample
NF2.cnv <- FilterMaf(cnv, Nf2.cnv.list, "Sample")
Nf2.avg <- nrow(NF2.cnv) / length(unique(NF2.cnv$Sample))
PlotMaf(APOBR.snv, "tumor_sample_name", 2, "Copy Number Variants in NF2 loss samples")
abline(h = Nf2.avg)
Nf2.avg


KMT2C.snv <- FilterMaf(val.snp.no.filter, "KMT2C", "Hugo_Symbol")
PlotMaf(KMT2C.snv, "tumor_sample_name", 2, "Copy Number Variants in NF2 loss samples")






## CNV in NF2 loss samples
NF2.cnv <- FilterMaf(cnv, "NF2", "Gene")
Nf2.cnv.list <- NF2.cnv$Sample
NF2.cnv <- FilterMaf(cnv, Nf2.cnv.list, "Sample")
Nf2.avg <- nrow(NF2.cnv) / length(unique(NF2.cnv$Sample))
PlotMaf(NF2.cnv, "Sample", 2, "Copy Number Variants in NF2 loss samples")
abline(h = Nf2.avg)
Nf2.avg

## CNV in samples without NF2 loss
nonf2.cnv <- FilterMaf(cnv, Nf2.cnv.list, "Sample", FALSE)
nonf2.avg <- nrow(nonf2.cnv) / length(unique(nonf2.cnv$Sample))
PlotMaf(nonf2.cnv, "Sample", 2, "Copy Number Variants in samples without NF2 loss")
abline(h = nonf2.avg)
nonf2.avg
length(unique(nonf2.cnv$Sample))
length(unique(NF2.cnv$Sample))

## Mutations in NF2 loss samples
NF2loss.muts <- FilterMaf(somatic, Nf2.cnv.list, "Tumor_Sample_Barcode")
NF2loss.muts.avg <- nrow(NF2loss.muts) / length(unique(Nf2.cnv.list ))
PlotMaf(NF2loss.muts, "Hugo_Symbol", 2, "SNPs/Indels in NF2 Loss samples")
NF2loss.muts.avg
length(unique(Nf2.cnv.list ))

## Mutations in Non NF2 loss samples
NF2no.loss.muts <- FilterMaf(somatic, Nf2.cnv.list, "Tumor_Sample_Barcode", FALSE)
NF2no.loss.muts.avg <- nrow(NF2no.loss.muts) / length(unique(nonf2.cnv$Sample))
PlotMaf(NF2no.loss.muts, "Hugo_Symbol", 2, "SNPs/Indels in non-NF2 loss samples")
NF2no.loss.muts.avg
length(unique(nonf2.cnv$Sample))

## CN Gain in NF2 mutated samples
Nf2muts.gain <- FilterMaf(cnv.gain, nf2.muts.list, "Sample")
PlotMaf(Nf2muts.gain, "Sample", 2, "Copy Number Gain in Samples with NF2 Mutations")

## Cn Loss in NF2 mutated samples
nf2muts.loss <- FilterMaf(cnv.loss, nf2.muts.list, "Sample")
PlotMaf(nf2muts.loss, "Sample", 2, "Copy number loss in samples with NF2 Mutations")

## Cn gain in non-nf2 mutated samples
nonf2muts.gain <- FilterMaf(cnv.gain, nf2.muts.list, "Sample", FALSE)
PlotMaf(nonf2muts.gain, "Sample", 2, "Copy Number Gain in Samples without NF2 Mutations")

## Cn Loss in non-NF2 mutated samples
nonf2muts.loss <- FilterMaf(cnv.loss, nf2.muts.list, "Sample", FALSE)
PlotMaf(nonf2muts.loss, "Sample", 2, "Copy number loss in samples without NF2 Mutations")

## Samples without Nf2 loss or NF2 mutations
list <- unique(nonf2.cnv$Sample)
nof2 <- FilterMaf(nonf2.muts, list, "Tumor_Sample_Barcode")
unique(nof2$Tumor_Sample_Barcode)
PlotMaf(nof2, "Hugo_Symbol", 2, "Genes mutated in samples without NF2 loss or NF2 mutations")
nof2.cnv.list <- nof2$Tumor_Sample_Barcode
nof2.cnv <- FilterMaf(cnv, nof2.cnv.list, "Gene")
PlotMaf(nof2.cnv, "Gene", 2, "Copy Number Alterations in Genes ")







## SNVs
## Fig 1. Barplot of mutations per sample
PlotMaf(snps, "Tumor_Sample_Barcode", 2, "Fig 1. SNP Mutations Per Sample")

length(unique(snps$Tumor_Sample_Barcode))
length(unique(val.snp.no.filter$Tumor_Sample_Barcode))
sum(snps$Tumor_Sample_Barcode == "M156-tumor")
sum(snps$Tumor_Sample_Barcode == "M204-tumor")
sum(snps$Tumor_Sample_Barcode == "M269-tumor")

## Remove MT-genes from dataset
mt.bool <- substr(snps$Hugo_Symbol, 1, 2) == "MT"
snps <- snps[!mt.bool, ]

## Mutations breakdown for invididual sample
snps.156 <- snps[snps$Tumor_Sample_Barcode == "M156-tumor", ]
PlotMaf(snps.156, "Hugo_Symbol", title = "Fig 2. SNPS per gene in Sample M156")

## Mutations per gene
PlotMaf(snps, "Hugo_Symbol", 2, "SNP Mutations Per Gene")

## Remove other problem genes from data set
bad.boys <- c("CDC27")
snps <- FilterMaf(snps, bad.boys, "Hugo_Symbol", FALSE)

## Fig 2. Average allelic fraction of mutations in a given sample
snps.204 <- snps[snps$Tumor_Sample_Barcode == "M204-tumor", ]
genelist.204 <- unique(snps.204$Hugo_Symbol)
avgs.204 <- AverageMaf(snps.204, genelist.204, "Hugo_Symbol","i_tumor_f")
barplot(avgs.204, names.arg = genelist.204, main = "Fig 2. Average Allelic Fraction of SNPs in M156", las = 2)

## Fig 3. Average allelic fraction of all mutated genes
genlist <- sort(unique(snps$Hugo_Symbol))
avgs <- AverageMaf(snps, genlist, "Hugo_Symbol", "i_tumor_f")
barplot(avgs, names.arg = genlist, main = "Fig 3. Allelic Fraction of SNPs", las = 2)

## Fig 4. Mutations in NF2 mutated samples
nf2.muts <- FilterMaf(snps, "NF2", "Hugo_Symbol")
nf2.muts.list <- nf2.muts$Tumor_Sample_Barcode
nf2.muts <- FilterMaf(snps, nf2.muts.list, "Tumor_Sample_Barcode")
PlotMaf(nf2.muts, "Tumor_Sample_Barcode", 2, "Fig 4. SNPs in NF2 Mutated Samples")

## Fig 5. Mutations in nonf2 samples
nonf2.muts <- FilterMaf(snps, nf2.muts.list, "Tumor_Sample_Barcode", FALSE)
PlotMaf(nonf2.muts, "Tumor_Sample_Barcode", 2, "Fig 5. SNPs in Non-NF2 Mutated Samples")

length(unique(nf2.muts$Tumor_Sample_Barcode))
length(unique(nonf2.muts$Tumor_Sample_Barcode))

## Divide samples into paired and unpaired. 148, 224, 255, 267, 274, 276 were not included in validation set. 9 inc23luded, 8 had mutations
paired.list <- c("M2-tumor", "M17-tumor", "M26-tumor", "M44-tumor", "M45-tumor", "M133-tumor", "M148-tumor", "M159-tumor", "M203-tumor", "M224-tumor", 
                 "M226-tumor", "M255-tumor", "M267-tumor", "M274-tumor", "M276-tumor")

paired <- FilterMaf(snps, paired.list, "Tumor_Sample_Barcode")
nonpaired <- FilterMaf(snps, paired.list, "Tumor_Sample_Barcode", FALSE)

## Number of paired samples vs unpaired samples with mutations
length(unique(paired$Tumor_Sample_Barcode))
length(unique(nonpaired$Tumor_Sample_Barcode))

mean(paired$i_tumor_f)
mean(nonpaired$i_tumor_f)

## Hotspot finder
snps2 <- ReccurentMaf(snps, "Hugo_Symbol")
list2 <- unique(snps2$Hugo_Symbol)
x <- FilterMaf(snps, list2[49], "Hugo_Symbol")
y <- x[, c("Start_position", "Variant_Classification", "Tumor_Sample_Barcode")]
y[order (y[,1]), ]
list2[42]
PlotMaf(y, "Start_position", title = "Mutations per location in CSMD3")




## Indels version
val.indel.no.filter <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/ValIndels.txt", 
                                  stringsAsFactors=FALSE, comment.char = "#")
# val.indel.no.filter <- run.exac(val.indel.no.filter)
# val.indel <- val.indel[which(!val.indel$germline), ]
val.indel <- FilterMaf(val.indel.no.filter, c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                                              "In_Frame_Ins", "Splice_Site"), "Variant_Classification")

for (i in 1:length(val.indel$Hugo_Symbol)){
  if (val.indel$dbSNP_RS[i] != ""){
    if (val.indel$COSMIC_overlapping_mutations[i] != ""){
      val.indel$db_filter[i] <- "flip.flopper"
    }else{
      val.indel$db_filter[i] <- "germline"
    } 
  }else{
    val.indel$db_filter[i] <- "somatic"
  }
}

indels <- val.indel
indels <- indels[indels$db_filter != "germline", ]

## Fig 1. Barplot of mutations per sample
PlotMaf(indels, "Tumor_Sample_Barcode", 2, "Fig 1. Indel Mutations Per Sample")
length(unique(indels$Tumor_Sample_Barcode))

## Mutations per gene
PlotMaf(indels, "Hugo_Symbol", 2, "Indels per Gene")
length(unique(indels$Hugo_Symbol))

indels3 <- ReccurentMaf(indels,"Hugo_Symbol", 2)
sort(unique(indels3$Hugo_Symbol))


## Fig 2. Average allelic fraction of mutations in a given sample
genelist.24 <- unique(indels.24$Hugo_Symbol)
avgs.24 <- AverageMaf(indels.24, genelist.24, "Hugo_Symbol", "tumor_f")
barplot(avgs.24, names.arg = genelist.24, main = "Fig 2. Average Allelic Fraction of Indels in M24", las = 2)

## Fig 3. Average allelic fraction of all mutated genes
genlist <- unique(indels$Hugo_Symbol)
avgs <- AverageMaf(indels, genlist, "Hugo_Symbol", "tumor_f")
barplot(avgs, names.arg = genlist, main = "Fig 3. Allelic Fraction of Indels", las = 2)
mean(indels$tumor_f)

## Fig 4. Mutations in NF2 mutated samples
nf2.muts <- FilterMaf(indels, "NF2", "Hugo_Symbol")
nf2.muts.list <- nf2.muts$Tumor_Sample_Barcode
nf2.muts <- FilterMaf(indels, nf2.muts.list, "Tumor_Sample_Barcode")
PlotMaf(nf2.muts, "Tumor_Sample_Barcode", 2, "Fig 4. Indels in NF2 Mutated Samples")

## Fig 5. Mutations in nonf2 samples
nonf2.muts <- FilterMaf(indels, nf2.muts.list, "Tumor_Sample_Barcode", FALSE)
PlotMaf(nonf2.muts, "Tumor_Sample_Barcode", 2, "Fig 5. Indels in Non-NF2 Mutated Samples")

length(unique(nonf2.muts$Tumor_Sample_Barcode))
length(unique(nf2.muts$Tumor_Sample_Barcode))

## Divide samples into paired and unpaired. 148, 224, 255, 267, 274, 276 were not included in validation set. 
paired.list <- c("M2-tumor", "M17-tumor", "M26-tumor", "M44-tumor", "M45-tumor", "M133-tumor", "M148-tumor", "M159-tumor", "M203-tumor", "M224-tumor", 
                 "M226-tumor", "M255-tumor", "M267-tumor", "M274-tumor", "M276-tumor")

paired <- FilterMaf(indels, paired.list, "Tumor_Sample_Barcode")
nonpaired <- FilterMaf(indels, paired.list, "Tumor_Sample_Barcode", FALSE)

## Number of paired samples vs unpaired samples with mutations
length(unique(paired$Tumor_Sample_Barcode))
length(unique(nonpaired$Tumor_Sample_Barcode))

## Average allelic fraction of paired vs unpaired samples
mean(paired$tumor_f)
mean(nonpaired$tumor_f)

list <- unique(indels3$Hugo_Symbol)
listt <- list[1:3]

# Hotspot finder
for (i in 1:length(list)){
x <- FilterMaf(somatic, "ATRX", "Hugo_Symbol")
y <- x[, c("Start_position", "Variant_Classification", "Tumor_Sample_Barcode")]
print(list[i])
print(y[order (y[,1]), ])
}
table(indels$Hugo_Symbol)

sort(unique(somatic$Hugo_Symbol))

x

