## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

unique(coding.snps[coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(coding.indels[coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)


unique(genes.all[genes.all$gene == "NF2", ]$sample)


## order master table for comut
master.table <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/coMut/mastertable.txt", stringsAsFactors = F)

## Hg only

analysis.set <- master.table[master.table$Analsysis.Set. == 1, ]

analysis.set <- analysis.set[order(-analysis.set$NF2.snp.indel, analysis.set$Grade), ]

## for heatmap generation

write.table(analysis.set[, 2], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/hg.unique.4.14/by_nf2_order.txt", sep = "\t", row.names = F, quote = F)

## Make comut friendly table

comut <- matrix(NA, 9, 40)
comut <- t(analysis.set[, c(2, 5,17, 19, 20)])

write.csv(comut, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/hg.unique.4.14/comut_by_nf2.csv")




## HG + low grade, sorted by grade -> nf2 -> chr

gistic.set <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]

gistic.set <- gistic.set[order(gistic.set$Grade, -gistic.set$NF2.snp.indel, -gistic.set$chr22.loss, -gistic.set$Chr1.loss), ]

write.table(gistic.set[, 2], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/gistic_set_table_order.txt", sep = "\t", row.names = F, quote = F)

comut <- matrix(NA, 9, 40)
comut <- t(gistic.set[, c(2, 5,6,17, 19, 20)])

write.csv(comut, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/total.gistic.comut_grade.csv")

## Hg + low grade, sorted by chr22 -> chr 1 -> grade

gistic.set <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]

gistic.set <- gistic.set[order(-gistic.set$chr22.loss, -gistic.set$Chr1.loss, gistic.set$Grade, -gistic.set$NF2.snp.indel), ]

write.table(gistic.set[, 2], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/gistic_set_by_chr22.txt", sep = "\t", row.names = F, quote = F)

comut <- matrix(NA, 9, 40)
comut <- t(gistic.set[, c(2, 5,6,17, 19, 20)])

write.csv(comut, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/total.gistic.bychr22.csv")



colnames(gistic.set)[18:20]

mutation.counts <- read.delim("C:/Users/Noah/Downloads/HG_unique.patients.counts_and_rates.txt")
new <- cbind(mutation.counts[-41, c(1, 5)], analysis.set)

loss22 <- new[new$chr22.loss == 1 & new$Analsysis.Set. == 1, ]
nf2and22 <- loss22[loss22$NF2.snp.indel == 1, ]
nof2and22 <- loss22[loss22$NF2.snp.indel == 0, ]



