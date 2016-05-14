## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

unique(discovery.coding.snps[discovery.coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(discovery.coding.indels[discovery.coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(ph.coding.indels[ph.coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(ph.coding.snps[ph.coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(genes.all[genes.all$gene == "NF2", ]$sample)

## Read in master table, add column for total mutations
master.table <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/mastertable_for_R.txt", stringsAsFactors = F)


## Create pair name lists for calculations
hg.list <- master.table[master.table$Analsysis.Set. == 1, ]$Pair.Name
hg.nf2.mutant.list <- master.table[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel == 1, ]$Pair.Name
hg.nf2.wt.list <- analysis.set[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel == 0, ]$Pair.Name
total.nf2.mutant.list <- master.table[(master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH") & 
                                          master.table$NF2.snp.indel == 1, ]$Pair.Name
total.nf2.wt.list <- master.table[(master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH") & 
                                          master.table$NF2.snp.indel == 0, ]$Pair.Name
ph.list <- master.table[master.table$Cohort == "PH",]$Pair.Name
ph.table <- master.table[master.table$Cohort == "PH" | (master.table$Cohort == "onc" & master.table$Grade == "I"), ]

## order master table for comut

## Hg only
hg.table <- master.table[master.table$Analsysis.Set. == 1 | (master.table$Cohort == "onc" & master.table$Grade != "I"), ]
hg.table <- hg.table[order(-hg.table$NF2.snp.indel, hg.table$Grade), ]

## for heatmap generation

write.table(hg.table[, 2], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/hg.unique.4.14/by_nf2_order.txt", sep = "\t", row.names = F, quote = F)

## Make comut friendly table

comut <- matrix(NA, 9, 40)
comut <- t(hg.table[, c(2, 5,17, 19, 20)])

write.csv(comut, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/hg.unique.4.14/comut_by_nf2.csv")




## HG + low grade, sorted by grade -> nf2 -> chr22 

total.table <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]
total.list <- total.table$Pair.Name

total.table <- total.table[order(total.table$Grade, -total.table$NF2.snp.indel.rearrangement, -total.table$chr22.loss, -total.table$Chr1.loss), ]

write.table(total.table[, 3], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.vnc/modified_table_order.txt", sep = "\t", row.names = F, quote = F)

comut <- t(total.table[, c(1,12:30)])

write.csv(comut, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/updated_by_grade_csv.csv")

## Hg + low grade, sorted by chr22 -> chr 1 -> grade

total.table <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]

total.table <- total.table[order(-total.table$chr22.loss, -total.table$Chr1.loss, total.table$Grade, -total.table$NF2.snp.indel), ]

write.table(total.table[, 2], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/gistic_set_by_chr22.txt", sep = "\t", row.names = F, quote = F)

comut <- matrix(NA, 9, 40)
comut <- t(gistic.set[, c(2, 5,6,17, 19, 20)])

write.csv(comut, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.414/total.gistic.bychr22.csv")



## Mutation incidence comparison
hg.unique.snps <- FilterMaf(discovery.coding.snps, c(hg.list, "MEN0093G-P2", "MEN0109-P", "MEN0110-P"), "Tumor_Sample_Barcode")
hg.unique.snps <- hg.unique.snps[hg.unique.snps$Tumor_Sample_Barcode > .1, ]
x <- table(hg.unique.snps$Tumor_Sample_Barcode)
dimnames(x) <- NULL
x <- x[-13]


temp <- ph.coding.snps[ph.coding.snps$i_tumor_f > .1, ]
y <- table(temp$Tumor_Sample_Barcode)
dimnames(y) <- NULL
y <- y[-(1:34)]
t.test(x, y)

## generate vals for prism plots
master.table[master.table$Analsysis.Set. == 1 & master.table$Grade == "II" & master.table$Subtype != "Rhabdoid", ]$nonsynoymous.mutations

master.table[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel.rearrangement == 0, ]$nonsynoymous.mutations
master.table[master.table$Cohort == "PH" & master.table$NF2.snp.indel.rearrangement == 1, ]$nonsynoymous.mutations


hg.coding.snps <- FilterMaf(discovery.snps, c("Silent", snp.variants),"Variant_Classification")
hg.coding.snps <- FilterMaf(hg.coding.snps, hg.list, "Tumor_Sample_Barcode")

lg.coding.snps <- FilterMaf(ph.snps, c("Silent", snp.variants), "Variant_Classification")

PlotMaf(disc.snindels, "Hugo_Symbol", percent = 5)


## fisher's tests
fisher.test(table(total.table$NF2.snp.indel, total.table$chr22.loss, dnn = c("NF2 status", "chr22 status")))

# check if statistically signficant difference in cohorts between concurrence of nf2/chr22
counts <- NULL
counts <- c(counts, sum(ph.table$NF2.snp.indel.rearrangement != 0 & ph.table$chr22.loss == 1))
counts <- c(counts, sum(hg.table$NF2.snp.indel.rearrangement == 1 & hg.table$chr22.loss == 1))
counts <- c(counts, sum(ph.table$NF2.snp.indel.rearrangement == 0 & ph.table$chr22.loss == 1))
counts <- c(counts, sum(hg.table$NF2.snp.indel.rearrangement == 0 & hg.table$chr22.loss == 1))

fisher.test(matrix(counts,2,2 ))

fisher.test(table(total.table$Chr1.loss, total.table$chr22.loss))

## for gene set enrichment analysis
total.snindels.2 <- total.snindels[total.snindels$i_tumor_f > .1, ]
total.snindels.2 <- ReccurentMaf(total.snindels.2, "Hugo_Symbol")
write.csv(total.snindels.2, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/genes_mutated_at_least_twice.csv", row.names = F)

disc.snindels.2 <- disc.snindels[disc.snindels$i_tumor_f > .1, ]
disc.snindels.2 <- ReccurentMaf(disc.snindels.2, "Hugo_Symbol")
write.csv(disc.snindels.2, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/genes_mutated_at_least_twice_hg.csv", row.names = F)

tbl <- (table(total.snindels$Hugo_Symbol))
tbl <- sort(tbl)

total.snindels[total.snindels$Hugo_Symbol== "AKT1", ]


## Plot allelic fraction of detected mutations in interesting cases

men039 <- total.snindels[total.snindels$Tumor_Sample_Barcode == "MEN0039G-P", ]
men039.af <- men039[, 2]
men039.af <- men039.af[as.numeric(men039.af) != 0]

men102 <- total.snindels[total.snindels$Tumor_Sample_Barcode == "MEN0102-P", ]
men102.af <- men102[, 2]
men102.af <- men102.af[as.numeric(men102.af) != 0]

men104 <- total.snindels[total.snindels$Tumor_Sample_Barcode == "MEN0104-P", ]
men104.af <- men104[, 2]
men104.af <- men104.af[as.numeric(men104.af) != 0]

men115 <- total.snindels[total.snindels$Tumor_Sample_Barcode == "MEN0115-P", ]
men115.af <- men115[, 2]
men115.af <- men115.af[as.numeric(men115.af) != 0]

men118 <- total.snindels[total.snindels$Tumor_Sample_Barcode == "MEN0118-P1", ]
men118.af <- men118[, 2]
men118.af <- men118.af[as.numeric(men118.af) != 0]

men016 <- total.snindels[total.snindels$Tumor_Sample_Barcode == "MEN0016-pair", ]
men016.af <- men016[, 2]
men016.af <- men016.af[as.numeric(men016.af) != 0]

plot(c(rep(1, length(men039.af)),rep(2, length(men102.af)), rep(3, length(men104.af)), rep(4, length(men118.af)), rep(5, length(men016.af))),
     c(men039.af, men102.af, men104.af, men118.af, men016.af), xlab = c("Men39, Men102, Men104, Men115, men118, men016"), ylab = "Allelic fraction")



## compare rhabdoid features with stuff
counts <- NULL

counts <- c(counts, sum(total.table$Subtype == "Rhabdoid" & total.table$Chr1.loss == 0))
counts <- c(counts, sum(total.table$Grade != "I" & total.table$Subtype != "Rhabdoid" & total.table$Chr1.loss == 0))
counts <- c(counts, sum(total.table$Subtype == "Rhabdoid" & total.table$Chr1.loss == 1))
counts <- c(counts, sum(total.table$Grade != "I" & total.table$Subtype != "Rhabdoid" & total.table$Chr1.loss == 1))
rhab.mtrx <- matrix(counts, nrow = 2)
fisher.test(rhab.mtrx)

