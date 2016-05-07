## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

unique(discovery.coding.snps[discovery.coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(discovery.coding.indels[discovery.coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(ph.coding.indels[ph.coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(ph.coding.snps[ph.coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(genes.all[genes.all$gene == "NF2", ]$sample)

master.table <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/mastertable_for_R.txt", stringsAsFactors = F)

## Create pair name lists for calculations
hg.list <- master.table[master.table$Analsysis.Set. == 1, ]$Pair.Name
hg.nf2.mutant.list <- master.table[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel == 1, ]$Pair.Name
hg.nf2.wt.list <- analysis.set[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel == 0, ]$Pair.Name
total.nf2.mutant.list <- master.table[(master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH") & 
                                          master.table$NF2.snp.indel == 1, ]$Pair.Name
total.nf2.wt.list <- master.table[(master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH") & 
                                          master.table$NF2.snp.indel == 0, ]$Pair.Name
ph.list <- master.table[master.table$Cohort == "PH", ]$Pair.Name
ph.table <- master.table[master.table$Cohort == "PH", ]

## order master table for comut

## Hg only
hg.table <- master.table[master.table$Analsysis.Set. == 1, ]
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
hg.unique.snps <- FilterMaf(discovery.coding.snps, hg.list, "Tumor_Sample_Barcode")
x <- table(hg.coding.snps$Tumor_Sample_Barcode)
dimnames(x) <- NULL
x <- x[-13]

y <- table(lg.coding.snps$Tumor_Sample_Barcode)
dimnames(y) <- NULL
y <- y[-(1:39)]
t.test(x, y)

hg.coding.snps <- FilterMaf(discovery.snps, c("Silent", snp.variants),"Variant_Classification")
hg.coding.snps <- FilterMaf(hg.coding.snps, hg.list, "Tumor_Sample_Barcode")

lg.coding.snps <- FilterMaf(ph.snps, c("Silent", snp.variants), "Variant_Classification")



PlotMaf(disc.snindels, "Hugo_Symbol", percent = 5)


## fisher's tests
fisher.test(table(total.table$NF2.snp.indel, total.table$chr22.loss, dnn = c("NF2 status", "chr22 status")))

# check if statistically signficant difference in cohorts between concurrence of nf2/chr22
counts <- NULL
counts <- c(counts, sum(ph.table$NF2.snp.indel.rearrangement == 1 & ph.table$chr22.loss == 1))
counts <- c(counts, sum(hg.table$NF2.snp.indel.rearrangement == 1 & hg.table$chr22.loss == 1))
counts <- c(counts, sum(ph.table$NF2.snp.indel.rearrangement == 0 & ph.table$chr22.loss == 1))
counts <- c(counts, sum(hg.table$NF2.snp.indel.rearrangement == 0 & hg.table$chr22.loss == 1))

fisher.test(matrix(counts,2,2 ))

fisher.test(table(total.table$Chr1.loss, total.table$chr22.loss))

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

## creates list of mutated genes in discovery
gene.list <- disc.snindels[disc.snindels$i_tumor_f > .049, ]
gene.list <- PerSampleMaf(gene.list, "Hugo_Symbol", "Tumor_Sample_Barcode")
gene.list <- ReccurentMaf(gene.list, "Hugo_Symbol")
temp <- sort(table(gene.list$Hugo_Symbol))
gene.table <- dimnames(temp)[[1]]
gene.table <- cbind(gene.table, temp)
gene.table <- as.data.frame(gene.table)
gene.table <- gene.table[order(gene.table[,1]), ]
onc <- gene.list[!duplicated(gene.list$Hugo_Symbol), ]
gene.table <- cbind(gene.table, onc[, 6])
colnames(gene.table)[2:3] <- c("number.times.mutated", "number.times.mutated.cosmic")
