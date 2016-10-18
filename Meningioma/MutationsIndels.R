## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

unique(discovery.coding.snps[discovery.coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(discovery.coding.indels[discovery.coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(ph.coding.indels[ph.coding.indels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(ph.coding.snps[ph.coding.snps$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)

unique(genes.all[genes.all$gene == "NF2", ]$sample)

## Read in master table, add column for total mutations
master.table <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/mastertable_for_R.txt", stringsAsFactors = F)


## GSEA
akt1.pathway.members <- c("FOXO4", "PDPK1", "FOXO1", "CASP9", "NR4A1",  "GSK3B" , "CDKN1B", "AKT1S1", "GSK3A" ,"MAPKAP1" , "AKT3",
                          "MLST8", "CDKN1A" ,"MDM2", "CHUK" ,"CREB1", "FOXO3", 	"AKT2" ,"BAD", "RPS6KB2", 
                          "GH1", "PTEN", "PIK3CA", "TSC2", "MTOR", "RICTOR", "IRS1", "AKT1")

swi.snf.members <- c("ARID1A", "ARID1B", "SMARCB1", "SMARCA2", "SMARCA4", "SMARCE1", "ARID2", "BRD7", "PBRM1",
                     "KMT2A", "KDM6A", "KDM5C", "KDM4E", "KDM6B")

hedgehog.members <- c("SHH", "GLI2", "GLI3", "KIF7", "STK36", "ADRBK1", "SAP18", "IHH", "PTCH1", 
                      "GLI1", "ARNTL", "DHH", "SUFU", "SMO", "SIN3A", "PTCH2", "CTNNB1", "FOXC1", "PTCH1", "FOXC1")

bwh.snindels <- rbind(total.snindels, ccgd.snindels)
bwh.snindels <- bwh.snindels[bwh.snindels$i_tumor_f > .09, ]
bwh.paired.sample.list <- master.table[master.table$Paired.normal == 1 & (master.table$Analsysis.Set. == 1 | master.table$Cohort %in% c("ccgd.lg", "PH", "ccgd.hg", "ccgd.tbd")), ]$Pair.Name
bwh.sample.list <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort %in% c("ccgd.lg", "PH", "ccgd.hg", "ccgd.tbd"), ]$Pair.Name
bwh.sample.list <- bwh.sample.list[!is.na(bwh.sample.list)]
hg.sample.list <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort %in% c("ccgd.hg", "ccgd.tbd"), ]$Pair.Name
lg.sample.list <- master.table[master.table$Cohort %in% c("ccgd.lg", "PH"), ]$Pair.Name

bwh.snindels <- FilterMaf(bwh.snindels, bwh.sample.list, "Tumor_Sample_Barcode")
bwh.snindels <- PerSampleMaf(bwh.snindels, "Hugo_Symbol")
bwh.paired.snindels <- FilterMaf(bwh.snindels, bwh.paired.sample.list, "Tumor_Sample_Barcode")

hg.snindels <- FilterMaf(bwh.snindels, hg.sample.list, "Tumor_Sample_Barcode")
hg.snindels.2 <- ReccurentMaf(hg.snindels, "Hugo_Symbol")
hg.snindels.3 <- ReccurentMaf(hg.snindels, "Hugo_Symbol", 2)
hg.paired.snindels <- FilterMaf(hg.snindels, bwh.paired.sample.list, "Tumor_Sample_Barcode")
hg.paired.snindels.2 <- ReccurentMaf(hg.paired.snindels, "Hugo_Symbol")
lg.snindels <- FilterMaf(bwh.snindels, lg.sample.list, "Tumor_Sample_Barcode")
lg.paired.snindels <- FilterMaf(lg.snindels, bwh.paired.sample.list, "Tumor_Sample_Barcode")

bwh.snindels.min2 <- ReccurentMaf(bwh.snindels, "Hugo_Symbol")
bwh.snindels.min3 <- ReccurentMaf(bwh.snindels, "Hugo_Symbol", 2)
bwh.snindels.min5 <- ReccurentMaf(bwh.snindels, "Hugo_Symbol", 4)
write.csv(bwh.snindels, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/bwh.normals.atleast2.csv")
write.csv(bwh.snindels.min3, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/bwh.atleast3.csv")
write.csv(bwh.snindels.min5, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/bwh.atleast5.csv")
write.csv(hg.paired.snindels.2, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/hg.paired.2.csv")
write.csv(hg.snindels.3, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/hg.3.csv")


bwh.table <- master.table[master.table$Pair.Name %in% bwh.sample.list, ]
bwh.table <- bwh.table[order(-as.numeric(bwh.table$Medium.Grade), bwh.table$chr22.loss, bwh.table$NF2.snp.indel.rearrangement, 
                             bwh.table$TRAF7, bwh.table$KLF4, bwh.table$AKT1, bwh.table$SMO,  decreasing = T), ]



bwh.table <- cbind(bwh.table, 0)
colnames(bwh.table)[69:72] <- c("mTOR", "hedgehog", "Cell.cycle", "Chromatin")

mTOR <- c("AKT1", "TSC2", "PIK3CA")
hedgehog <- c("CTNNB1", "SUFU", "SMO")
cell.cycle <- c("CDKN2C", "ANAPC2", "CDC27", "TP53", "CHEK2")
Chromatin <- c("KMT2D", "ARID1A", "SMARCA2", "ARID1B", "SMARCB1", "KMT2C", "PRDM9", "ARID2")

target <- swi.snf.members

## fill in pathway mutations
for (i in 1:length(target)){
    bwh.table[bwh.table$Pair.Name %in% bwh.snindels[bwh.snindels$Hugo_Symbol == target[i],]$Tumor_Sample_Barcode, 72] <- target[i]
}

bwh.table <- bwh.table[order(bwh.table$Simple.Histopath.Grade, -bwh.table$NF2.snp.indel.rearrangement, 
                             bwh.table$mTOR == "0", bwh.table$hedgehog == "0", bwh.table$Cell.cycle == "0", bwh.table$Chromatin == "0"), ]


write.csv(bwh.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/BWH.Comut/bhw.all_members.csv")


## check overlap between genes and NF2
bwh.snindels[bwh.snindels$Tumor_Sample_Barcode %in% bwh.snindels[bwh.snindels$Hugo_Symbol == "RICTOR", ]$Tumor_Sample_Barcode, ]

## for gene set enrichment analysis
total.snindels.2 <- total.snindels[total.snindels$i_tumor_f > .1, ]
total.snindels.2 <- ReccurentMaf(total.snindels.2, "Hugo_Symbol")
write.csv(total.snindels.2, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/genes_mutated_at_least_twice.csv", row.names = F)

disc.snindels.2 <- disc.snindels[disc.snindels$i_tumor_f > .1, ]
disc.snindels.2 <- PerSampleMaf(disc.snindels.2, "Hugo_Symbol")
disc.snindels.2 <- ReccurentMaf(disc.snindels.2, "Hugo_Symbol")
write.csv(disc.snindels.2, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/genes_mutated_at_least_twice_disc.csv", row.names = F)




## statistics for paper


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
##
##

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


## Hg + Lg, with angiomatous/rhaboid separated

total.table <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]

total.table <- total.table[order(total.table$Heatmap.Grade, -total.table$chr22.loss), ]

write.table(total.table[, 3], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Heatmap/finalized_clustering_order.txt",
            sep = "\t", row.names = F, quote = F)

write.csv(total.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Heatmap/finalized_clustering_order_comut.txt")

disc.table <- master.table[master.table$Analsysis.Set. == 1, ]
disc.table <- disc.table[!is.na(disc.table$Analsysis.Set.), ]


## Massive plot: multiple grades
massive.table <- master.table[master.table$Master.Cohort == 1, ]
massive.table <- massive.table[order(massive.table$Grade, -massive.table$chr22.loss, -massive.table$NF2.snp.indel.rearrangement,
                                     -massive.table$TRAF7, -massive.table$KLF4, -massive.table$AKT1, -massive.table$SMO), ]
write.csv(massive.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Master.Comut/master.comut.csv")

## Simple grades, missing data compacted: sort by TRAF7 plotting
massive.table <- massive.table[order(massive.table$Simple.Grade, -massive.table$chr22.loss, -massive.table$NF2.snp.indel.rearrangement,
                                     -massive.table$TRAF7, -massive.table$KLF4, -massive.table$AKT1, -massive.table$SMO), ]
write.csv(massive.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Master.Comut/master.comut.new_yale.csv")

## all grades, unclear at end
massive.table <- massive.table[order(massive.table$Grade, -massive.table$chr22.loss, -massive.table$NF2.snp.indel.rearrangement,
                                     -massive.table$TRAF7, -massive.table$KLF4, -massive.table$AKT1, -massive.table$SMO), ]
write.csv(massive.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Master.Comut/master.comut.v3_molec.path.comp.csv")

## all grades, location included
massive.table <- massive.table[order(massive.table$Simple.Grade, -massive.table$chr22.loss, -massive.table$NF2.snp.indel.rearrangement,
                                     -massive.table$TRAF7, -massive.table$KLF4, -massive.table$AKT1, -massive.table$SMO), ]
write.csv(massive.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Master.Comut/master.comut.v4_verticle_location.csv")

## all grades, histopath grade used
massive.table <- massive.table[order(massive.table$Simple.Histopath.Grade, -massive.table$chr22.loss, -massive.table$NF2.snp.indel.rearrangement,
                                     -massive.table$TRAF7, -massive.table$KLF4, -massive.table$AKT1, -massive.table$SMO), ]
write.csv(massive.table, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heatmap.Comut/Master.Comut/master.comut.updateder.csv")


## p values for major mutation classes
fisher.test(table(massive.table$chr22.loss, massive.table$TKAS))
fisher.test(table(massive.table$NF2.snp.indel.rearrangement, massive.table$TKAS))
fisher.test(table(massive.table$Simple.Grade, massive.table$chr22.loss))
fisher.test(table(massive.table$Simple.Grade, massive.table$TKAS))
fisher.test(table(massive.table$Simple.Grade, massive.table$SMO))
fisher.test(table(massive.table$Simple.Grade, massive.table$AKT1))
fisher.test(table(massive.table$Simple.Grade, massive.table$KLF4))
fisher.test(table(massive.table$Simple.Grade, massive.table$TRAF7))
nrow(massive.table[massive.table$Simple.Grade == 2 & massive.table$KLF4 == 1, ])

## comparison of non-nf2 driver events controlling for NF2
fisher.test(table(massive.table[massive.table$chr22.loss == F & massive.table$NF2.snp.indel.rearrangement == F, ]$Simple.Grade, 
                  massive.table[massive.table$chr22.loss == F & massive.table$NF2.snp.indel.rearrangement == F, ]$TKAS))

## mutual exclusivity of NF2 and non-nf2 alterations, even among samples with chr 22 loss
fisher.test(table(massive.table[massive.table$chr22.loss == T,]$NF2.snp.indel.rearrangement, massive.table[massive.table$chr22.loss == T, ]$TKAS))

## mutual exclusivity of TKAS and NF2 in high grades
fisher.test(table(massive.table[massive.table$Simple.Grade == 2, ]$chr22.loss, massive.table[massive.table$Simple.Grade == 2, ]$TKAS))

# p value for difference in NF2 inactivation rates

fisher.test(table(disc.snindels$Hugo_Symbol == "NF2", disc.snindels$Variant_Classification %in% c("Frame_Shift_Del", "Nonsense_Mutation","Splice_Site")))

table(disc.snindels.clean$Hugo_Symbol)[order(table(disc.snindels.clean$Hugo_Symbol), decreasing = T)][1:62]



## investigate suspicious samples

fishy <- massive.table[massive.table$chr22.loss == F & massive.table$NF2.snp.indel.rearrangement == T, 1:10]
samples <- c("MEN_PH_LG_59-pair","MEN_PH_LG_60-pair", "MEN_PH_LG_77-pair", "MEN_PH_LG_85-pair", "MG-133-tumor")

ph.snps[ph.snps$Tumor_Sample_Barcode %in% samples[4], c(1,5,9,10, 16, 91)]
ph.indel[ph.indel$Tumor_Sample_Barcode %in% samples[4], c(1,5,9,10, 16, 268)]

validation.snps[validation.snps$Tumor_Sample_Barcode %in% samples[5], c(1,5,9,10, 16, 266)]
validation.indels[validation.indels$Tumor_Sample_Barcode %in% samples[5], c(1,5,9,10, 16, 268)]

sum(massive.table$Simple.Grade == "II" & (massive.table$NF2.snp.indel.rearrangement == T | massive.table$chr22.loss == T))
sum(massive.table$Simple.Grade == "II" & massive.table$TKAS == 1)



## generate csv with all mutation comparisons for master table
hg.snps <- discovery.duplicate.snps
hg.snps.low.af <- hg.snps[hg.snps$i_tumor_f < .1, ]
hg.snps.coding <- hg.snps[hg.snps$Variant_Classification != "Silent", ]
hg.snps.coding.high.af <- hg.snps.coding[hg.snps.coding$i_tumor_f > .09999, ]

hg.names <- sort(unique(hg.snps$Tumor_Sample_Barcode))
all.muts <- as.numeric(table(hg.snps$Tumor_Sample_Barcode))
all.muts.low.af <-as.numeric(table(hg.snps.low.af$Tumor_Sample_Barcode))
coding.muts.high.af <- as.numeric(table(hg.snps.coding.high.af$Tumor_Sample_Barcode))

## same for low-grade
lg.snps <- ph.snps[ph.snps$Variant_Classification %in% c(snp.variants, "Silent"), ]
lg.snps <- lg.snps[lg.snps$Tumor_Sample_Barcode %in% sort(unique(lg.snps$Tumor_Sample_Barcode))[-(1:39)], ]
lg.snps.low.af <- lg.snps[lg.snps$i_tumor_f < .1, ]
lg.snps.coding <- lg.snps[lg.snps$Variant_Classification != "Silent", ]
lg.snps.coding.high.af <- lg.snps.coding[lg.snps.coding$i_tumor_f > .09999, ]

lg.names <- sort(unique(lg.snps$Tumor_Sample_Barcode))
all.muts.lg <- as.numeric(table(lg.snps$Tumor_Sample_Barcode))
all.muts.low.af.lg <-as.numeric(table(lg.snps.low.af$Tumor_Sample_Barcode))
coding.muts.high.af.lg <- as.numeric(table(lg.snps.coding.high.af$Tumor_Sample_Barcode))


## MEN0045-P3 has no low allelic fraction mutations
## MEN0045-P4 has no coding snps
## MEN0014 has no low.af muts

comparison.csv <- data.frame(hg.names, all.muts, c(all.muts.low.af, 0), c(coding.muts.high.af, 0))
write.csv(comparison.csv, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/comparison.csv")

comparison.csv.lg <- data.frame(lg.names, all.muts.lg, c(all.muts.low.af.lg, 0), coding.muts.high.af.lg)
write.csv(comparison.csv.lg, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/comparison.lg.csv")

## Mutation incidence comparison
hg.unique.snps <- FilterMaf(discovery.coding.snps, c(hg.list, "MEN0093G-P2", "MEN0109-P", "MEN0110-P"), "Tumor_Sample_Barcode")
hg.unique.snps <- hg.unique.snps[hg.unique.snps$i_tumor_f < .1, ]
x <- table(hg.unique.snps$Tumor_Sample_Barcode)
dimnames(x) <- NULL
x <- x[-13]


temp <- ph.coding.snps[ph.coding.snps$i_tumor_f > .1, ]
y <- table(ph.coding.snps$Tumor_Sample_Barcode)
dimnames(y) <- NULL
y <- y[-(1:38)]
t.test(x, y)

## calculations for signficant difference in mutation rates
t.test(total.table[total.table$Simple.Grade == 1, ]$nonsynoymous.mutations.high.af,
       total.table[total.table$Simple.Grade == 2, ]$nonsynoymous.mutations.high.af )

rad.table <- master.table[master.table$Radiation.analysis.set == 1, ]
rad.table <- rad.table[!is.na(rad.table$Radiation.analysis.set), ]

wilcox.test(rad.table[rad.table$XRT..preop.radiation == 1, ]$nonsynoymous.mutations.high.af, 
       rad.table[rad.table$XRT..preop.radiation == 0, ]$nonsynoymous.mutations.high.af)

wilcox.test(rad.table[rad.table$XRT..preop.radiation == 1, ]$total.muts, 
       rad.table[rad.table$XRT..preop.radiation == 0, ]$total.muts)

wilcox.test(disc.table[disc.table$NF2.snp.indel.rearrangement == 1, ]$nonsynoymous.mutations.high.af,
       disc.table[disc.table$NF2.snp.indel.rearrangement == 0, ]$nonsynoymous.mutations.high.af, )

wilcox.test(disc.table[disc.table$NF2.snp.indel.rearrangement == 1, ]$nonsynoymous.mutations.high.af,
            disc.table[disc.table$NF2.snp.indel.rearrangement == 0, ]$nonsynoymous.mutations.high.af)

## copy number

wilcox.test(disc.table[disc.table$NF2.snp.indel.rearrangement == 1, ]$percent.disruption,
            disc.table[disc.table$NF2.snp.indel.rearrangement == 0, ]$percent.disruption)

wilcox.test(disc.table[disc.table$chr22.loss == 1, ]$percent.disruption,
            disc.table[disc.table$chr22.loss == 0, ]$percent.disruption)

wilcox.test(disc.table[disc.table$NF2.snp.indel.rearrangement == 1, ]$percent.disruption,
            disc.table[disc.table$NF2.snp.indel.rearrangement == 0, ]$percent.disruption)

wilcox.test(total.table[total.table$Simple.Grade == 1, ]$percent.disruption,
            total.table[total.table$Simple.Grade == 2, ]$percent.disruption)

wilcox.test(disc.table[disc.table$XRT..preop.radiation == 1, ]$percent.disruption,
            disc.table[disc.table$XRT..preop.radiation == 0, ]$percent.disruption)
## rearrangements

wg.table <- master.table[master.table$Sequencing == "WGS", ]

fisher.test(table(wg.table$complex.event, wg.table$chr22.loss))
wilcox.test(wg.table[wg.table$NF2.snp.indel.rearrangement == 1, ]$rearrangement.burden, wg.table[wg.table$NF2.snp.indel.rearrangement == 0, ]$rearrangement.burden )


fisher.test(table(wg.table$complex.event, wg.table$XRT..preop.radiation))
wilcox.test(wg.table[wg.table$XRT..preop.radiation == 1, ]$rearrangement.burden, wg.table[wg.table$XRT..preop.radiation == 0, ]$rearrangement.burden )


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


## generate vals for prism plots
master.table[master.table$Analsysis.Set. == 1 & master.table$Grade == "II" & master.table$Subtype != "Rhabdoid", ]$nonsynoymous.mutations

master.table[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel.rearrangement == 0, ]$nonsynoymous.mutations
master.table[master.table$Cohort == "PH" & master.table$NF2.snp.indel.rearrangement == 1, ]$nonsynoymous.mutations


hg.coding.snps <- FilterMaf(discovery.snps, c("Silent", snp.variants),"Variant_Classification")
hg.coding.snps <- FilterMaf(hg.coding.snps, hg.list, "Tumor_Sample_Barcode")

lg.coding.snps <- FilterMaf(ph.snps, c("Silent", snp.variants), "Variant_Classification")

PlotMaf(disc.snindels, "Hugo_Symbol", percent = 5)






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


## generate mutation counts for dotplot
pancan <- read.delim("C:/Users/Noah/Documents/Big Files/dbs/compact.data.v3.maf", stringsAsFactors = F)
pancan.filtered <- FilterMaf(pancan, snp.variants, "type")
unique.cancers <- unique(pancan.filtered$ttype)
output.folder <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/mutations/mutation_counts"
for (i in 1:length(unique.cancers)){
    temp <- table(pancan.filtered[pancan.filtered$ttype == unique.cancers[i], ]$patient)
    write.csv(temp, paste(output.folder, paste(unique.cancers[i], "csv", sep = "."), sep = "_"))
}


## calculate mean difference from expected for allelic fraction of NF2

targets <- master.table[(master.table$Analsysis.Set. == 1) & master.table$NF2.snp.indel == 1,"Pair.Name"]
targets <- targets[!is.na(targets)]
output <- matrix(NA, length(targets), 3)
rownames(output) <- targets
colnames(output) <- c("NF2_AF", "Avg_AF", "T_Stat")
for (i in 1:length(targets)){

    muts <- FilterMaf(total.snindels, targets[i], "Tumor_Sample_Barcode")
    nf2.af <- muts[muts$Hugo_Symbol == "NF2", "i_tumor_f"]
    expected.af <- nf2.af / 2
    muts <- muts[muts$Hugo_Symbol != "NF2", ]
    all.af <- muts$i_tumor_f
    for (j in 1:length(all.af)){
        if (all.af[j] > .7*nf2.af){
            all.af[j] <- all.af[j] / 2
        } 
    }
    output[i, ] <- c(nf2.af, mean(all.af), TestStatistic(all.af, expected.af))
}

output
