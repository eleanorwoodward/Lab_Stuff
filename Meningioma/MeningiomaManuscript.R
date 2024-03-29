## Meningioma Paper Calcs Generator
## Replicates all calculations made in paper for easy repeatability/updatedability

## Load excel table with categories to be correlated
big.daddy <- read.delim("C:/Users/Noah/OneDrive/Meningioma/Analysis...", stringsAsFactors = FALSE)
mutations.per.sample <- table(snindels$Tumor_Sample_Barcode)
sort(names(mutations.per.sample))
strsplit("M103-tumor", "-")

nf2 <- c(0,1,1,1,0,0,0)
muts <- c(10, 100, 20, 40, 4, 10, 23)

mtx <- cbind(nf2, muts)
colnames(mtx) <- c("NF2 status", "Number of mutations")
cor(mtx)
cor.test(mtx)


## Screen for real mutations at given cutoff
snindels.elim <- ReccurentMaf(snindels, "Hugo_Symbol", 4)

snindels.tot <- rbind(snindels, snindels.disc)
snindels.tot <- ReccurentMaf(snindels.tot, "Hugo_Symbol", 4)
snindels.elim <- snindels.tot
snindels.nf <- rbind(snindels.nf, snindels.disc.silent)


coding.variants <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Nonstop_Mutation", 
                     "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")
noncoding.variants <- "Silent"

silent.maf <- FilterMaf(snindels.nf, c(coding.variants, noncoding.variants), "Variant_Classification")

PlotMaf(snindels.elim, "Hugo_Symbol", 2, "Mutations in Genes with at least 5 hits")

small <- FilterMaf(silent.maf, snindels.elim$Hugo_Symbol, "Hugo_Symbol")
silent <- FilterMaf(small, noncoding.variants, "Variant_Classification")
loud <- FilterMaf(small, coding.variants, "Variant_Classification")
silent.muts <- table(silent$Hugo_Symbol)
loud.muts <- table(loud$Hugo_Symbol)
table3 <- EqualizeTable(loud.muts, silent.muts)
table3 <- rbind(table3, rep(1, ncol(table3)))
rownames(table3) <- c("Loud", "Silent", "Ratio")


## Find ratio of silent to nonsilent

for(i in 1:length(colnames(table3))){
  silent <- table3[2, i]
  coding <- table3[1, i]
  if (silent == 0){
    table3[3, i] <- 0
  }else{
    table3[3, i] <- round(silent / coding, 2)
  }
}
df <- data.frame(table3)
df <- df[, order(df[3, ], -df[1, ])]

## Take only those with given ratio
idx1 <- df[3, ] > .2
idx2 <- df[3, ] > .3
idx3 <- df[3, ] > .4
idx4 <- df[3, ] > .5
cutoff1 <- which(idx1)[1]
cutoff2 <- which(idx2)[1]
cutoff3 <- which(idx3)[1]
cutoff4 <- which(idx4)[1]
cols <- rep(c("red", "grey"), cutoff1 - 1)
cols <- c(cols, rep(c("lightsalmon3", "grey"), cutoff2 -  cutoff1))
cols <- c(cols, rep(c("lightblue", "grey"), cutoff3 - cutoff2))
colors <- c(cols, rep("grey", (ncol(df) - cutoff3 + 1) * 2 ))
barplot(as.matrix(df[1:2, 1: cutoff4]), beside = TRUE, main = "Silent Mutations vs Coding Mutations", 
        las = 2, col = colors)
legend("topright", c("Fewer than 20% silent", "Fewer than 30% silent", "Fewer than 40% Silent", "Silent Mutations"), 
       col = c("red", "lightsalmon3", "lightblue", "grey"), pch = 15)


## Look at damaging mutations only, with cutoff from previous step
snindels.elim.silent <- FilterMaf(snindels.elim, names(df)[1:cutoff1], "Hugo_Symbol")
table2 <- table(snindels.elim.silent$Hugo_Symbol)
table1 <- table(FilterMaf(snindels.elim.silent, c("Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "Frame_Shift_Del", "De_novo_Start_OutOfFrame"), 
                          "Variant_Classification")$Hugo_Symbol)
table3 <- EqualizeTable(table1,table2)
table3 <-  rbind(table3, rep(0, ncol(table3)))

## Get ratio of damaging to total
for(i in 1:ncol(table3)){
  damaging <- table3[1, i]
  total <- table3[2, i]
  table3[3, i] <- round(damaging / total, 2)
}

df <- data.frame(table3)
df <- df[, order(df[3, ], df[1, ], decreasing = TRUE)]
idx1 <- df[3, ] < .8
idx2 <- df[3, ] < .4
cutoff1 <- which(idx1)[1]
cutoff2 <- which(idx2)[1]
cols <- rep(c("red", "grey"), cutoff1 - 1)
cols <- c(cols, rep(c("lightsalmon3", "grey"), cutoff2 -  cutoff1))
colors <- c(cols, rep("grey", (ncol(df) - cutoff2 + 1) * 2 ))
barplot(as.matrix(df[1:2, 1: (cutoff2 + 2)]), main = "Damaging Mutations vs Total Coding", las = 2, beside = TRUE, col = colors)
legend("topright", c("Greater than 80% damaging", "Greater than 40% Damaging", "Total Coding Mutations"), 
       col = c("red", "lightsalmon3", "grey"), pch = 15)
final.list <- names(df)[1:(cutoff1 - 1)]

## Look for hotspot mutations
gene.list <- sort(unique(snindels.elim.silent$Hugo_Symbol))

## Construct matrix
mtrx <- matrix(0, 41, length(gene.list))
colnames(mtrx) <- gene.list
for(i in 1:length(gene.list)){
  temp <- FilterMaf(snindels.elim.silent, gene.list[i], "Hugo_Symbol")
  x <- table(temp$Start_position)
  x <- sort(unname(x), decreasing = TRUE)
  if (length(x) >= 40){
    x <- x[1:40]
  }
  
  if (length(x) < 40){
    temp <- 40 - length(x)
    x <- c(x, rep(0, temp))
  }
  mtrx[,i] <- c(x, 0)
}

## Percent in first 2 spots
for(i in 1: ncol(mtrx)){
  mtrx[41, i] <- round(sum(mtrx[1:2, i]) / sum(mtrx[1:40, i]), 3)
}

df <- data.frame(mtrx)
df <- df[, order(df[41, ], df[1, ], decreasing = TRUE)]
barplot(as.matrix(df[1:40, ]), main = "Mutations grouped by base start position", las = 2, col = c("gray8", "gray47", "gray70", "gray100"))
# 12 for validation, 95 for combined
abline(v = 95)
## corresponds to 10 for validation, 79 for combined
final.list <- c(final.list, names(df)[1:79])
final.list <- unique(final.list)
## Fig 3. Average allelic fraction of genes mutated at least twice
avgs <- AverageMaf(snindels.elim.silent, final.list, "Hugo_Symbol", "tumor_f")
names(avgs) <- final.list
avgs <- sort(avgs, decreasing = TRUE)
barplot(avgs, main = "Allelic Fraction of Screened Genes", las = 2)


