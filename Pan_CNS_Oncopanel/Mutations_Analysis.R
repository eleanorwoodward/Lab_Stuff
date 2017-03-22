## analysis of mutation data from Oncopanel calls

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

setwd("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/Analysis")



## read in cleaned master sheet with decoding
master.sheet <- read.delim("Master_Sheet_R_Upload.txt", stringsAsFactors = F)

table(master.sheet$Cancer_Type_Broad)

## generate disease lists
chondrosarcoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Chondrosarcoma", ]$SAMPLE_ACCESSION_NBR)
craniopharyngioma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Craniopharyngioma", ]$SAMPLE_ACCESSION_NBR)    
ependymoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Ependymoma", ]$SAMPLE_ACCESSION_NBR)
glioma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Glioma", ]$SAMPLE_ACCESSION_NBR)
lymphoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Lymphoma", ]$SAMPLE_ACCESSION_NBR)
medulloblastoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Medulloblastoma", ]$SAMPLE_ACCESSION_NBR)
meningioma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Meningioma", ]$SAMPLE_ACCESSION_NBR)
metastasis <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Metastasis", ]$SAMPLE_ACCESSION_NBR)
neuroblastoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Neuroblastoma", ]$SAMPLE_ACCESSION_NBR)
pineal_parenchymal_tumor <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Pineal_parenchymal_tumor", ]$SAMPLE_ACCESSION_NBR)
pituitary_adenoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Pituitary_adenoma", ]$SAMPLE_ACCESSION_NBR)
pituitary_other <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Pituitary_other", ]$SAMPLE_ACCESSION_NBR)
schwannoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Schwannoma", ]$SAMPLE_ACCESSION_NBR)


pathologies <- list(chondrosarcoma, craniopharyngioma, ependymoma, glioma, lymphoma, medulloblastoma, meningioma, metastasis, neuroblastoma, 
                    pineal_parenchymal_tumor, pituitary_adenoma, pituitary_other,schwannoma)

names(pathologies) <- c("chondrosarcoma", "craniopharyngioma", "ependymoma", "glioma", "lymphoma", "medulloblastoma", "meningioma", "metastasis", "neuroblastoma", 
                            "pineal_parenchymal_tumor", "pituitary_adenoma", "pituitary_other","schwannoma")



########### exploratory analysis plots

## Generates bar plots

pdf("all_tier_1_to_3_mutations.pdf", width = 14)
temp <- all.mutations.tier1.3
temp$BEST_EFF_GENE <- factor(temp$BEST_EFF_GENE, levels = names(sort(table(temp$BEST_EFF_GENE), decreasing = T)))
max <- 80
ggplot(data = subset(temp, BEST_EFF_GENE %in% levels(temp$BEST_EFF_GENE)[1:max]), aes(x = BEST_EFF_GENE)) + geom_bar() + 
labs(title = "All Tier 1 to 3 mutations") + theme(axis.text.x=element_text(angle=90, hjust=1)) + rameen_theme
dev.off()

## generates graph for specific subtypes of most frequently mutated genes
file.path.base <- "_all_tier_1_to_4_mutations.pdf"
title.base.count <- "All Tier 1 to 4 mutations in "
title.base.percent <- 

for (i in 1:length(pathologies)){
    file.name <- paste(names(pathologies)[i], file.path.base, sep = "")
    pdf(file.name, width = 14)
    title <- paste(title.base, names(pathologies)[i])
    temp <- all.mutations.tier1.4[all.mutations.tier1.4$SAMPLE_ACCESSION_NBR %in% pathologies[[i]], ]
    temp.deduped <- PerSampleMaf(temp, "BEST_EFF_GENE", "SAMPLE_ACCESSION_NBR")
    temp.deduped$BEST_EFF_GENE <- factor(temp.deduped$BEST_EFF_GENE, levels = names(sort(table(temp.deduped$BEST_EFF_GENE), decreasing = T)))
    print(ggplot(data = subset(temp.deduped, BEST_EFF_GENE %in% levels(temp$BEST_EFF_GENE)[1:80]), aes(x = BEST_EFF_GENE)) + geom_bar() + labs(title = title) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) + rameen_theme)
    dev.off()
}

## generates comuts

cutoffs <- c(1, 1, 1, .1, 1, 1, .5, .5, 1, 1, 1, 1, 1)
for (i in 1:length(pathologies)){
    
    PlotComut(maf = subset(all.mutations.tier1.3, SAMPLE_ACCESSION_NBR %in% pathologies[[i]]), 
              samples = subset(master.sheet, SAMPLE_ACCESSION_NBR %in% pathologies[[i]], select = SAMPLE_ACCESSION_NBR), input.samples = "SAMPLE_ACCESSION_NBR", 
              input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification", gene.cutoff = cutoffs[i], file.path = paste(names(pathologies)[i], "1-3 comut.pdf"), 
              unique.muts = unique(all.mutations$variant_classification), phenotypes = subset(master.sheet, SAMPLE_ACCESSION_NBR %in% pathologies[[i]], 
              select = c("SAMPLE_ACCESSION_NBR","PANEL_VERSION", "CNV_ONC")), title = paste(names(pathologies)[i], " comut"))
}



## recurrent hotspots across the cohort
## Hotspot finder

all.mutations.hotspot <- all.mutations

#excludes major categories
all.mutations.hotspot <- all.mutations[!(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["glioma"]],pathologies[["meningioma"]], pathologies[["metastasis"]])), ] 

## individual tumor types
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["glioma"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["meningioma"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["metastasis"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["glioma"]])), ]

gene.list <- sort(unique(all.mutations.hotspot$BEST_EFF_GENE))
mtrx <- matrix(0, 41, length(gene.list))
colnames(mtrx) <- gene.list

## loops through all genes present in recurrent maf
for(i in 1:length(gene.list)){
    temp <- FilterMaf(all.mutations.hotspot, gene.list[i], "BEST_EFF_GENE")
    x <- table(temp$POSITION)
    x <- sort(unname(x), decreasing = TRUE)
    
    ## if there are more than 40 start positions, takes only the first 40
    if (length(x) > 40){
        x <- x[1:40]
    ## otherwise, appends trailing zeros to produce 40 total entries for each    
    }else if (length(x) < 40){
        temp <- 40 - length(x)
        x <- c(x, rep(0, temp))
    }
    
    ## adds information back onto matrix
    mtrx[,i] <- c(x, 0)
}

## Sliding weighted average to determine relative importance
for(i in 1: ncol(mtrx)){
    if (mtrx[1,i] == 1){
        mtrx[41, i] <- 0
    }else if (sum(mtrx[, i]) > 100){
        mtrx[41, i] <- round(sum(mtrx[1:5, i]) / sum(mtrx[1:40, i]), 3)
    }else if (sum(mtrx[, i]) > 50){
        mtrx[41, i] <- round(sum(mtrx[1:4, i]) / sum(mtrx[1:40, i]), 3)
    }else if (sum(mtrx[, i]) > 20){
        mtrx[41, i] <- round(sum(mtrx[1:3, i]) / sum(mtrx[1:40, i]), 3)
    }else if (sum(mtrx[, i]) > 5){
        mtrx[41, i] <- round(sum(mtrx[1:2, i]) / sum(mtrx[1:40, i]), 3)
    }else {
        mtrx[41, i] <- 0
    }
}

## sort mutations based on weighted concentration of hotspots
df <- data.frame(mtrx)
df <- df[, order(df[41, ], df[1, ], decreasing = TRUE)]

pdf("hotspots_in_all.pdf", width = 25)
pdf("hotspots_in_all_min_10.pdf", width = 25)
pdf("hotspots_in_non_met_men_glio.pdf", width = 25)
pdf("hotspots_in_glioma.pdf", width = 25)
pdf("hotspots_in_meningioma.pdf", width = 25)
pdf("hotspots_in_metastasis.pdf", width = 25)

barplot(as.matrix(df[1:40, c(1:100)]), main = "Mutations grouped by base start position in all tumors, min 10 total", las = 2, col = c("gray8", "gray47", "gray70", "gray100"))
dev.off()



### look for tumor supressor genes


## Look at damaging mutations only, with cutoff from previous step
all.mutations.tsg <- all.mutations
table2 <- table(all.mutations.tsg$BEST_EFF_GENE)
table1 <- table(FilterMaf(all.mutations.tsg, c("Frameshift", "Inframe_Ins", "Nonsense_Mutation", "Frameshift_mutation", "Nonstop_Mutation", "In_frame_del", 
                                               "Stop_Lost", "In_Frame_Del", "Frame_Shift_Del", "In_Frame_Ins", "Nonsense", "Frame_Shift_Ins", "Inframe_Del",
                                               "Nonsense_mutation"),"CANONICAL_VARIANT_CLASS")$BEST_EFF_GENE)
table3 <- EqualizeTable(table1,table2)
table3 <-  rbind(table3, rep(0, ncol(table3)))

## Get ratio of damaging to total
for(i in 1:ncol(table3)){
    damaging <- table3[1, i]
    total <- table3[2, i]
    if (total < 5){
        table3[3, i] <- 0
    }else{
        table3[3, i] <- round(damaging / total, 2)
    }
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

pdf("TSG_in_all.pdf")
barplot(as.matrix(df[1:2, 1: (cutoff2 + 2)]), main = "Damaging Mutations vs Total Coding", las = 2, beside = TRUE, col = colors)
legend("topright", c("Greater than 80% damaging", "Greater than 40% Damaging", "Total Coding Mutations"), 
       col = c("red", "lightsalmon3", "grey"), pch = 15)
dev.off()
final.list <- names(df)[1:(cutoff1 - 1)]



## add preprocessing to consolidate mutation class
## turn comut code into function





## plot mutations per gene with color based on underlying pathology
temp <- merge(all.mutations.tier1.4, master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad")], "SAMPLE_ACCESSION_NBR")
temp$BEST_EFF_GENE <- factor(temp$BEST_EFF_GENE, levels = names(sort(table(temp$BEST_EFF_GENE), decreasing = T)))
print(ggplot(data = subset(temp, BEST_EFF_GENE %in% levels(temp$BEST_EFF_GENE)[1:80]), aes(x = BEST_EFF_GENE, fill = Cancer_Type_Broad)) + geom_bar() + labs(title = "title") + 
    rameen_theme)


## plot percent of samples with a given alteration, controlling for # in which gene is covered





