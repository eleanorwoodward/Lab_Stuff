## analysis of mutation data from Oncopanel calls

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

setwd("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/Analysis")

## read in cleaned master sheet with decoding
master.sheet <- read.delim("Master_Sheet_R_Upload.txt", stringsAsFactors = F)

table(master.sheet$Cancer_Type_Broad)

## generate disease lists
master.sheet <- master.sheet[master.sheet$PANEL_VERSION > 0, ]
master.sheet <- master.sheet[master.sheet$exclude != "no_tumor", ]

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
file.path.base <- " 1-3 mutation"
title.base <- "All Tier 1 to 3 mutations in"

for (i in 1:length(pathologies)){
    ## percent barplot
    file.name <- paste(names(pathologies)[i], file.path.base, " percent.pdf", sep = "")
    title <- paste(title.base, names(pathologies)[i], " percent")
    temp <- subset(all.mutations.tier1.3, SAMPLE_ACCESSION_NBR %in% pathologies[[i]])
    temp.deduped <- PerSampleMaf(temp, "BEST_EFF_GENE", "SAMPLE_ACCESSION_NBR")
    pdf(file.name, width = 14)
    PlotHistogramCoverage(maf = subset(temp.deduped, BEST_EFF_GENE %in% names(sort(table(temp.deduped$BEST_EFF_GENE), decreasing = T))[1:40]),
                                       gene.column = "BEST_EFF_GENE", samples = subset(master.sheet, SAMPLE_ACCESSION_NBR %in% pathologies[[i]], 
                                        select = SAMPLE_ACCESSION_NBR)$SAMPLE_ACCESSION_NBR, title)
    dev.off()
    
    ## count barplot
    file.name <- paste(names(pathologies)[i], file.path.base, " count.pdf", sep = "")
    title <- paste(title.base, names(pathologies)[i], " count")
    pdf(file.name, width = 14)
    temp.deduped$BEST_EFF_GENE <- factor(temp.deduped$BEST_EFF_GENE, levels = names(sort(table(temp.deduped$BEST_EFF_GENE), decreasing = T)))
    print(ggplot(data = subset(temp.deduped, BEST_EFF_GENE %in% levels(temp.deduped$BEST_EFF_GENE)[1:40]), aes(x = BEST_EFF_GENE)) + geom_bar() + labs(title = title) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) + rameen_theme)
    dev.off()
}


master.sheet$Cancer_Type_Details[master.sheet$Cancer_Type_Details == "-"] <- 2




## generates comuts with relevant data
for (i in 1:length(pathologies)){
    ## determine threshold for inclusion
    if (length(pathologies[[i]]) < 30){
        min.val <- 2
    }else{
        min.val <- length(pathologies[[i]]) / 20
    }

    ## only include tier 4 variants if same gene is labeled as tier 1-3 variant in cancer type
    temp.snv <- subset(all.mutations.tier1.4, SAMPLE_ACCESSION_NBR %in% pathologies[[i]])
    temper.snv <- subset(all.mutations.tier1.3, SAMPLE_ACCESSION_NBR %in% pathologies[[i]])
    temp.snv.1 <- temp.snv[temp.snv$BEST_EFF_GENE %in% temper.snv$BEST_EFF_GENE, ]
    
    genes <- sort(table(temp.snv.1$BEST_EFF_GENE), decreasing = T)

    alterations <- temp.snv.1[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
    alterations.cnv <- subset(all.cnv.high, SAMPLE_ACCESSION_NBR %in% pathologies[[i]], 
                       select = c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD"))
    alterations.sv <- all.svs.formatted[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
    
    colnames(alterations.cnv) <- colnames(alterations)
    colnames(alterations.sv) <- colnames(alterations)
    alterations <- rbind(alterations, alterations.cnv, alterations.sv)
    genes <- sort(table(alterations$BEST_EFF_GENE))
    alterations <- alterations[alterations$BEST_EFF_GENE %in% names(genes)[genes > min.val], ]
    alterations$variant_classification[!alterations$variant_classification %in% c("2DEL", "HA", "missense", "rearrangement")] <- "damaging mutation"
    
    temp.cnv.broad <- all.cnvs.broad[, colnames(all.cnvs.broad) %in% pathologies[[i]]]
    idx <- apply(temp.cnv.broad, 1, function(x){sum(x != 0)})
    idx <- as.numeric(idx) > min.val
    temp.cnv.broad <- temp.cnv.broad[idx, ]
    if (nrow(temp.cnv.broad) == 0){
        temp.cnv.broad <- NA
    }
    
    PlotComut(
              mut.maf1 = alterations, 
              broad.cnv.maf = temp.cnv.broad,
              samples = master.sheet[master.sheet$SAMPLE_ACCESSION_NBR %in% pathologies[[i]], "SAMPLE_ACCESSION_NBR"], 
              input.samples = "SAMPLE_ACCESSION_NBR", 
              input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification", 
              file.path = paste("Comut ", names(pathologies)[i], ".pdf", sep = ""),
              #file.path = paste("Comut by age.pdf"),
              #y.axis.font = 15, 
              dimensions = c(18,10), 
              col.vector = c("frameshift_indel" = "red", "missense" ="skyblue" , "nonsense" = "orange", "splice_site" = "yellow", "stop_codon" = "purple", "in_frame_indel" = "gold",
                             "other" = "cyan", "arm-level gain" = "red", "arm-level loss" = "blue", "2DEL" = "blue", "HA" = "red", "damaging mutation" = "orange", "rearrangement" = "green",
                             "Glioblastoma" = "black", "AnaplasticAstro" = "green", "Angiocentric" = "orange", "Astro" = "green", 
                             "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                             "0" = "gray90", "1" = "gray70", "2" = "gray50", "3" = "gray25", "4" = "gray1"),
              #fixed.order = clust$labels[clust$order],
              #manual.order = c("7p", "7q", "10q", "MDM2-cnv", "MDM2", "MDM4-cnv", "MDM4"),
              #manual.order = c("1p", "19q", "7p", "10q", "TP53", "IDH1"), 
              #manual.order = c("Primary", "7p"), 
              phenotypes = subset(master.sheet, SAMPLE_ACCESSION_NBR %in% pathologies[[i]], 
              select = c("SAMPLE_ACCESSION_NBR","Cancer_Type_Specific")), title = paste(names(pathologies)[i], " comut"))

}


## Glioma comut

## traiditonal criteria of tier 1-3 with added 4 for x most common genes
temp.snv <- subset(all.mutations.tier1.4, SAMPLE_ACCESSION_NBR %in% pathologies[["glioma"]])
temper.snv <- subset(all.mutations.tier1.3, SAMPLE_ACCESSION_NBR %in% pathologies[["glioma"]])
temp.snv.1 <- temp.snv[temp.snv$BEST_EFF_GENE %in% temper.snv$BEST_EFF_GENE, ]
temp.snv.1 <- temp.snv.1[temp.snv.1$BEST_EFF_GENE %in% names(sort(table(temp.snv.1$BEST_EFF_GENE), decreasing = T))[1:20], ]


temp.cnv.focal <- subset(all.cnv.high, SAMPLE_ACCESSION_NBR %in% pathologies[["glioma"]], 
                         select = c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD"))

temp.cnv.focal <- temp.cnv.focal[temp.cnv.focal$GENE %in% names(sort(table(temp.cnv.focal$GENE), decreasing = T))[1:7], ]

temp.cnv.broad <- all.cnvs.broad[, colnames(all.cnvs.broad) %in% pathologies[["glioma"]]]

## Rare comut
## traiditonal criteria of tier 1-3 with added 4 for x most common genes
temp.snv <- all.mutations.tier1.4[!(all.mutations.tier1.4$SAMPLE_ACCESSION_NBR %in% unlist(pathologies[c("glioma", "meningioma", "pituitary", "metastasis")])), ]
temper.snv <- subset(all.mutations.tier1.3, SAMPLE_ACCESSION_NBR %in% temp.snv$SAMPLE_ACCESSION_NBR)
temp.snv.1 <- temp.snv[temp.snv$BEST_EFF_GENE %in% temper.snv$BEST_EFF_GENE, ]

temp.cnv.focal <- all.cnv.high[!(all.cnv.high$SAMPLE_ACCESSION_NBR %in% unlist(pathologies[c("glioma", "meningioma", "pituitary", "metastasis")])), ]

temp.cnv.broad <- all.cnvs.broad[, !(colnames(all.cnvs.broad) %in% unlist(pathologies[c("glioma", "meningioma", "pituitary", "metastasis")])), ]

## combine all mutation classes together
alterations <- temp.snv.1[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
alterations.cnv <- all.cnv.high[, c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD")]
colnames(alterations.cnv) <- colnames(alterations)
alterations.svs <- all.svs.formatted[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
alterations <- rbind(alterations, alterations.cnv, alterations.svs)
alterations$BEST_EFF_GENE[alterations$BEST_EFF_GENE %in% c("CDKN2A", "CDKN2B")] <- "CDKN2A/B"

alterations$variant_classification_simple <- "damaging mutation"
alterations$variant_classification_simple[alterations$variant_classification == "2DEL"] <- "focal loss"
alterations$variant_classification_simple[alterations$variant_classification == "HA"] <- "focal gain"
alterations$variant_classification_simple[alterations$variant_classification == "rearrangement"] <- "rearrangement"
alterations$variant_classification_simple[alterations$variant_classification == "missense"] <- "missense"
alterations <- alterations[alterations$BEST_EFF_GENE %in% names(sort(table(alterations$BEST_EFF_GENE), decreasing = TRUE)[1:30]), ]
clinical <- master.sheet[!(master.sheet$SAMPLE_ACCESSION_NBR %in% unlist(pathologies[c("glioma", "meningioma", "pituitary", "metastasis")])),
                         c("SAMPLE_ACCESSION_NBR","Cancer_Type_Broad", "Primary", "Grade")]

PlotComut(
    mut.maf1 = alterations, 
    #mut.maf2 = temp.snv.2, 
    #focal.cnv.maf = temp.cnv.focal,
    broad.cnv.maf = temp.cnv.broad,
    samples = master.sheet[master.sheet$SAMPLE_ACCESSION_NBR %in% clinical$SAMPLE_ACCESSION_NBR, "SAMPLE_ACCESSION_NBR"],
    input.samples = "SAMPLE_ACCESSION_NBR", 
    input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification_simple",
    dimensions = c(18,10), 
    col.vector = c("damaging mutation" = "orange", "focal loss" = "blue", "focal gain" = "red", "rearrangement" = "green", "missense" = "skyblue",
                   "arm-level gain" = "red", "arm-level loss" = "blue","Glioblastoma" = "black", "AnaplasticAstro" = "green", "Pediatric" = "purple",
                   "Angiocentric" = "orange", "Astro" = "green", "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                   "0" = "gray90", "1" = "gray70", "2" = "gray50", "3" = "gray25", "4" = "gray1"),
    #fixed.order = clust$labels[clust$order],
    manual.order = c("Cancer_Type_Broad"),
    #file.path = "../Analysis/GBM Comut All Alterations by Cancer Type TP53 EGFR.pdf",
    phenotypes = clinical)



## recurrent hotspots across the cohort
## Hotspot finder

all.mutations.hotspot <- all.mutations.tier1.4

#excludes major categories
all.mutations.hotspot <- all.mutations[!(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["glioma"]],pathologies[["meningioma"]], pathologies[["metastasis"]])), ] 

## individual tumor types
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["pituitary_adenoma"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["meningioma"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$SAMPLE_ACCESSION_NBR %in% c(pathologies[["metastasis"]])), ]

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

barplot(as.matrix(df[1:40, c(1:10)]), las = 2, col = c("gray8", "gray47", "gray70", "gray100"))
dev.off()

write.csv(all.mutations[all.mutations$SAMPLE_ACCESSION_NBR %in% pathologies[[11]], ], "C:/Users/Noah/OneDrive/Work/Presentations/Linda's stuff/WFNS/adenoma.mutations.csv")

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



## plot mutations per gene with color based on underlying pathology
temp <- merge(all.mutations.tier1.4, master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad")], "SAMPLE_ACCESSION_NBR")
temp$BEST_EFF_GENE <- factor(temp$BEST_EFF_GENE, levels = names(sort(table(temp$BEST_EFF_GENE), decreasing = T)))
pdf("../Analysis/barplot mutations by tumor type.pdf", width = 15)
print(ggplot(data = subset(temp, BEST_EFF_GENE %in% levels(temp$BEST_EFF_GENE)[1:80]), aes(x = BEST_EFF_GENE, fill = Cancer_Type_Broad)) + geom_bar() + labs(title = "title") + 
    rameen_theme)
dev.off()
## chi-squared test to determine which genes are differentially mutated

top.genes <- names(sort(table(all.mutations.tier1.4$BEST_EFF_GENE), decreasing = T))





## plot 2D grid with size corresponding to incidence
temp <- data.frame(c("Class1", "Class2", "Class3", "Class1", "Class2", "Class3"), 
                   c("Gene1", "Gene1", "Gene1", "Gene2", "Gene2", "Gene2"), 
                   c(1,2,3,4,5,6))
colnames(temp) <- c("Class", "Gene", "Freq")

ggplot(data = temp, aes(x = Class, y = Gene)) + geom_point(aes(size = Freq))

## use real data
## generate gene list of top hits

temp <- merge(all.mutations.tier1.3, master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad")])
genes <- names(sort(table(temp$BEST_EFF_GENE), decreasing = TRUE)[1:20])
genes <- c(genes, "CTNNB1", "MYD88", "CD79B", "PIM1", "PTCH1", "KMT2D", "AKT1", "SMO", "PTPRD", "PRKDC", "DICER1", "GNAS")
temp <- temp[temp$BEST_EFF_GENE %in% genes, ]
temp <- temp[temp$Cancer_Type_Broad != "", ]
temp2 <- table(temp$Cancer_Type_Broad, temp$BEST_EFF_GENE)

## option to convert table to percentage
vals <- as.numeric(table(master.sheet[!(master.sheet$Cancer_Type_Broad %in% c("", "Schwannoma check", "Pituitary_other")), ]$Cancer_Type_Broad))
temp2 <- apply(temp2, 2, function(x){x / vals})

melted <- melt(temp2)
colnames(melted) <- c("Tumor", "Gene", "Percent")

## re order levels based on number of genes mutated per tumor type
orders <- sort(apply(temp2,1, function(x){sum(x!=0)}), decreasing = T)
melted$Tumor <- factor(melted$Tumor, levels = names(orders))
orders <- sort(apply(temp2, 2, function(x){sum(x!=0)}))
melted$Gene <- factor(melted$Gene, levels = names(orders))
melted <- melted[!(melted$Count == 0), ]
pdf("../Analysis/2D Circle Plot Mutations.pdf")
ggplot(data = melted, aes(x = Tumor, y = Gene)) + geom_point(aes(size = Percent)) + rameen_theme
dev.off()


## look at permutations
temp <- df
temp.glioma <- df[df$samples %in% df[df$mutations == "Glioma", ]$samples,] 
temp.glioma <- temp.glioma[temp.glioma$mutations != 0, ]


## hypermutators analysis
hyper <- master.sheet$SAMPLE_ACCESSION_NBR[as.numeric(master.sheet$SNV_COUNT) > 50]
hyper <- hyper[!is.na(hyper)]
hyper.info <- master.sheet[master.sheet$SAMPLE_ACCESSION_NBR %in% hyper, ]
hyper.info <- hyper.info[, c("MRN", "Cancer_Type_Broad", "Cancer_Type_Specific", "SNV_COUNT", "CHILDRENS_MRN", "SAMPLE_ACCESSION_NBR")]
hyper.mutations <- all.mutations.tier1.4[all.mutations.tier1.4$BEST_EFF_GENE %in% c("MSH2", "MSH6", "MLH1"), ]
hyper.info.combined <- merge(hyper.info, hyper.mutations[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "BEST_EFF_PROTEIN_CHANGE", "BEST_EFF_CDNA_CHANGE", "variant_classification")])
write.csv(hyper.info.combined, "Hypermutator_sample_info.csv", quote = FALSE, row.names = FALSE)
