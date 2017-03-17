## analysis of mutation data from Oncopanel calls

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

setwd("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/Analysis")


## read in raw data 
all.mutations <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_MUTATION_RESULTS.csv", stringsAsFactors = F)
all.mutations.tier1.3 <- all.mutations[all.mutations$TIER_ID < 4, ]
all.mutations.tier1.4 <- all.mutations[all.mutations$TIER_ID < 5, ]
all.cnvs <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_CNV_RESULTS.csv", stringsAsFactors = F)
all.svs <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_SV_RESULTS.csv", stringsAsFactors = F)

## read in cleaned master sheet with decoding
master.sheet <- read.delim("Master_Sheet_R_Upload.txt", stringsAsFactors = F)

table(master.sheet$Cancer_Type_Broad)

## generate disease lists
chondrosarcoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Chondrosarcoma", ]$PATIENT_ID)
craniopharyngioma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Craniopharyngioma", ]$PATIENT_ID)    
ependymoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Ependymoma", ]$PATIENT_ID)
glioma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Glioma", ]$PATIENT_ID)
lymphoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Lymphoma", ]$PATIENT_ID)
medulloblastoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Medulloblastoma", ]$PATIENT_ID)
meningioma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Meningioma", ]$PATIENT_ID)
metastasis <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Metastasis", ]$PATIENT_ID)
neuroblastoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Neuroblastoma", ]$PATIENT_ID)
pineal_parenchymal_tumor <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Pineal_parenchymal_tumor", ]$PATIENT_ID)
pituitary_adenoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Pituitary_adenoma", ]$PATIENT_ID)
pituitary_other <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Pituitary_other", ]$PATIENT_ID)
schwannoma <- unique(master.sheet[master.sheet$Cancer_Type_Broad == "Schwannoma", ]$PATIENT_ID)


pathologies <- list(chondrosarcoma, craniopharyngioma, ependymoma, glioma, lymphoma, medulloblastoma, meningioma, metastasis, neuroblastoma, 
                    pineal_parenchymal_tumor, pituitary_adenoma, pituitary_other,schwannoma)

names(pathologies) <- c("chondrosarcoma", "craniopharyngioma", "ependymoma", "glioma", "lymphoma", "medulloblastoma", "meningioma", "metastasis", "neuroblastoma", 
                            "pineal_parenchymal_tumor", "pituitary_adenoma", "pituitary_other","schwannoma")



########### exploratory analysis plots

## generates graph for all samples
pdf("all_tier_1_to_3_mutations.pdf", width = 14)
temp <- all.mutations.tier1.3
temp$BEST_EFF_GENE <- factor(temp$BEST_EFF_GENE, levels = names(sort(table(temp$BEST_EFF_GENE), decreasing = T)))
max <- 80
ggplot(data = subset(temp, BEST_EFF_GENE %in% levels(temp$BEST_EFF_GENE)[1:max]), aes(x = BEST_EFF_GENE)) + geom_bar() + 
labs(title = "All Tier 1 to 3 mutations") + theme(axis.text.x=element_text(angle=90, hjust=1)) + rameen_theme
dev.off()

## generates graph for specific subtypes of most frequently mutated genes
file.path.base <- "all_tier_1_to_3_mutations_"
title.base <- "All Tier 1 to 3 mutations in"

for (i in 1:length(pathologies)){
    file.name <- paste(file.path.base, names(pathologies)[i], ".pdf", sep = "")
    pdf(file.name, width = 14)
    title <- paste(title.base, names(pathologies)[i])
    temp <- all.mutations.tier1.3[all.mutations.tier1.3$PATIENT_ID %in% pathologies[[i]], ]
    temp$BEST_EFF_GENE <- factor(temp$BEST_EFF_GENE, levels = names(sort(table(temp$BEST_EFF_GENE), decreasing = T)))
    print(ggplot(data = subset(temp, BEST_EFF_GENE %in% levels(temp$BEST_EFF_GENE)[1:80]), aes(x = BEST_EFF_GENE)) + geom_bar() + labs(title = title) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) + rameen_theme)
    dev.off()
}


## recurrent hotspots across the cohort
## Hotspot finder

all.mutations.hotspot <- all.mutations

#excludes major categories
all.mutations.hotspot <- all.mutations[!(all.mutations$PATIENT_ID %in% c(pathologies[["glioma"]],pathologies[["meningioma"]], pathologies[["metastasis"]])), ] 

## individual tumor types
all.mutations.hotspot <- all.mutations[(all.mutations$PATIENT_ID %in% c(pathologies[["glioma"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$PATIENT_ID %in% c(pathologies[["meningioma"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$PATIENT_ID %in% c(pathologies[["metastasis"]])), ]
all.mutations.hotspot <- all.mutations[(all.mutations$PATIENT_ID %in% c(pathologies[["glioma"]])), ]

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

## create comut
## adapted from code at https://benchtobioinformatics.wordpress.com/2015/05/25/how-to-make-a-co-mutation-plot/
samples <- unique(all.mutations$PATIENT_ID)[1:50]
genes <- unique(all.mutations$BEST_EFF_GENE)[1:10]
df <- expand.grid(genes, samples, stringsAsFactors = F)
colnames(df) <- c("genes", "sample")
df$mutation <- NA
## loops through each row in df, determining if given gene is mutated in given sample, and if so, what type
temp <- all.mutations.tier1.3

for (i in 1:nrow(df)){
    tmp.sample <- df$sample[i]
    tmp.gene <- df$genes[i]
    matches <- temp[temp$PATIENT_ID == tmp.sample & temp$BEST_EFF_GENE == tmp.gene, ]
    if (nrow(matches) == 0){
        ## do nothing
    }else{
        df$mutation[i] <- temp$BEST_EFF_VARIANT_CLASS[i]
    }
}


df_sub <- subset(df, !is.na(df$mutation))

ord <- names(sort(table(df_sub$genes), decreasing = T))

# re-order rows with most frequently mutated gene on top?re
df$genes <- factor(df$genes, levels = ord)
df <- df[order(df$genes), ]
df <- df[!is.na(df$genes), ]

## reshape to end result comut in order to properly order columns
df.wide <- reshape(df, v.names = "mutation", idvar = "sample", timevar = "genes", direction = "wide")
for (i in 2:ncol(df.wide)){
    df.wide[, i][is.na(df.wide[, i])] <- "z"
}
df.wide$sample %in% df$sample

## uses all columns
df.wide <- df.wide[do.call(order, df.wide[, -1]), ]

## takes given order, and sorts samples based on ordering
missing <- df$sample[!(df$sample %in% df.wide$sample)]
df$sample <- factor(df$sample, levels = c(df.wide$sample, missing))
str(df)

## switches factor order back
df$genes <- factor(df$genes, rev(levels(df$genes)))

# now for a Comut plot with ggplot2
mut <- ggplot(df, aes(x=sample, y=genes, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=mutation)) +
    scale_fill_brewer(palette = "Set1", na.value="Grey90") +
    xlab("Subject") +
    ggtitle("Example Comut plot") +
    theme(
        legend.key = element_rect(fill='NA'),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size=8, face="bold"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.text.y=element_text(colour="Black"),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.background=element_blank()
    ))



