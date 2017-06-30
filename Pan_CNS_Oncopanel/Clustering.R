## Clustering analysis
## use NMF to cluster
## adapted from: https://cran.r-project.org/web/packages/NMF/vignettes/heatmaps.pdf 
## https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
## https://cran.r-project.org/web/packages/NMF/NMF.pdf
## generate fake data for clustering
num.samples <- 200
num.genes <- 20
set.seed(50)
group1 <- t(sapply(sample(seq(0.01, .8, 0.01), num.genes, TRUE), function(x){rbinom(num.samples, 1, x)}))
group2 <- t(sapply(sample(seq(0.01, .8, 0.01), num.genes, TRUE), function(x){rbinom(num.samples, 1, x)}))
group3 <- t(sapply(sample(seq(0.01, .8, 0.01), num.genes, TRUE), function(x){rbinom(num.samples, 1, x)}))
group4 <- t(sapply(sample(seq(0.01, .8, 0.01), num.genes, TRUE), function(x){rbinom(num.samples, 1, x)}))
group5 <- t(sapply(sample(seq(0.01, .8, 0.01), num.genes, TRUE), function(x){rbinom(num.samples, 1, x)}))

input.group <- cbind(group1, group2, group3, group4, group5)

##possibility: nmfEstimateRank

test.2 <- nmf(testing, 2)

testing.r.3 <- nmf(input.group, 2:7, nrun = 20, seed = 123456)
plot(hello)
plot(testing.r.5)

saveRDS(testing.r.3, "../Upload to xchip/testing.object.RDS")
hello <- readRDS("../Upload to xchip/NMF_OUPUT.RDS")

coefmap(minfit(testing.r.3), annCol = colnames(input.group), annColors = list(c("red", "yellow", "blue", "green4", "purple", "orange")), Colv = "basis")


consensusmap(testing.r.5)

colnames(input.group) <- c(rep("A", num.samples), rep("B",num.samples), rep("C",num.samples), rep("D",num.genes), rep("D",num.genes))
write.table(input.group, "../Upload to xchip/example_nmf_input.txt", sep = "\t")



## first generate combined alterations list
alterations <- all.mutations.tier1.4[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
alterations.cnv <- all.cnv.high[, c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD")]
colnames(alterations.cnv) <- colnames(alterations)
alterations.svs <- all.svs.formatted[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
alterations <- rbind(alterations, alterations.cnv, alterations.svs)
alterations$BEST_EFF_GENE[alterations$BEST_EFF_GENE %in% c("CDKN2A", "CDKN2B")] <- "CDKN2A/B"
alterations$BEST_EFF_GENE[alterations$BEST_EFF_GENE == "MYCN"] <- "MYC"
alterations$variant_classification <- "altered"

## next restrict to only those present in at least 20% of tumor types that have at least 100 samples each
trimmed <- c()
for (i in 1:length(pathologies)){
    if (length(pathologies[[i]]) > 100){
        current.alterations <- alterations[alterations$SAMPLE_ACCESSION_NBR %in% pathologies[[i]], ]
        candidates <- sort(table(current.alterations$BEST_EFF_GENE), decreasing = TRUE)
        finalists <- names(candidates[candidates > floor(length(pathologies[[i]]) / 5)])
        trimmed <- c(trimmed, finalists)
    }
}

trimmed <- unique(trimmed)

## figure out which samples have no mutations
alterations.trimmed <- alterations[alterations$BEST_EFF_GENE %in% trimmed, ]
names <- unique(master.sheet$SAMPLE_ACCESSION_NBR)
quiet.names <- names[!(names %in% alterations.trimmed$SAMPLE_ACCESSION_NBR)]
quiet.samples <- cbind(quiet.names, "quiet.mut", "altered")
colnames(quiet.samples) <- colnames(alterations.trimmed)
alterations.trimmed <- rbind(alterations.trimmed, quiet.samples)

## Add meta features
clustering.cnvs.gain <- all.cnvs.broad
clustering.cnvs.gain[clustering.cnvs.gain < 0] <- 0
clustering.cnvs.loss <- all.cnvs.broad
clustering.cnvs.loss[clustering.cnvs.loss > 0] <- 0
clustering.cnvs.loss <- abs(clustering.cnvs.loss)
rownames(clustering.cnvs.gain) <- paste(unique(gene.to.arm$ID), "_gain", sep = "")
rownames(clustering.cnvs.loss) <- paste(unique(gene.to.arm$ID), "_loss", sep = "")

clustering.cnvs <- rbind(clustering.cnvs.gain, clustering.cnvs.loss)
clustering.cnvs["count", ] <- as.numeric(apply(clustering.cnvs, 2, function(x){sum(x != 0)}))
quantile(as.numeric(clustering.cnvs[69, ]))
clustering.cnvs["cnv.disrupted", clustering.cnvs["count", ] > 6 ] <- 1
clustering.cnvs["cnv.disrupted", is.na(clustering.cnvs["cnv.disrupted", ])] <- 0
clustering.cnvs["cnv.quiet", clustering.cnvs["count", ] == 0] <- 1
clustering.cnvs["cnv.quiet", is.na(clustering.cnvs["cnv.quiet", ])] <- 0
clustering.cnvs["cnv.middle", clustering.cnvs["count", ] %in% 1:5] <- 1
clustering.cnvs["cnv.middle", is.na(clustering.cnvs["cnv.middle", ])] <- 0
clustering.cnvs <- clustering.cnvs[-69, ]

## same as above, figure out which features are present in at least 20% of cases
trimmed.broad <- c()
for (i in 1:length(pathologies)){
    if (length(pathologies[[i]]) > 100){
        current.samples <- clustering.cnvs[, names.fixed %in% pathologies[[i]]]
        candidates <- rowSums(current.samples)
        finalists <- rownames(current.samples)[candidates > floor(length(pathologies[[i]]) / 5)]
        trimmed.broad <- c(trimmed.broad, finalists)
    }
}

trimmed.broad <- unique(trimmed.broad)

clustering.cnvs.trimmed <- clustering.cnvs[trimmed.broad, ]


clustering.clinical <- master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad", "Cancer_Type_Specific")]

clustering.input.cnvs <- PlotComut(
    mut.maf1 = alterations.trimmed,
    broad.cnv.maf = clustering.cnvs.trimmed,
    samples = master.sheet$SAMPLE_ACCESSION_NBR,
    input.samples = "SAMPLE_ACCESSION_NBR", 
    input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification", 
    return.matrix = TRUE,
    dimensions = c(18,10),
    #manual.order = c("quiet.mut", "cnv.quiet", "cnv.disrupted"), 
    col.vector = c("frameshift_indel" = "red", "missense" ="red" , "nonsense" = "red", "splice_site" = "red", "stop_codon" = "red", "in_frame_indel" = "red",
                   "other" = "red", "arm-level gain" = "red", "arm-level loss" = "red", "2DEL" = "red", "HA" = "red", "indel" = "red", "rearrangement" = "red",
                   "Glioblastoma" = "black", "AnaplasticAstro" = "green", "Angiocentric" = "orange", "Astro" = "green", "altered" = "red", 
                   "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                   "0" = "gray90", "1" = "gray70", "2" = "gray50", "3" = "gray25", "4" = "gray1")
    ,phenotypes = clustering.clinical
)


## change to non-negative matrix, features as rows and samples as columns
clustering.input.gliomas <- clustering.input[clustering.input$samples %in% master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Broad == "Glioma"], ]
nmf.input <- reshape(clustering.input, idvar = "samples", timevar = "genes", direction = "wide")

nmf.input[nmf.input == "wt"] <- 0
nmf.input[nmf.input == "nc"] <- 0
nmf.input[, -1][nmf.input[, -1] != 0] <- 1
nmf.input <- t(nmf.input)
colnames(nmf.input) <- nmf.input[1, ]
nmf.input <- nmf.input[-1, ]
row.labels <- rownames(nmf.input)
nmf.input <- apply(nmf.input, 2, as.numeric)
rownames(nmf.input) <- gsub("(mutations.)(*)", "\\2", row.labels)

## remove empty columns
## figure out number of covariates incldued in clinical clustering, remove that many from end
end <- nrow(nmf.input)
remove <- ncol(clustering.clinical) - 2
start <- end - remove
nmf.input <- nmf.input[-c(start:end), ]
nmf.input <- nmf.input[, colSums(nmf.input) != 0]


column.annotations <- master.sheet
column.annotations <- column.annotations[column.annotations$SAMPLE_ACCESSION_NBR %in% colnames(nmf.input), ]
column.annotations$Cancer_Type_Broaddd <- column.annotations$Cancer_Type_Broad
column.annotations$Cancer_Type_Glioma[column.annotations$Cancer_Type_Broad == "Glioma"] <- 
    column.annotations$Cancer_Type_Specific[column.annotations$Cancer_Type_Broad == "Glioma"]
column.annotations$Cancer_Type_Glioma[!(column.annotations$Cancer_Type_Glioma %in% 
                                            c("Astro", "Oligo", "Glioblastoma", "OligoAstro", "AnaplasticAstro", "DiffuseAstro", "Pediatric"))] <- "glioma_other"
column.annotations$Cancer_Type_Glioma[column.annotations$Cancer_Type_Glioma %in% c("AnaplasticAstro", "DiffuseAstro")] <- "Astro"
column.annotations$Cancer_Type_Glioma[column.annotations$Cancer_Type_Glioma %in% c("Oligo", "OligoAstro")] <- "Oligo"
column.annotations$Cancer_Type_Broaddd[column.annotations$Cancer_Type_Broaddd %in% c("Craniopharyngioma", "Paraganglioma",  "Pineal_parenchymal_tumor", "Chondrosarcoma",
                                                                                     "Schwannoma", "Chordoma", "Schwannoma check", "Neuroblastoma", "Lymphoma",
                                                                                     "Medulloblastoma", "Pituitary_other", "Ependymoma", "Paraganglioglioma", "SFT")] <- "other"
column.annotations$Cancer_Type_Broaddd_Glioma <- column.annotations$Cancer_Type_Broaddd
column.annotations$Cancer_Type_Broaddd_Glioma[column.annotations$Cancer_Type_Broad == "Glioma"] <- column.annotations$Cancer_Type_Glioma[column.annotations$Cancer_Type_Broad == "Glioma"]

column.annotations$Cancer_Type_Broad_Glioma <- column.annotations$Cancer_Type_Broad
column.annotations$Cancer_Type_Broad_Glioma[column.annotations$Cancer_Type_Broad == "Glioma"] <- column.annotations$Cancer_Type_Glioma[column.annotations$Cancer_Type_Broad == "Glioma"]
column.annotations$Cancer_Type_Glioma_Grade <- paste(column.annotations$Cancer_Type_Glioma, column.annotations$Grade, sep = "_")
column.annotations$Cancer_Type_Glioma_Grade[!(column.annotations$Cancer_Type_Glioma_Grade %in% c("Astro_2", "Astro_3", "Glioblastoma_4", 
                                                                                                 "glioma_other_", "Oligo_2", "Oligo_3", "Pediatric_"))] <- "other"
## re order based on nmf order
column.annotations <- column.annotations[order(match(column.annotations$SAMPLE_ACCESSION_NBR, colnames(nmf.input))), ]

write.table(nmf.input, "../Upload to xchip/all_samples_nmf_input.txt", sep = "\t")

testing <- read.delim("../Upload to xchip/all_samples_nmf_input.txt", stringsAsFactors = FALSE)
## if multiple runs
nmf.consensus <- nmf(nmf.input, 1:8)
nmf.output.4 <- nmf(nmf.input, 2, nrun = 2)
coefmap(nmf.output.4, Colv = "consensus")


## read in nmf data generated from broad server
nmf.glioma.broad <- readRDS("c:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/Analysis/6.12-6.22/NMF_OUPUT_4_200.RDS")


## if single run
nmf.output.glioma.4 <- nmf(nmf.input, 4, seed = 123456789)
tiff("../Analysis/NMF K4 heatmap.tiff", width = 800, height = 800)
pdf("../Analysis/NMF K4 heatmap.pdf")
coefmap(minfit(nmf.glioma.broad), annCol = column.annotations$Cancer_Type_Broaddd, annColors = list(c("red", "yellow", "blue", "green4", "purple", "orange")), Colv = "basis")
dev.off()

## calculate overlap
classifier <- nmf.glioma.broad@fit@H
classifier <- as.data.frame(t(classifier))
classifier$subtype <- column.annotations$Cancer_Type_Broaddd_Glioma
classifier$basis <- sapply(1:nrow(classifier), function(x){which(classifier[x, 1:(ncol(classifier) - 1)] == max(classifier[x, 1:(ncol(classifier) -1)]))[1]}, 
                           USE.NAMES = FALSE)
classifier$subtype[classifier$subtype == ""] <- "unlabeled"
classifier$sample <- column.annotations$SAMPLE_ACCESSION_NBR

## percentage based classifier
percent.classifier <- table(classifier$basis, classifier$subtype)
for (i in 1:ncol(percent.classifier)){
    percent.classifier[, i] <- percent.classifier[, i] / sum(percent.classifier[, i])
}
percent.classifier <- as.data.frame(percent.classifier)
colnames(percent.classifier) <- c("basis", "subtype", "percent")
percent.classifier$basis <- paste("factor", percent.classifier$basis, sep = "_")


pdf("NMF K4 summary by basis.pdf", width = 10, height = 7)
colors <- distinctColorPalette(length(unique(classifier$subtype)))
names(colors) <- unique(classifier$subtype)
ggplot(data = classifier, aes(x = basis, fill = factor(subtype))) + geom_bar(position = "dodge") + scale_fill_manual(values = colors)
dev.off()

pdf("NMF K4 summary by basis percent.pdf", width = 10, height = 7)
ggplot(data = percent.classifier, aes(x = basis, y = percent, fill = subtype)) + geom_bar(stat = "Identity", position = "dodge") + scale_fill_manual(values = colors)
dev.off()



pdf("NMF Glioma K4 summary by subtype with grade.pdf", width = 10, height = 7)
pdf("NMF K4 summary by subtype flipped.pdf")
colors.factor <- distinctColorPalette(length(unique(classifier$basis)))
ggplot(data = classifier[!(classifier$subtype == "Oligo" & classifier$basis != 1), ], aes(x = subtype, fill = factor(basis))) + geom_bar(position = "dodge") + coord_flip() + rameen_theme + scale_fill_manual(values = colors.factor)
+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
+ coord_flip()
dev.off()

pdf("NMF Glioma K4 summary by subtype percent.pdf", width = 10, height = 7)
tiff("NMF K4 summary by subtype percent.tiff")
ggplot(data = percent.classifier, aes(x = subtype, y = percent, fill = basis)) + geom_bar(stat = "Identity", position = "dodge") 
+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()


## grey out factors that aren't significant enriched
incidence.table <- table(classifier$subtype, classifier$basis)
grey.table <- incidence.table

for (i in 1:nrow(incidence.table)){
    total <- sum(incidence.table[i, ])
    base.prop <- matrix(c(total/4, total*3/4), ncol = 2)
    grey.idx <- rep(1, ncol(incidence.table))
    
    for (j in 1: ncol(incidence.table)){
        current.prop <- matrix(c(incidence.table[i, j], total - incidence.table[i, j]), ncol = 2)
        if (prop.test(rbind(current.prop, base.prop))$p.value < 0.05 & current.prop[1,1] > base.prop[1,1]){
            grey.idx[j] <- 0
        }
    }
    
    grey.table[i, ] <- grey.idx
}

classifier.grey <- classifier

for (i in 1:nrow(grey.table)){
    for (j in 1:ncol(grey.table)){
        if (grey.table[i, j]){
            classifier.grey$basis[classifier.grey$basis == colnames(grey.table)[j] & classifier.grey$subtype == rownames(grey.table)[i]] <- 
                classifier.grey$basis[classifier.grey$basis == colnames(grey.table)[j] & classifier.grey$subtype == rownames(grey.table)[i]] + 0.5
        }
    }
}

classifier.grey$basis <- factor(classifier.grey$basis, levels = sort(unique(classifier.grey$basis)))

pdf("NMF K4 summary by subtype flipped grey.pdf")
colors.grey <- distinctColorPalette(4)
colors.grey <- c(colors.grey, colors.grey)
names(colors.grey)[5:8] <- c(1.5, 2.5, 3.5, 4.5)
colors.grey[5:8] <- "grey"

ggplot(data = classifier.grey, aes(x = subtype, fill = factor(basis))) + geom_bar(position = "dodge") + coord_flip() + rameen_theme + scale_fill_manual(values = colors.grey)
+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
+ coord_flip()
dev.off()



## plot percentage of samples within a given cancer type residing in most dominant factor
max.values <- data.frame(sapply(as.character(unique(percent.classifier$subtype)), function(x){max(percent.classifier$percent[percent.classifier$subtype == x])}))
max.values$subtype <- rownames(max.values)
colnames(max.values) <- c("percent", "subtype")
pdf("NMF CNVs K4 summary max basis.pdf", width = 10, height = 7)
ggplot(data = max.values, aes(x = subtype, y = percent, fill = subtype)) + geom_bar(stat = "identity") + scale_fill_manual(values = colors) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()


## plot contributing factors
tiff("../Analysis/NMF K4 components.tiff")
basismap(nmf.glioma.broad)
dev.off()
## show relationship between input variables and metagenes, with each variable coded by the dominant metagene it contributes to


## show a consensus map of the clustering based on multiple iterations
temp <- readRDS("NMF_OUPUT_1-15.RDS")
tiff("Consensus_map.jpg",  = 10000, height = 10000)
consensusmap(temp, annCol = column.annotations$Cancer_Type_Broad)
dev.off()

## multiple ranks at multiple runs each time
output.2_10 <- nmf(nmf.input, 2:10, nrun = 10)

consensusmap(output)




## pull out specific tumor types, order based on which cluster they were assigne

## meningioma vs metastasis
input.order <- c(classifier[classifier$basis == 2 & classifier$subtype == "Meningioma", "sample"] ,classifier[classifier$basis == 2 & classifier$subtype == "Metastasis", "sample"], 
                 classifier[classifier$basis == 4 & classifier$subtype == "Meningioma", "sample"], classifier[classifier$basis == 4 & classifier$subtype == "Metastasis", "sample"])

## GBM traditional vs other
input.order <- c(classifier[classifier$basis == 3 & classifier$subtype == "Glioblastoma", "sample"], classifier[classifier$basis != 3 & classifier$subtype == "Glioblastoma", "sample"])


PlotComut(
    mut.maf1 = alterations.trimmed,
    broad.cnv.maf = clustering.cnvs.trimmed[rownames(clustering.cnvs.trimmed) != "cnv.middle", ],
    samples = input.order,
    input.samples = "SAMPLE_ACCESSION_NBR", 
    input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification", 
    dimensions = c(18,10),
    fixed.order = input.order,
    col.vector = c("frameshift_indel" = "red", "missense" ="red" , "nonsense" = "red", "splice_site" = "red", "stop_codon" = "red", "in_frame_indel" = "red",
                   "other" = "red", "arm-level gain" = "red", "arm-level loss" = "red", "2DEL" = "red", "HA" = "red", "indel" = "red", "rearrangement" = "red",
                   "Glioblastoma" = "black", "AnaplasticAstro" = "green", "Angiocentric" = "orange", "Astro" = "green", "altered" = "red", 
                   "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                   "0" = "gray90", "1" = "gray70", "2" = "gray50", "3" = "gray25", "4" = "gray1")
    ,phenotypes = clustering.clinical[, -3]
)


