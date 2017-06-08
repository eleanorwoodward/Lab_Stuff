## Clustering analysis

## generate distances matrix using dist(). Distances are computed between rows: need to transpose if samples in columns


matrix <- rbind(c(1, rep(0,7),1,1), c(rep(0,10)), c(rep(0,10)), c(0, rep(1,7),0,0), c(1, rep(0,7),1,1))
chromosomes <- paste("chr", 1:10, "p", sep = "")
rownames(matrix) <- c(paste("sample", 1:5, sep =  ""))
colnames(matrix) <- chromosomes

d <- dist(matrix)
clust <- hclust(d)
plot(clust)
matrix.df <- data.frame(matrix)
matrix.df$samples <-  rownames(matrix)

samples <- rownames(matrix)
long.df <- melt(matrix.df, id = "samples")

long.df$samples <- factor(long.df$samples, levels = clust$labels)
long.df$value[abs(long.df$value) < 3] <- NA
long.df$value[long.df$value > 3] <- 1
long.df$value[long.df$value < -3] <- -1

mut <- ggplot(long.df, aes(x=samples, y=variable, height=0.9, width=0.9)) + 
    geom_tile(aes(fill=value)) +
    scale_fill_gradient(low = "blue", high = "red", na.value = "white", limits = c(-1, 1)) +
    xlab("Subject") +
    ggtitle("Trial Heatmap") +
    theme(
        legend.key = element_rect(fill='NA'),
        legend.key.size = unit(1, 'cm'),
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
    )
print(mut) 


## create legend
legend.df <- data.frame(rep(samples, 2), c(rep("tumor", length(samples)), rep("sex", length(samples))), 
                        c(rep(c("glioma", "meingioma"), length(samples)/2), rep(c("male", "female"), length(samples)/2)))
colnames(legend.df) <- c("samples", "phenotype", "values")
l <- ggplot(legend.df, aes(x=samples, y = phenotype)) + coord_fixed(ratio = .2) + geom_tile(aes(fill = values))


fd=data.frame(x = rep(c("x","y","z"),3), 
              y=c("a","b","c","b","c","a","c","a","b"),
              z=c(0,1,0,1,1,1,0,0,1))

# plot
p <- ggplot(fd, aes(x, y)) + geom_tile(aes(fill = z)) + 
scale_fill_gradient(low = "white",high = "steelblue", limits=c(0,1)) + 
theme_grey() + 
labs(x = "", y= "") + 
coord_fixed(ratio= 1)
scale_x_discrete(expand = c(0, 0)) + 
scale_y_discrete(expand = c(0, 0)) + 
theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size=12, angle=90, hjust=0, colour="black"))


## use NMF to cluster
## adapted from: https://cran.r-project.org/web/packages/NMF/vignettes/heatmaps.pdf and 

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
nmf.input <- reshape(clustering.input.gliomas, idvar = "samples", timevar = "genes", direction = "wide")

nmf.input[nmf.input == "wt"] <- 0
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


## if multiple runs
nmf.output.4 <- nmf(nmf.input, 2, nrun = 2)
coefmap(nmf.output.4, Colv = "consensus")

## if single run
nmf.output.glioma.4 <- nmf(nmf.input, 4, seed = 123456789)
tiff("../Analysis/NMF Glioma K4 heatmap.tiff", width = 1000, height = 1000)
pdf("../Analysis/NMF K4 heatmap.pdf")
coefmap(minfit(nmf.output.glioma.4), annCol = column.annotations$Cancer_Type_Glioma, annColors = list(c("red", "yellow", "blue", "green4", "purple", "orange")), Colv = "basis")
dev.off()

## calculate overlap
classifier <- nmf.output.glioma.4@fit@H
classifier <- as.data.frame(t(classifier))
classifier$subtype <- column.annotations$Cancer_Type_Glioma_Grade
classifier$basis <- sapply(1:nrow(classifier), function(x){which(classifier[x, 1:(ncol(classifier) - 1)] == max(classifier[x, 1:(ncol(classifier) -1)]))[1]}, 
                           USE.NAMES = FALSE)
classifier$subtype[classifier$subtype == ""] <- "unlabeled"

## percentage based classifier
percent.classifier <- table(classifier$basis, classifier$subtype)
for (i in 1:ncol(percent.classifier)){
    percent.classifier[, i] <- percent.classifier[, i] / sum(percent.classifier[, i])
}
percent.classifier <- as.data.frame(percent.classifier)
colnames(percent.classifier) <- c("basis", "subtype", "percent")
percent.classifier$basis <- paste("factor", percent.classifier$basis, sep = "_")


pdf("NMF Glioma K4 summary by basis.pdf", width = 10, height = 7)
colors <- distinctColorPalette(length(unique(classifier$subtype)))
names(colors) <- unique(classifier$subtype)
ggplot(data = classifier, aes(x = basis, fill = factor(subtype))) + geom_bar(position = "dodge") + scale_fill_manual(values = colors)
dev.off()

pdf("NMF Glioma K4 summary by basis percent.pdf", width = 10, height = 7)
ggplot(data = percent.classifier, aes(x = basis, y = percent, fill = subtype)) + geom_bar(stat = "Identity", position = "dodge") + scale_fill_manual(values = colors)
dev.off()



pdf("NMF Glioma K4 summary by subtype with grade.pdf", width = 10, height = 7)
ggplot(data = classifier, aes(x = subtype, fill = factor(basis))) + geom_bar(position = "dodge") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
+ coord_flip()
dev.off()

pdf("NMF Glioma K4 summary by subtype percent.pdf", width = 10, height = 7)
ggplot(data = percent.classifier, aes(x = subtype, y = percent, fill = basis)) + geom_bar(stat = "Identity", position = "dodge") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
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
tiff("../Analysis/NMF CNVs K4 components.tiff", width = 1000, height = 1000)
basismap(nmf.output.cnvs.4)
dev.off()
## show relationship between input variables and metagenes, with each variable coded by the dominant metagene it contributes to


## show a consensus map of the clustering based on multiple iterations
consensusmap(output, annCol = column.annotations$Cancer_Type_Specific)


## multiple ranks at multiple runs each time
output.2_10 <- nmf(nmf.input, 2:10, nrun = 10)

consensusmap(output)
