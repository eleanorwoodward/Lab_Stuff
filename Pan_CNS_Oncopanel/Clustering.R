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
alterations <- all.mutations.tier1.4[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
alterations.cnv <- all.cnv.high[, c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD")]
colnames(alterations.cnv) <- colnames(alterations)
alterations.svs <- all.svs.formatted[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
alterations <- rbind(alterations, alterations.cnv, alterations.svs)
alterations$BEST_EFF_GENE[alterations$BEST_EFF_GENE %in% c("CDKN2A", "CDKN2B")] <- "CDKN2A/B"
alterations$BEST_EFF_GENE[alterations$BEST_EFF_GENE == "MYCN"] <- "MYC"
alterations$variant_classification <- "altered"
alterations.trimmed <- alterations[alterations$BEST_EFF_GENE %in% c("TP53", "IDH1", "NF2", "MYC", "KMT2D", "PTPRD", "PRKDC", "KRAS", 
                                                                    "PTCH1", "PTEN", "EGFR", "ATRX", "CDKNA/B", "NF1", "SMARC4A"), ]
names <- unique(master.sheet$SAMPLE_ACCESSION_NBR)
quiet.idx <- sapply(names, function(x){if(sum(alterations.trimmed$SAMPLE_ACCESSION_NBR == x) == 0){return(TRUE)}else{return(FALSE)}}, USE.NAMES = FALSE)
quiet.names <- names[quiet.idx]
quiet.samples <- cbind(quiet.names, "quiet.mut", "altered")
colnames(quiet.samples) <- colnames(alterations.trimmed)
alterations.trimmed <- rbind(alterations.trimmed, quiet.samples)

## feature selection time! consolidate arm-level CNVs
clustering.cnvs <- all.cnvs.broad
clustering.cnvs["count", ] <- as.numeric(apply(clustering.cnvs, 2, function(x){sum(x != 0)}))
clustering.cnvs["1p.19q", clustering.cnvs["1p", ] == -1 & clustering.cnvs["19q", ] == -1] <- 1
clustering.cnvs["1p.19q", is.na(clustering.cnvs["1p.19q", ])] <- 0
clustering.cnvs["7pq", clustering.cnvs["7p", ] == 1 & clustering.cnvs["7q", ] == 1] <- 1
clustering.cnvs["7pq", is.na(clustering.cnvs["7pq", ])] <- 0
clustering.cnvs["disrupted", clustering.cnvs["count", ] > 5 ] <- 1
clustering.cnvs["disrupted", is.na(clustering.cnvs["disrupted", ])] <- 0
clustering.cnvs["quiet", clustering.cnvs["count", ] == 0] <- 1
clustering.cnvs["quiet", is.na(clustering.cnvs["quiet", ])] <- 0


## isolate only direction we're interested in for specific chromosomes
clustering.cnvs["1p", clustering.cnvs["1p", ] == 1] <- 0
clustering.cnvs["1p", clustering.cnvs["1p", ] == -1] <- 1

clustering.cnvs["10q", clustering.cnvs["10q", ] == 1] <- 0
clustering.cnvs["10q", clustering.cnvs["10q", ] == -1] <- 1

clustering.cnvs["17p", clustering.cnvs["17p", ] == 1] <- 0
clustering.cnvs["17p", clustering.cnvs["17p", ] == -1] <- 1

clustering.cnvs["22q", clustering.cnvs["22q", ] == 1] <- 0
clustering.cnvs["22q", clustering.cnvs["22q", ] == -1] <- 1
clustering.cnvs.trimmed <- clustering.cnvs[c("1p","7pq", "10q", "17p", "1p.19q", "22q", "disrupted", "quiet"), ]


clustering.clinical <- master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad", "exclude")]
clustering.clinical$Cancer_Type_Broad[clustering.clinical$exclude == "bch"] <- "pediatric"

clustering.input <- PlotComut(
    mut.maf1 = alterations.trimmed,
    broad.cnv.maf = clustering.cnvs.trimmed,
    samples = master.sheet$SAMPLE_ACCESSION_NBR,
    input.samples = "SAMPLE_ACCESSION_NBR", 
    input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification", gene.cutoff = 60, 
    return.matrix = TRUE,
    dimensions = c(18,10), 
    col.vector = c("frameshift_indel" = "red", "missense" ="red" , "nonsense" = "red", "splice_site" = "red", "stop_codon" = "red", "in_frame_indel" = "red",
                   "other" = "red", "arm-level gain" = "red", "arm-level loss" = "red", "2DEL" = "red", "HA" = "red", "indel" = "red", "rearrangement" = "red",
                   "Glioblastoma" = "black", "AnaplasticAstro" = "green", "Angiocentric" = "orange", "Astro" = "green", "altered" = "red", 
                   "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                   "0" = "gray90", "1" = "gray70", "2" = "gray50", "3" = "gray25", "4" = "gray1")
    ,phenotypes = clustering.clinical[, -3]
)


## make fake data for testing
samples <- c(rep("A", 40), rep("B", 40), rep("C", 40))
row1 <- c(rbinom(40, 1, .8), rbinom(40, 1, .5), rbinom(40, 1, .2))
row2 <- c(rbinom(40, 1, .9), rbinom(40, 1, .5), rbinom(40, 1, .1))
row3 <- c(rbinom(40, 1, .9), rbinom(40, 1, .1), rbinom(40, 1, .9))
row4 <- c(rbinom(40, 1, .1), rbinom(40, 1, .9), rbinom(40, 1, .4))
row5 <- c(rbinom(40, 1, .9), rbinom(40, 1, .3), rbinom(40, 1, .9))
row6 <- c(rbinom(40, 1, .1), rbinom(40, 1, .1), rbinom(40, 1, .9))
test.matrix <- rbind(row1, row2, row3, row4, row5, row6)
test.matrix <- test.matrix[, colSums(test.matrix) != 0]

test.output <- nmf(test.matrix, 4)
aheatmap(test.matrix, annCol = samples, annColors = list(c("red", "blue", "green")))
coefmap(minfit(test.output), annCol = samples, Colv = "basis", annColors = list(c("red", "orange", "purple", "green", "black", "blue")))

## change to non-negative matrix, features as rows and samples as columns (still not sure why not reversed??)
nmf.input <- reshape(clustering.input, idvar = "samples", timevar = "genes", direction = "wide")

nmf.input[nmf.input == "wt"] <- 0
nmf.input[, -1][nmf.input[, -1] != 0] <- 1
nmf.input <- t(nmf.input)
colnames(nmf.input) <- nmf.input[1, ]
nmf.input <- nmf.input[-1, ]
nmf.input <- apply(nmf.input, 2, as.numeric)

## remove empty columns
nmf.input <- nmf.input[, colSums(nmf.input) != 0]
nmf.input <- nmf.input[-23, ]


column.annotations <- master.sheet
column.annotations <- column.annotations[column.annotations$SAMPLE_ACCESSION_NBR %in% colnames(nmf.input), ]
column.annotations$Cancer_Type_Broaddd <- column.annotations$Cancer_Type_Broad
column.annotations$Cancer_Type_Broaddd[column.annotations$Cancer_Type_Broaddd %in% c("Craniopharyngioma", "Paraganglioma",  "Pineal_parenchymal_tumor", "Chondrosarcoma",
                                                                                     "Schwannoma", "Chordoma", "Schwannoma check", "Neuroblastoma", "Lymphoma",
                                                                                     "Medulloblastoma", "Pituitary_other", "Ependymoma", "Paraganglioglioma", "SFT")] <- "other"
## re order based on nmf order
column.annotations <- column.annotations[order(match(column.annotations$SAMPLE_ACCESSION_NBR, colnames(nmf.input))), ]
## order columns by consensus classes defined from given number of iterations with set k value
aheatmap(nmf.input, annCol = column.annotations$Cancer_Type_Broaddd, annColors = list(c("red", "yellow", "blue", "green", "purple", "orange")))


## if multiple runs
nmf.output.4 <- nmf(nmf.input, 2, nrun = 2)
coefmap(nmf.output.4, Colv = "consensus")

## if single run
nmf.output.3 <- nmf(nmf.input, 3)
pdf("../Analysis/NMF K2 heatmap.pdf")
coefmap(minfit(nmf.output.2), annCol = column.annotations$Cancer_Type_Broaddd, annColors = list(c("red", "yellow", "blue", "green", "purple", "orange")), Colv = "basis")
dev.off()

## calculate overlap
classifier <- nmf.output.2@fit@H
classifier <- as.data.frame(t(classifier))
classifier$subtype <- column.annotations$Cancer_Type_Broaddd
classifier$group <- "A"
classifier$group[classifier$V2 > classifier$V1] <- "B"
pdf("NMF K2 summary.pdf")
ggplot(data = classifier, aes(x = subtype, fill = group)) + geom_bar(position = "dodge")
dev.off()

## show relationship between input variables and metagenes, with each variable coded by the dominant metagene it contributes to
basismap(output)


## show a consensus map of the clustering based on multiple iterations
consensusmap(output, annCol = column.annotations$Cancer_Type_Specific)


## multiple ranks at multiple runs each time
output.2_10 <- nmf(nmf.input, 2:10, nrun = 10)

consensusmap(output)
