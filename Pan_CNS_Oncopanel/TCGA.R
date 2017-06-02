## analyze TCGA mutation data

## first import data
## reformat data for ggplot
TCGA.LGG.focal <- read.delim("../TCGA Data/LGG_focal_CNV.txt", stringsAsFactors = FALSE)
TCGA.LGG.focal <- t(TCGA.LGG.focal)
colnames(TCGA.LGG.focal) <- as.character(TCGA.LGG.focal[1, ])
TCGA.LGG.focal <- as.data.frame(TCGA.LGG.focal)
TCGA.LGG.focal$Tumor_Sample_Barcode <- rownames(TCGA.LGG.focal)
TCGA.LGG.focal <- TCGA.LGG.focal[-c(1:3), ]
TCGA.LGG.focal <- melt(TCGA.LGG.focal, id.vars = "Tumor_Sample_Barcode")

TCGA.GBM.focal <- read.delim("../TCGA Data/GBM_focal_CNV.txt", stringsAsFactors = FALSE)
TCGA.GBM.focal <- t(TCGA.GBM.focal)
colnames(TCGA.GBM.focal) <- as.character(TCGA.GBM.focal[1, ])
TCGA.GBM.focal <- as.data.frame(TCGA.GBM.focal)
TCGA.GBM.focal$Tumor_Sample_Barcode <- rownames(TCGA.GBM.focal)
TCGA.GBM.focal <- TCGA.GBM.focal[-c(1:3), ]
TCGA.GBM.focal <- melt(TCGA.GBM.focal, id.vars = "Tumor_Sample_Barcode")

## broad data
TCGA.LGG.broad <- read.delim("../TCGA Data/LGG_broad_CNV.txt", stringsAsFactors = FALSE)
TCGA.LGG.broad <- t(TCGA.LGG.broad)
colnames(TCGA.LGG.broad) <- as.character(TCGA.LGG.broad[1, ])
TCGA.LGG.broad <- as.data.frame(TCGA.LGG.broad)
TCGA.LGG.broad$Tumor_Sample_Barcode <- rownames(TCGA.LGG.broad)
TCGA.LGG.broad <- TCGA.LGG.broad[-1, ]
TCGA.LGG.broad <- melt(TCGA.LGG.broad, id.vars = "Tumor_Sample_Barcode")


TCGA.GBM.broad <- read.delim("../TCGA Data/GBM_broad_CNV.txt", stringsAsFactors = FALSE)
TCGA.GBM.broad <- t(TCGA.GBM.broad)
colnames(TCGA.GBM.broad) <- as.character(TCGA.GBM.broad[1, ])
TCGA.GBM.broad <- as.data.frame(TCGA.GBM.broad)
TCGA.GBM.broad$Tumor_Sample_Barcode <- rownames(TCGA.GBM.broad)
TCGA.GBM.broad <- TCGA.GBM.broad[-1, ]
TCGA.GBM.broad <- melt(TCGA.GBM.broad, id.vars = "Tumor_Sample_Barcode")




## mutations data
TCGA.LGG.mut <- read.delim("../TCGA Data/LGG_Mutations.txt", stringsAsFactors = FALSE)
TCGA.LGG.mut <- melt(TCGA.LGG.mut[, -c(2:6)], id.vars = "Tumor")
TCGA.LGG.mut <- TCGA.LGG.mut[TCGA.LGG.mut$value == "Mutant", ]

TCGA.GBM.mut <- read.delim("../TCGA Data/GBM_Mutations.txt", stringsAsFactors = FALSE, comment.char = "#")



## DFCI samples
gbms <- master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Specific == "Glioblastoma"]
gliomas <- master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Broad == "Glioma" & !(master.sheet$Cancer_Type_Specific %in% c("Glioblastoma", "Pediatric"))]


## compare rates of specific alterations
## Mutations for GBM
genes <- names(sort(table(all.mutations.tier1.3$BEST_EFF_GENE), decreasing = TRUE))[1:20]
short.TCGA.GBM.mut <- TCGA.GBM.mut[TCGA.GBM.mut$Hugo_Symbol %in% genes, c("Tumor_Sample_Barcode", "Hugo_Symbol")]

DFCI.GBM.mut <- all.mutations.tier1.4[all.mutations.tier1.4$SAMPLE_ACCESSION_NBR %in% master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Specific == "Glioblastoma"], ]
short.DFCI.GBM.mut <- DFCI.GBM.mut[DFCI.GBM.mut$BEST_EFF_GENE %in% genes, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE")]
colnames(short.TCGA.GBM.mut) <- c("SAMPLE_ACCESSION_NBR", "GENE")
colnames(short.DFCI.GBM.mut) <- c("SAMPLE_ACCESSION_NBR", "GENE")

## Mutations for LGG
genes <- names(sort(table(all.mutations.tier1.3$BEST_EFF_GENE), decreasing = TRUE))[1:20]
short.TCGA.LGG.mut <- TCGA.LGG.mut[TCGA.LGG.mut$variable %in% genes, ]

DFCI.LGG.mut <- all.mutations.tier1.4[all.mutations.tier1.4$SAMPLE_ACCESSION_NBR %in% master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Broad == "Glioma"], ]
DFCI.LGG.mut <- DFCI.LGG.mut[DFCI.LGG.mut$SAMPLE_ACCESSION_NBR %in% master.sheet$SAMPLE_ACCESSION_NBR[!(master.sheet$Cancer_Type_Specific %in% c("Glioblastoma", "Pediatric"))], ]
short.DFCI.LGG.mut <- DFCI.LGG.mut[DFCI.LGG.mut$BEST_EFF_GENE %in% genes, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE")]
short.TCGA.LGG.mut <- short.TCGA.LGG.mut[, -3]
colnames(short.TCGA.LGG.mut) <- c("SAMPLE_ACCESSION_NBR", "GENE")
colnames(short.DFCI.LGG.mut) <- c("SAMPLE_ACCESSION_NBR", "GENE")


## focal CNVs GBM
genes.focal <- names(sort(table(all.cnv.high$GENE), decreasing = TRUE))[1:10]
short.TCGA.GBM.focal <- TCGA.GBM.focal[TCGA.GBM.focal$value %in% c(-2, " 2", "2"), ]
short.TCGA.GBM.focal <- short.TCGA.GBM.focal[short.TCGA.GBM.focal$variable %in% genes.focal, 1:2]

DFCI.GBM.focal <- all.cnv.high[all.cnv.high$SAMPLE_ACCESSION_NBR %in% master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Specific == "Glioblastoma"], ]
short.DFCI.GBM.focal <- DFCI.GBM.focal[DFCI.GBM.focal$GENE %in% genes.focal, c("SAMPLE_ACCESSION_NBR", "GENE")]
colnames(short.TCGA.GBM.focal) <- c("SAMPLE_ACCESSION_NBR", "GENE")
colnames(short.DFCI.GBM.focal) <- c("SAMPLE_ACCESSION_NBR", "GENE")


## focal CNVs LGG 
genes.focal <- names(sort(table(all.cnv.high$GENE), decreasing = TRUE))[1:10]
short.TCGA.LGG.focal <- TCGA.LGG.focal[TCGA.LGG.focal$value %in% c(-2, " 2", "2"), ]
short.TCGA.LGG.focal <- short.TCGA.LGG.focal[short.TCGA.LGG.focal$variable %in% genes.focal, 1:2]

DFCI.LGG.focal <- all.cnv.high[all.cnv.high$SAMPLE_ACCESSION_NBR %in% master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Broad == "Glioma"], ]
DFCI.LGG.focal <- DFCI.LGG.focal[DFCI.LGG.focal$SAMPLE_ACCESSION_NBR %in% master.sheet$SAMPLE_ACCESSION_NBR[!(master.sheet$Cancer_Type_Specific %in% c("Glioblastoma", "Pediatric"))], ]
short.DFCI.LGG.focal <- DFCI.LGG.focal[DFCI.LGG.focal$GENE %in% genes.focal, c("SAMPLE_ACCESSION_NBR", "GENE")]
colnames(short.TCGA.LGG.focal) <- c("SAMPLE_ACCESSION_NBR", "GENE")
colnames(short.DFCI.LGG.focal) <- c("SAMPLE_ACCESSION_NBR", "GENE")


## Broad CNVs
## first fix DFCI samples
DFCI.broad <- t(all.cnvs.broad)
DFCI.broad <- as.data.frame(DFCI.broad)
DFCI.broad$SAMPLE_ACCESSION_NBR <- rownames(DFCI.broad)
DFCI.broad <- melt(DFCI.broad, id.vars = "SAMPLE_ACCESSION_NBR")
DFCI.GBM.broad <- DFCI.broad[DFCI.broad$SAMPLE_ACCESSION_NBR %in% gbms, ]
DFCI.LGG.broad <- DFCI.broad[DFCI.broad$SAMPLE_ACCESSION_NBR %in% gliomas, ]

## Broad CNV GBM
genes.broad <- c("1p", "7q", "7p", "10q", "19q")
short.TCGA.GBM.broad <- TCGA.GBM.broad[abs(as.numeric(TCGA.GBM.broad$value)) > 0.3, ]
short.TCGA.GBM.broad <- short.TCGA.GBM.broad[short.TCGA.GBM.broad$variable %in% genes.broad, 1:2]

short.DFCI.GBM.broad <- DFCI.GBM.broad[DFCI.GBM.broad$variable %in% genes.broad, ]
short.DFCI.GBM.broad <- short.DFCI.GBM.broad[short.DFCI.GBM.broad$value != 0, 1:2]
colnames(short.TCGA.GBM.broad) <- c("SAMPLE_ACCESSION_NBR", "GENE")
colnames(short.DFCI.GBM.broad) <- c("SAMPLE_ACCESSION_NBR", "GENE")

## Broad CNV LGG
genes.broad <- c("1p", "7q", "7p", "10q", "19q")
short.TCGA.LGG.broad <- TCGA.LGG.broad[abs(as.numeric(TCGA.LGG.broad$value)) > 0.3, ]
short.TCGA.LGG.broad <- short.TCGA.LGG.broad[short.TCGA.LGG.broad$variable %in% genes.broad, 1:2]

short.DFCI.LGG.broad <- DFCI.LGG.broad[DFCI.LGG.broad$variable %in% genes.broad, ]
short.DFCI.LGG.broad <- short.DFCI.LGG.broad[short.DFCI.LGG.broad$value != 0, 1:2]
colnames(short.TCGA.LGG.broad) <- c("SAMPLE_ACCESSION_NBR", "GENE")
colnames(short.DFCI.LGG.broad) <- c("SAMPLE_ACCESSION_NBR", "GENE")



## generate table and collapse to dataframe for plotting
TCGA.input <- short.TCGA.LGG.broad
TCGA.original <- TCGA.LGG.broad
DFCI.input <- short.DFCI.LGG.broad
DFCI.original <- DFCI.LGG.broad
genes <- genes.broad

TCGA.df <- data.frame(genes, as.numeric(table(TCGA.input$GENE)[genes]), stringsAsFactors = FALSE)
colnames(TCGA.df)[2] <- "count"
TCGA.df$count <- TCGA.df$count / length(unique(TCGA.original$Tumor_Sample_Barcode))

DFCI.df <- data.frame(genes, as.numeric(table(DFCI.input$GENE)[genes]), stringsAsFactors = FALSE)
colnames(DFCI.df)[2] <- "count"
DFCI.df$count <-DFCI.df$count / length(unique(DFCI.original$SAMPLE_ACCESSION_NBR))

TCGA.df$cohort <- "TCGA"
DFCI.df$cohort <- "DFCI"
short.combined <- rbind(TCGA.df, DFCI.df)
pdf("TCGA Comparison LGG Broad CNV.pdf")
ggplot(data = short.combined, aes(x = genes, y = count, fill = cohort)) + geom_bar(position = "dodge", stat = "identity") + rameen_theme
dev.off()






