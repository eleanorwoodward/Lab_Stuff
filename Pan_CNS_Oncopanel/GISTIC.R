## Comparison of aCGH data
acgh <- read.delim("../GISTIC/broad_values_by_arm.txt", stringsAsFactors = FALSE)
colnames(acgh)[-1] <- gsub("US83003520_", "", colnames(acgh)[-1])
colnames(acgh)[-1] <- gsub("_[A-z0-9]*_[A-z]*_[0-9]*_[A-z0-9]*.txt", "", colnames(acgh)[-1])
annotation <- read.delim("../GISTIC/ACGH_Sample_Key.txt", stringsAsFactors = FALSE)
annotation <- annotation[, 1:7]

annotation <- merge(annotation, master.sheet[, c("SAMPLE_ACCESSION_NBR", "MRN")], all.x = TRUE)
annotation <- annotation[!duplicated(annotation$MRN), ]
annotation <- annotation[annotation$Barcode %in% colnames(acgh) & !is.na(annotation$SAMPLE_ACCESSION_NBR), ]
acgh.comparison <- acgh[, c(TRUE, colnames(acgh)[-1] %in% annotation$Barcode)]

acgh.comparison <- acgh.comparison[, c(TRUE, colnames(acgh.comparison)[-1] %in% annotation$Barcode)]

acgh.comparison <- t(acgh.comparison)
colnames(acgh.comparison) <- acgh.comparison[1, ]
acgh.comparison <- acgh.comparison[-1, ]

cnv.comparison <- t(all.cnvs.broad)
cnv.comparison <- cnv.comparison[rownames(cnv.comparison) %in% annotation$SAMPLE_ACCESSION_NBR, ]

## reorder acgh comparison based on order of 
idx <- match(rownames(cnv.comparison), annotation$SAMPLE_ACCESSION_NBR)
annotation <- annotation[idx, ]
idx.2 <- match(annotation$Barcode, rownames(acgh.comparison))
acgh.comparison <- acgh.comparison[idx.2, ]
acgh.comparison <- acgh.comparison[, colnames(acgh.comparison) %in% colnames(cnv.comparison)]

acgh.comparison<- apply(acgh.comparison, 2, as.numeric)

acgh.comparison[acgh.comparison < -.3] <- -1
acgh.comparison[acgh.comparison > .3] <- 1
acgh.comparison[abs(acgh.comparison) != 1] <- 0
percent <- c()
for(i in 1:ncol(acgh.comparison)){
    tbl <- table(acgh.comparison[, i], cnv.comparison[, i])
    agree <- sum(tbl[1,1], tbl[2,2], tbl[3,3])
    total <- sum(tbl[1:length(tbl)])
    percent <- c(percent, agree / total)
    
}

percent.df <- data.frame(colnames(acgh.comparison), percent)
ggplot(data = percent.df, aes(x = colnames.acgh.comparison., y = percent)) + geom_bar(stat = "identity")


