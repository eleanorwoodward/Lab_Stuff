## rearrangements analysis


## generate gene list by duplicating gene1 and gene2
all.svs <- merge(all.svs, master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad")], all.x = TRUE)
all.svs.genes <- all.svs[, -7]
temp <- all.svs[, -6]
colnames(all.svs.genes)[6] <- "gene"
colnames(temp)[6] <- "gene"
all.svs.genes <- rbind(all.svs.genes, temp)
all.svs$Cancer_Type_Broad[all.svs$Cancer_Type_Broad == ""] <- "other"

## plot most frequently rearranged gene, 1d style
temp <- all.svs.genes[all.svs.genes$gene != "", ]
temp$gene <- factor(temp$gene, levels = names(sort(table(temp$gene), decreasing = TRUE)))
temp <- temp[temp$gene %in% levels(temp$gene)[1:30], ]
pdf("Rearrangement_1D_by_tumor_type.pdf")
ggplot(data = temp, aes(x=gene, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme
dev.off()


## calculate most frequently rearranged gene pairs

gene.rearrangements <- all.svs[all.svs$rearrangement == TRUE, ]
gene.rearrangements$Gene1 <- sapply(gene.rearrangements$Gene1, trimws, USE.NAMES = F)
gene.rearrangements$Gene1 <- sapply(gene.rearrangements$Gene1, trimws, USE.NAMES = F)

for (i in 1:nrow(gene.rearrangements)){
    gene.rearrangements[i, 6:7] <- sort(as.character(gene.rearrangements[i, 6:7]))
}
    
matches <- NULL
for (i in 1:nrow(gene.rearrangements)){
    # Gets first gene-gene pair
    current.gene <- gene.rearrangements$Gene1[i]
    current.target <- gene.rearrangements$Gene2[i]
    friends <- FilterMaf(gene.rearrangements[-i,], current.gene, "Gene1")
    
    #checks to see if any other instances of same starter gene have same target
    if (nrow(friends) > 0){
        idx <- friends$Gene2 %in% current.target
        close.friends <- friends[idx,]
        if (nrow(close.friends) > 0){
            temp <- c(rownames(gene.rearrangements)[i], current.gene, current.target, gene.rearrangements$SAMPLE_ACCESSION_NBR[i], gene.rearrangements$Cancer_Type_Broad[i])
            matches <- c(matches, temp)
            
        }
        
    } 
}

matches.table <-matrix(matches, length(matches)/5, 5, byrow = T)
pairwise.matches <- data.frame(matches.table, stringsAsFactors = FALSE)
colnames(pairwise.matches) <- c("Original Index", "Gene1", "Gene2", "sample", "Cancer_Type_Broad")
View(pairwise.matches)
pairwise.matches$name <- paste(pairwise.matches$Gene1, pairwise.matches$Gene2, sep = " - ")

unique(pairwise.matches[, 2:3])

temp <- pairwise.matches
temp$name <- factor(temp$name, levels = names(sort(table(temp$name), decreasing = TRUE)))
pdf("Rearrangement_2D_by_tumor_type.pdf")
ggplot(data = temp, aes(x=name, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme
dev.off()

