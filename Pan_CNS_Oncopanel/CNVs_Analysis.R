## Oncopanel CNV analysis

## looks genes most frequently suffering high-level gains or complete loss
all.cnvs <- merge(all.cnvs, master.sheet[, c("Cancer_Type_Broad", "SAMPLE_ACCESSION_NBR")], "SAMPLE_ACCESSION_NBR")
all.cnv.high <- subset(all.cnvs, CNV_TYPE_CD %in% c("2DEL", "HA"))
temp <- all.cnv.high
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all high-level CNVs.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = CNV_TYPE_CD)) + geom_bar() + rameen_theme
dev.off()


## look just at deletions, colored by tumor type
temp <- subset(all.cnv.high, CNV_TYPE_CD == "2DEL")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all high-level CNVs deletions.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme + ggtitle(label = "2 Copy deletions")
dev.off()


## look just at amplifications, colored by tumor type
temp <- subset(all.cnv.high, CNV_TYPE_CD == "HA")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
temp <- merge(temp, master.sheet[, c("Cancer_Type_Broad", "SAMPLE_ACCESSION_NBR")], "SAMPLE_ACCESSION_NBR")
pdf("all high-level CNVs amplifications.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme +ggtitle(label = "High level amplifications")
dev.off()

temp <- all.cnv.high
temp$Cancer_Type_Broad <- factor(temp$Cancer_Type_Broad, levels = names(sort(table(temp$Cancer_Type_Broad), decreasing = T)))
pdf("all high-level CNVs by tumor.pdf", width = 14)
ggplot(data = temp, aes(x = Cancer_Type_Broad, fill = CNV_TYPE_CD)) + geom_bar(position = "dodge") + rameen_theme +
    ggtitle(label = "high level CNVs stratified by tumor type")
dev.off()




## look at low level amplifications and deletions
all.cnv.low <- subset(all.cnvs, CNV_TYPE_CD %in% c("1DEL", "LA"))

temp <- all.cnv.low
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all low-level CNVs.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:100]),
       aes(x = GENE, fill = CNV_TYPE_CD)) + geom_bar() + rameen_theme
dev.off()


## look just at deletions, colored by tumor type
temp <- subset(all.cnv.low, CNV_TYPE_CD == "1DEL")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all low-level CNVs deletions.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme + ggtitle(label = "1 Copy deletions")
dev.off()


## look just at amplifications, colored by tumor type
temp <- subset(all.cnv.low, CNV_TYPE_CD == "LA")
temp$GENE <- factor(temp$GENE, levels = names(sort(table(temp$GENE), decreasing = T)))
pdf("all low-level CNVs amplifications.pdf", width = 14)
ggplot(data = subset(temp, GENE %in% levels(temp$GENE)[1:50]),
       aes(x = GENE, fill = Cancer_Type_Broad)) + geom_bar() + rameen_theme +ggtitle(label = "Low level amplifications")
dev.off()

temp <- all.cnv.low
temp$Cancer_Type_Broad <- factor(temp$Cancer_Type_Broad, levels = names(sort(table(temp$Cancer_Type_Broad), decreasing = T)))
pdf("all low-level CNVs by tumor.pdf", width = 14)
ggplot(data = temp, aes(x = Cancer_Type_Broad, fill = CNV_TYPE_CD)) + geom_bar(position = "dodge") + rameen_theme +
    ggtitle(label = "low level CNVs stratified by tumor type")
dev.off()


## Generate arm-level CNV data
gene.to.arm <- read.delim("../OncDRS data/Gene_to_arm_level.txt", stringsAsFactors = F)

## remove x + y chromosomes, along with those with fewer than 5 genes genes
gene.to.arm <- gene.to.arm[!(gene.to.arm$ID %in% c("0Xp", "0Xq", "0Yp", "5p", "10p", "20p", "21q")), ]
all.cnvs.broad <- matrix(nrow = length(unique(gene.to.arm$ID)), ncol = length(unique(master.sheet$SAMPLE_ACCESSION_NBR)))
all.cnvs.broad <- as.data.frame((all.cnvs.broad))
rownames(all.cnvs.broad) <- unique(gene.to.arm$ID)
colnames(all.cnvs.broad) <- unique(master.sheet$SAMPLE_ACCESSION_NBR)
all.cnvs.broad.gain <- all.cnvs.broad
all.cnvs.broad.loss <- all.cnvs.broad
## produce QC  plots
## first generate data structure to hold gene lists
gene.lists <- list()
for (i in 1:length(unique(gene.to.arm$ID))){
    arm <- unique(gene.to.arm$ID)[i]
    v1 <- gene.to.arm[gene.to.arm$ID == arm & gene.to.arm$v1 == "X", 1]
    v2 <- gene.to.arm[gene.to.arm$ID == arm & gene.to.arm$v2 == "X", 1]
    v3 <- gene.to.arm[gene.to.arm$ID == arm & gene.to.arm$v3 == "X", 1]
    gene.lists[[i]] <- list(v1, v2, v3)
    
}

## also produce not covered cnv data
not.covered.cnv.1 <- gene.to.arm$Gene[gene.to.arm$v1 == ""]
not.covered.cnv.2 <- gene.to.arm$Gene[gene.to.arm$v2 == ""]
not.covered.cnv.3 <- gene.to.arm$Gene[gene.to.arm$v3 == ""]
not.covered.cnv <- list(not.covered.cnv.1, not.covered.cnv.2, not.covered.cnv.3)
## generate list of ordered oncopanel versions
panel.version <- sapply(colnames(all.cnvs.broad), function(x){master.sheet[master.sheet$SAMPLE_ACCESSION_NBR == x, ]$PANEL_VERSION}, USE.NAMES = F)

for (i in 1:ncol(all.cnvs.broad)){
    ## for each sample
    
    for (j in 1:nrow(all.cnvs.broad)){
        ## loop through each chromosome arm
        all.cnvs.broad.loss[j, i] <- sum(gene.lists[[j]][[panel.version[i]]] %in% all.cnvs[all.cnvs$SAMPLE_ACCESSION_NBR == colnames(all.cnvs.broad)[i] &
                                                                                    all.cnvs$CNV_TYPE_CD %in% c("1DEL", "2DEL"), ]$GENE) / 
            length(gene.lists[[j]][[panel.version[i]]])
    }    
}

for (i in 1:nrow(all.cnvs.broad)){
    pdf(paste("../Analysis/CNV Histogram ", rownames(all.cnvs.broad)[i], ".pdf", sep = ""))
    hist(as.numeric(all.cnvs.broad[i, ]), main = paste("Histogram for chromosome ", rownames(all.cnvs.broad)[i]), breaks = (-10:10)/10)
    dev.off()
}

all.cnvs.broad <- all.cnvs.broad.gain - all.cnvs.broad.loss

write.table(all.cnvs.broad.gain, "../OncDRS data/all.cnvs.broad.gain.txt", sep = "\t", row.names = FALSE, quote = F)
write.table(all.cnvs.broad.loss, "../OncDRS data/all.cnvs.broad.loss.txt", sep = "\t", row.names = FALSE, quote = F)

all.cnvs.broad.gain <- read.delim("../OncDRS data/all.cnvs.broad.gain.txt", stringsAsFactors = F)
all.cnvs.broad.loss <- read.delim("../OncDRS data/all.cnvs.broad.loss.txt", stringsAsFactors = F)


## turns everything into binary
all.cnvs.broad[all.cnvs.broad < -0.499] <- -1
all.cnvs.broad[all.cnvs.broad > 0.499] <- 1
all.cnvs.broad[abs(all.cnvs.broad) != 1] <- 0

temp <- t(all.cnvs.broad[, colnames(all.cnvs.broad) %in% pathologies[["glioma"]]])
d <- dist(temp)
clust <- hclust(d)
order <- clust$labels[clust$order]
temp <- as.data.frame(temp)
temp$samples <- rownames(temp)


## add phenotype data

temp$subtype <- sapply(temp$samples, function(x){master.sheet[master.sheet$SAMPLE_ACCESSION_NBR == x, ]$Cancer_Type_Specific})
glioma.colors <- c("-1" = "blue",  "1" = "red", "AnaplasticAstro" = "green", "Angiocentric" = "orange", "Astro" = "green", 
                          "DiffuseAstro" = "green", "Glioblastoma" = "black", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green")

## generate plot
temp <- melt(temp, id = "samples")
temp$samples <- factor(temp$samples, levels = order)
temp$value[temp$value == 0] <- NA
mut <- ggplot(temp, aes(x=samples, y=variable, height=0.9, width=0.9)) + 
    geom_tile(aes(fill=value)) +
    scale_fill_manual(values = glioma.colors,  na.value = "white") +
    xlab("Subject") +
    ggtitle("Glioma unsupervised clustering") +
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

rowSums(temp)

