source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")


total.rearrangements <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/Rearrangements/v117.events.tsv", stringsAsFactors = F, header = T)
annotated <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/Rearrangements/events.under500cut.txt", stringsAsFactors = F)

## Keeps only those that pass QC filter
## If span is greater than 1000 bases, allows for 5% allelic fraction and 2 reads
## if span is less than 1000 bases, requires 10% af and at least 5 reads
good.rearrangements <- total.rearrangements[(total.rearrangements$SPAN > 1000 | total.rearrangements$SPAN == -1) & total.rearrangements[, "TUMALT"] > 1 & 
                                                (total.rearrangements[, "TUMALT"] / total.rearrangements[, "TUMCOV"] > .05 )
                                            | total.rearrangements[, "TUMCOV"] > 5 & (total.rearrangements[,"TUMALT"] / total.rearrangements[, "TUMCOV"] > .10 ) , ]

## combine annotation of event type from poorly formatted sheet
good.rearrangements$event.type <- NA
good.rearrangements$complex <- NA
for (i in 1:nrow(good.rearrangements)){
    chr1 <- good.rearrangements$chr1[i]
    chr2 <- good.rearrangements$chr2[i]
    pos1 <- good.rearrangements$pos1[i]
    pos2 <- good.rearrangements$pos2[i]
    sample <- good.rearrangements$Sample[i]
    matches <- annotated[annotated$chr1 == chr1 & annotated$pos1 == pos1 & annotated$chr2 == chr2 & annotated$pos2 == pos2 & annotated$Sample == sample,]
    if (nrow(matches) == 1){
        good.rearrangements$event.type[i] <- matches$mechanism
        good.rearrangements$complex[i] <- matches$complex.event
    }else if (nrow(matches) > 1){
        stop("error")
    }else{
        ## do nothing
    }
}

good.rearrangements$event.type.simple <- good.rearrangements$event.type
good.rearrangements$event.type.simple[good.rearrangements$SPAN == -1] <- "translocation"
good.rearrangements$complex.simple <- good.rearrangements$complex
good.rearrangements$complex.simple[good.rearrangements$complex.simple != ""] <- "complex.event"
good.rearrangements$complex.simple[good.rearrangements$complex.simple == ""] <- "simple.event"
rearrangements <- good.rearrangements

rearrangements <- meerdog(rearrangements)

write.table(rearrangements, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/rearrangements_final.txt",
            quote = F, row.names = F, sep = "\t")

## removes NAs, look for multiple hits
rearrangements[is.na(rearrangements$gene1), "gene1"] <- ""
rearrangements[is.na(rearrangements$gene2), "gene2"] <- ""
rearrangements[is.na(rearrangements$gene1_100kb), "gene1_100kb"] <- ""
rearrangements[is.na(rearrangements$gene2_100kb), "gene2_100kb"] <- ""

## creates column that has gene name if within gene, otherwise genes within 100kb
rearrangements$gene1_or_100kb <- rearrangements$gene1
gene1.idx <- rearrangements$gene1 == "" & rearrangements$gene1_100kb != ""
rearrangements$gene1_or_100kb[gene1.idx] <- rearrangements$gene1_100kb[gene1.idx]

rearrangements$gene2_or_100kb <- rearrangements$gene2
gene2.idx <- rearrangements$gene2 == "" & rearrangements$gene2_100kb != ""
rearrangements$gene2_or_100kb[gene2.idx] <- rearrangements$gene2_100kb[gene2.idx]

rearrangements.duplicates <- rearrangements
rearrangements.duplicates$remove <- 0

## takes any sample with multiple genes within range in gene1, converts to duplicate row with each gene as its own row
for (i in 1:nrow(rearrangements)){
    genes <- strsplit(rearrangements$gene1_or_100kb[i], ",")[[1]]
    if (length(genes) > 1){
        row <- rearrangements.duplicates[i, ]
        for (j in 1:length(genes)){
            row$gene1_or_100kb <- genes[j]
            rearrangements.duplicates <- rbind(rearrangements.duplicates, row)
        }
        rearrangements.duplicates$remove[i] <- 1
    }
}

rearrangements.duplicates <- rearrangements.duplicates[rearrangements.duplicates$remove != 1, ]
## same thing for gene 2

rearrangements.duplicates.duplicates <- rearrangements.duplicates
rearrangements.duplicates.duplicates$remove <- 0

## takes any sample with multiple genes within range in gene2, converts to duplicate row with each gene as its own row
for (i in 1:nrow(rearrangements.duplicates)){
    genes <- strsplit(rearrangements.duplicates$gene2_or_100kb[i], ",")[[1]]
    if (length(genes) > 1){
        row <- rearrangements.duplicates.duplicates[i, ]
        for (j in 1:length(genes)){
            row$gene2_or_100kb <- genes[j]
            rearrangements.duplicates.duplicates <- rbind(rearrangements.duplicates.duplicates, row)
        }
        rearrangements.duplicates.duplicates$remove[i] <- 1
    }
}

rearrangements.duplicates.duplicates <- rearrangements.duplicates.duplicates[rearrangements.duplicates.duplicates$remove != 1, ]

## gets all genes, using gene or  gene 100k for gene category if intergenic
genes1 <- rearrangements.duplicates.duplicates[,c("sample", "chr1", "pos1", "gene1_or_100kb", "chr2", "pos2", "gene2_or_100kb")]
genes2 <- rearrangements.duplicates.duplicates[,c("sample", "chr2", "pos2", "gene2_or_100kb", "chr1", "pos1", "gene1_or_100kb")]
colnames(genes2) <- c("sample", "chr", "pos", "gene", "partner.chr", "partner.pos", "partner.gene")
colnames(genes1) <- c("sample", "chr", "pos", "gene", "partner.chr", "partner.pos", "partner.gene")
genes.all <- rbind(genes1, genes2)

## Plot recurrence-unique rearrangements per gene
genes.modified <- genes.all
genes.modified <- FilterMaf(genes.modified, "", "gene", F)
genes.modified <- PerSampleMaf(genes.modified, "gene", "sample")
genes.modified.2 <- ReccurentMaf(genes.modified, "gene")
genes.modified.2 <- genes.modified.2[order(genes.modified.2$sample), ]
PlotMaf(genes.modified.2, "gene")

genes.modified.2 <- genes.modified.2[order(genes.modified.2$gene), ]
##Investigate multiple hits
write.csv(genes.modified.2, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/multiple_hits.csv", row.names = F)


## Gets all recurrent rearrangments between pairs of genes. Replaces intergenic regions with nearby cancer genes
gene.rearrangements <- FilterMaf(rearrangements, "", "gene1_or_100kb", F)
gene.rearrangements <- FilterMaf(gene.rearrangements, "", "gene2_or_100kb", F)
matches <- NULL
for (i in 1:nrow(gene.rearrangements)){
    # Gets first gene-gene pair
    current.gene <- gene.rearrangements$gene1_or_100kb[i]
    gene.chr <- gene.rearrangements$chr1[i]
    gene.pos <- gene.rearrangements$pos1[i]
    current.target <- gene.rearrangements$gene2_or_100kb[i]
    target.chr <- gene.rearrangements$chr2[i]
    target.pos <- gene.rearrangements$pos2[i]
    friends <- FilterMaf(gene.rearrangements[-i,], current.gene, "gene1_or_100kb")
    
    #checks to see if any other instances of same starter gene have same target
    if (nrow(friends) > 0){
        idx <- friends$gene2_or_100kb %in% current.target
        close.friends <- friends[idx,]
        if (nrow(close.friends) > 0){
            temp <- c(rownames(gene.rearrangements)[i], current.gene, gene.chr, gene.pos, current.target, 
                      target.chr, target.pos, gene.rearrangements$sample[i])
            matches <- c(matches, temp)
            
        }
        
    } 
}

matches.table <-matrix(matches, length(matches)/8, 8, byrow = T)
pairwise.matches <- data.frame(matches.table, stringsAsFactors = FALSE)
colnames(pairwise.matches) <- c("Original Index", "Gene1", "Chr1", "Pos1", "Gene2", "Chr2", "pos2", "sample")
x <- pairwise.matches

## gets rid of duplicates from same sample
for (i in 1:nrow(x)){
    ## Sets current row info
    current.gene <- x[, 2][i]
    current.target <- x[, 5][i]
    current.individual <- x[, 8][i]
    
    # Finds all samples with same starting gene
    friends <- FilterMaf(x[-i,], current.gene, 2)
    
    if (nrow(friends) > 0){
        
        #keeps those with same target
        idx <- friends[, 5] %in% current.target
        close.friends <- friends[idx,]
        
        # marks duplicates if all examples come from same individual, not if more than 1 sample
        if (nrow(close.friends) > 0){
            ## If only 1 total hit, check if same as original
            if (nrow(close.friends) == 1){
                if (close.friends[, 8] == current.individual){
                    x[as.numeric(rownames(close.friends)[1]), 9] <- 1
                }
                
                ## check if at least 2 unique hits        
            }else if (length(unique(close.friends[, 8])) > 1){
                ## Do nothing
                
                ## check and see if non-unique current samples are different from original    
            }else if (close.friends[1, 8] !=current.individual){
                # Do nothing
                
                #Otherwise, must be redundant hits    
            }else{
                for (j in 1:nrow(close.friends)){
                    x[as.numeric(rownames(close.friends)[j]), 9] <- 1
                }
                
            }
            
        }
        
    }
    
} 


## Figures
#remove NAs
plotting.matrix <- rearrangements
plotting.matrix <- plotting.matrix[plotting.matrix$event.type != "", ]
plotting.matrix <- plotting.matrix[!is.na(plotting.matrix$event.type), ]

## reclassify based on grade, as well as hyper-rearranged sample
lowgrade.num <- sum(plotting.matrix$Sample %in% unique(plotting.matrix$Sample)[1:9])
grade <- c(rep("LG", lowgrade.num), rep("HG", nrow(plotting.matrix) - lowgrade.num))
outlier <- plotting.matrix$Sample == "MEN0048G-P1"
grade.and.outlier <- grade
grade.and.outlier[outlier] <- "MEN0048"
plotting.matrix <- cbind(plotting.matrix, grade, grade.and.outlier)

low.grade <- table(plotting.matrix[plotting.matrix$grade == "LG", ]$Sample)
high.grade <- table(plotting.matrix[plotting.matrix$grade == "HG", ]$Sample)

mean(low.grade)
mean(high.grade)
median(low.grade)
median(high.grade)
high.grade.cleaned <- high.grade[-4]

wilcox.test(low.grade, high.grade)

## for fixed analysis, taking percentage of translocations separately
translocation.percent <- sum(plotting.matrix$event.type.simple == "translocation") / nrow(plotting.matrix)
inversion.percent <- sum(plotting.matrix$event.type.simple == "inversion") / (nrow(plotting.matrix)) - sum(plotting.matrix$event.type.simple == "translocation"))
deletion.percent <- sum(plotting.matrix$event.type.simple == "simple_deletion") / (nrow(plotting.matrix)) - sum(plotting.matrix$event.type.simple == "translocation"))
duplication.percent <- sum(plotting.matrix$event.type.simple == "tandeum_dup") / (nrow(plotting.matrix)) - sum(plotting.matrix$event.type.simple == "translocation"))
event.type.df <- data.frame(c("Translocation", "Deletion", "Inversion", "Duplication"), c(translocation.percent, inversion.percent, deletion.percent, duplication.percent))
colnames(event.type.df) <- c("event.classification", "vals")
ggplot(data=event.type.df, aes(x=event.classification, y=vals)) + geom_bar(stat="identity") + labs(title = "Translocations vs Intrachromosomal Rearrangements", x = "Event Type", y = "Count")


## non-fixed analysis
ggplot(plotting.matrix[plotting.matrix$event.type.simple != "balanced_translocation_intra", ], 
       aes(x="", fill=event.type.simple))+ geom_bar(width = 1) + coord_polar("y") + theme(axis.text.x=element_blank())


## event type by grade
plotting.matrix$event.type.simple <- factor(plotting.matrix$event.type.simple, levels = c("translocation", "simple_deletion", "inversion", "tandeum_dup", 
                                                                                          "balanced_translocation_intra"))
ggplot(data=plotting.matrix[plotting.matrix$event.type.simple != "balanced_translocation_intra", ], 
       aes(x=event.type.simple)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")



## complex events percentages
x <- table(plotting.matrix$complex.simple)
x[1] / (x[1] + x[2])
complex.names.idx <- unique(plotting.matrix$sample) %in% plotting.matrix[plotting.matrix$complex.simple == "complex.event", ]$sample

mean(table(plotting.matrix$sample)[complex.names.idx])
mean(table(plotting.matrix$sample)[!complex.names.idx])
median(table(plotting.matrix$sample)[complex.names.idx])
median(table(plotting.matrix$sample)[!complex.names.idx])
wilcox.test(table(plotting.matrix$sample)[complex.names.idx], table(plotting.matrix$sample)[!complex.names.idx])

mean(table(plotting.matrix[plotting.matrix$complex.simple == "complex.event", ]$sample) / 
         table(plotting.matrix[plotting.matrix$sample %in% unique(plotting.matrix$sample)[complex.names.idx], ]$sample))

## Comparison of number of complex vs noncomplex events
ggplot(data=plotting.matrix, aes(x=complex.simple)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count") + scale_fill_grey()

## plot of events per sample that are complex
plotting.matrix$Sample <- factor(plotting.matrix$Sample, levels=names(sort(table(plotting.matrix$Sample), decreasing = T)))
ggplot(data=plotting.matrix[plotting.matrix$Sample %in% unique(plotting.matrix$Sample)[complex.names.idx], ], 
       aes(x=Sample, fill = complex.simple)) + geom_bar(width = .5)+ scale_fill_grey() + labs(title= "Event type comparison", x = "Event Type", y = "count")



## log scale adjacent graph sorted by presence of complex event or not
plotting.matrix$Sample <- as.character(plotting.matrix$Sample)
low.grades <- sort(unique(plotting.matrix$Sample))[c(1:9, 17)]
high.grades <- sort(unique(plotting.matrix$Sample))[!(unique(plotting.matrix$Sample) %in% low.grades)]

hg.complex.order <- names(sort(table(plotting.matrix[plotting.matrix$complex.simple == "complex.event" & plotting.matrix$Sample %in% high.grades, ]$Sample), decreasing = T))
hg.simple.order <- names(sort(table(plotting.matrix[!(plotting.matrix$Sample %in% hg.complex.order) & plotting.matrix$Sample %in% high.grades, ]$Sample), decreasing = T))

lg.complex.order <- names(sort(table(plotting.matrix[plotting.matrix$complex.simple == "complex.event" & plotting.matrix$Sample %in% low.grades, ]$Sample), decreasing = T))
lg.simple.order <- names(sort(table(plotting.matrix[!(plotting.matrix$Sample %in% lg.complex.order) & plotting.matrix$Sample %in% low.grades, ]$Sample), decreasing = T))


plotting.matrix$Sample <- factor(plotting.matrix$Sample, levels = c(hg.complex.order, hg.simple.order, lg.complex.order, lg.simple.order))
ggplot(data=plotting.matrix, aes(x=Sample, fill = complex.simple)) + geom_bar(width = .5, position = "dodge")+  
    labs(title= "Event type comparison", x = "Samples", y = "Rearrangement count")  + rameen_theme
scale_y_log10()

## percentage instead of absolute count
complex.num <- c()
simple.num <- c()
samples <- unique(plotting.matrix$Sample) 
for (i in 1:length(samples)){
    complex.num <- c(complex.num, sum(plotting.matrix[plotting.matrix$complex.simple == "complex.event", ]$Sample == samples[i]))
    simple.num <- c(simple.num, sum(plotting.matrix[plotting.matrix$complex.simple == "simple.event", ]$Sample == samples[i]))
}

complex.df <- data.frame(c(samples, samples), c(complex.num / (complex.num + simple.num), simple.num / (complex.num + simple.num)), 
                         c(rep("complex",length(complex.num)), rep("simple", length(simple.num))))
colnames(complex.df) <- c("sample", "percent", "event")
complex.df$sample <- factor(complex.df$sample, levels = c(complex.order, simple.order))


ggplot(data=complex.df, aes(x=sample, fill = event, y = percent)) + geom_bar(width = .5, stat= "identity")+  
    labs(title= "Event type comparison", x = "Samples", y = "Rearrangement count")


## meerdog event mechanism by grade
plotting.matrix$meerdog <- factor(plotting.matrix$meerdog, levels=c("NHEJ", "MMEJ", "MMBIR", "SSR"))
ggplot(data=plotting.matrix, aes(x=meerdog, fill=grade)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")


## look at possible fills and covariates
ggplot(data=plotting.matrix, aes(x=meerdog, fill=complex.simple)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")
ggplot(data=plotting.matrix, aes(x=event.type.simple, fill=complex.simple)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")
ggplot(data=plotting.matrix, aes(x=event.type.simple, fill=grade.and.outlier)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")



fisher.test(table(plotting.matrix[plotting.matrix$Sample == "MEN0048G-P1", ]$event.type.simple == "translocation", 
                  plotting.matrix[plotting.matrix$Sample == "MEN0048G-P1", ]$complex.simple))

ggplot(data=plotting.matrix, aes(fill=SPAN==-1, x=complex.simple)) + geom_bar(width = .5) + scale_fill_grey() +
    labs(title= "Event type comparison", x = "Event Type", y = "count") + rameen_theme

## meerdog event mechanism with grade separate
ggplot(data=plotting.matrix[plotting.matrix$grade == "LG",], aes(x=meerdog)) + geom_bar(width = .5)+ labs(title= "Low grade samples", x = "Event Type", y = "count")
ggplot(data=plotting.matrix[plotting.matrix$grade == "HG",], aes(x=meerdog)) + geom_bar(width = .5)+ labs(title= "High grade samples", x = "Event Type", y = "count")


## compare incidence we find to incidence across cancer types.

pancan.rearrangements <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Meerkat/pancan_rearrangement_classification.txt", 
                                    stringsAsFactors = F)
pancan.names <- unique(pancan.rearrangements$patient.ID)

pancan.tumor.type <- c(rep("breast", sum(pancan.rearrangements$patient.ID %in% pancan.names[1:35])),
                       rep("crc", sum(pancan.rearrangements$patient.ID %in% pancan.names[36:49])),
                       rep("gbm", sum(pancan.rearrangements$patient.ID %in% pancan.names[50:65])),
                       rep("kirc", sum(pancan.rearrangements$patient.ID %in% pancan.names[66:68])),
                       rep("lusc", sum(pancan.rearrangements$patient.ID %in% pancan.names[69:86])),
                       rep("hcc", sum(pancan.rearrangements$patient.ID %in% pancan.names[87:105])),
                       rep("mm", sum(pancan.rearrangements$patient.ID %in% pancan.names[106:112])),
                       rep("ovo", sum(pancan.rearrangements$patient.ID %in% pancan.names[113:121])),
                       rep("pro", sum(pancan.rearrangements$patient.ID %in% pancan.names[122:128])),
                       rep("ucec", sum(pancan.rearrangements$patient.ID %in% pancan.names[129:138])), "")
pancan.rearrangements <- cbind(pancan.rearrangements, pancan.tumor.type)
table(pancan.rearrangements$mechanism)
cleaned.pancan <- pancan.rearrangements[pancan.rearrangements$mechanism %in% c("alt-EJ", "FoSTeS", "NHEJ", "VNTR"),]
cleaned.pancan <- cleaned.pancan[, c(1,3,4,26)]
cleaned.meningioma <- plotting.matrix[, c("Sample", "event.type.simple", "meerdog")]
cleaned.meningioma <- cbind(cleaned.meningioma, rep("meningioma", nrow(cleaned.meningioma)))
colnames(cleaned.meningioma) <- colnames(cleaned.pancan)

cleaned.pancan[cleaned.pancan$mechanism == "alt-EJ", ]$mechanism <- "MMEJ"
cleaned.pancan[cleaned.pancan$mechanism == "FoSTeS", ]$mechanism <- "MMBIR"
cleaned.combined <- rbind(cleaned.pancan, cleaned.meningioma)

ggplot(data=cleaned.combined, aes(x=pancan.tumor.type, fill=mechanism)) + geom_bar(position = "dodge") 
scale_y_log10()


## plot as percent
cleaned.combined.raw <- table(cleaned.combined$mechanism, cleaned.combined$pancan.tumor.type)

cleaned.percent.mtrx <- as.matrix(t(cleaned.combined.raw))
rownames(cleaned.percent.mtrx)

cleaned.defactor.df <- cbind(as.numeric(cleaned.percent.mtrx[,1]), 
                             as.numeric(cleaned.percent.mtrx[,2]), as.numeric(cleaned.percent.mtrx[,3]))


cleaned.defactor.df <- cleaned.defactor.df[-1, ]
cleaned.defactor.df <- cbind(cleaned.defactor.df, as.numeric(rowSums(cleaned.defactor.df)))
cleaned.defactor.df[,1] <- cleaned.defactor.df[,1] / cleaned.defactor.df[,4]
cleaned.defactor.df[,2] <- cleaned.defactor.df[,2] / cleaned.defactor.df[,4]
cleaned.defactor.df[,3] <- cleaned.defactor.df[,3] / cleaned.defactor.df[,4]
cleaned.defactor.df <- cleaned.defactor.df[, -4]
rownames(cleaned.defactor.df) <- colnames(cleaned.combined.raw)[-1]
colnames(cleaned.defactor.df) <- rownames(cleaned.combined.raw)[-(4:5)]

names <- c()
vals <- c()
type <- rep(c("MMBIR", "MMEJ", "NHEJ"), 11)
for (i in 1:(nrow(cleaned.defactor.df))){
    names <- c(names, rep(rownames(cleaned.defactor.df)[i], 3))
    vals <- c(vals, as.numeric(cleaned.defactor.df[i, ]))
}

small.df <- data.frame(names, vals, type, stringsAsFactors = F)
small.df$names <- factor(small.df$names, levels = unique(small.df$names)[order(small.df[small.df$type == "MMEJ", ]$vals, decreasing = T)])
small.df$type <- factor(small.df$type, levels = c("MMEJ", "NHEJ", "MMBIR"))
color.theme <- c("#edf8b1", "#7fcdbb", "#2c7fb8")
color.theme2 <- c("#ece2f0", "#a6bddb", "#1c9099")
ggplot(data=small.df, aes(x=names, y = vals, fill = type)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = color.theme2) +
    rameen_theme
