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
        good.rearrangements$even.typet[i] <- matches$mechanism
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
ggplot(data=plotting.matrix, aes(x=Sample, fill = complex.simple)) + geom_bar(width = .5, position = "dodge")+ scale_fill_grey() + 
labs(title= "Event type comparison", x = "Samples", y = "Rearrangement count") + theme(axis.text.x=element_blank()) +
scale_y_log10()


## create stacked bar with percentages instead of absolute counts
plotting.matrix$Sample <- as.character(plotting.matrix$Sample)
complex.num <- c()
simple.num <- c()
samples <- unique(plotting.matrix$Sample) 
for (i in 1:length(samples)){
    complex.num <- c(complex.num, sum(plotting.matrix[plotting.matrix$complex.simple == "complex.event", ]$Sample == samples[i]))
    simple.num <- c(simple.num, sum(plotting.matrix[plotting.matrix$complex.simple == "simple.event", ]$Sample == samples[i]))
}

complex.df <- data.frame(samples, complex.num / (complex.num + simple.num), simple.num / (complex.num + simple.num))

## meerdog event mechanism by grade
plotting.matrix$meerdog <- factor(plotting.matrix$meerdog, levels=c("NHEJ", "MMEJ", "MMBIR", "SSR"))
ggplot(data=plotting.matrix, aes(x=meerdog, fill=grade)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")


## look at possible fills and covariates
ggplot(data=plotting.matrix, aes(x=meerdog, fill=complex.simple)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")
ggplot(data=plotting.matrix, aes(x=event.type.simple, fill=complex.simple)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")
ggplot(data=plotting.matrix, aes(x=event.type.simple, fill=grade.and.outlier)) + geom_bar(width = .5)+ labs(title= "Event type comparison", x = "Event Type", y = "count")



fisher.test(table(plotting.matrix[plotting.matrix$Sample == "MEN0048G-P1", ]$event.type.simple == "translocation", 
                  plotting.matrix[plotting.matrix$Sample == "MEN0048G-P1", ]$complex.simple))

ggplot(data=plotting.matrix, aes(x=SPAN==-1, fill=complex.simple)) + geom_bar(width = .5) + scale_fill_grey() +
    labs(title= "Event type comparison", x = "Event Type", y = "count")

## meerdog event mechanism with grade separate
ggplot(data=plotting.matrix[plotting.matrix$grade == "LG",], aes(x=meerdog)) + geom_bar(width = .5)+ labs(title= "Low grade samples", x = "Event Type", y = "count")
ggplot(data=plotting.matrix[plotting.matrix$grade == "HG",], aes(x=meerdog)) + geom_bar(width = .5)+ labs(title= "High grade samples", x = "Event Type", y = "count")














## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#old snowman calcs
## Gets centromere position
centromere.positions <- read.delim("C:/Users/Noah/Onedrive/Work/Coding/R/dbs/centromere_positions.txt", stringsAsFactors = F, comment.char = "", header = T)
centromere.positions <- centromere.positions[order(centromere.positions$chrom), ]
mtrx <- matrix(NA, 24, 3)
colnames(mtrx) <- c("chr", "start", "end")



pretty <- x[is.na(x[[9]]), ]
pretty <- pretty[, c(2:8)]
pretty <- pretty[order(pretty$Gene1), ]
write.csv(pretty, file = "C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Plots/pretty.csv", row.names = F)



Gistic <- read.delim("C:/Users/Noah/Onedrive/Work/Meningioma/Analysis/GISTIC_Arm_Level.txt", stringsAsFactors = F, header = T)
helper <- function(input){
    if (input < 0){
        input <- -1
    }else{
        if(input >0){
            input <- 1
        }else{
            input <- 0
        }
    }
}

wgs <- Gistic[, c(1,2,5, 6, 15, 18, 19, 20, 21)]
wgs[,2] <- sapply(wgs[, 2], helper)
wgs[,3] <- sapply(wgs[, 3], helper)
wgs[,4] <- sapply(wgs[, 4], helper)
wgs[,5] <- sapply(wgs[, 5], helper)
wgs[,6] <- sapply(wgs[, 6], helper)
wgs[,7] <- sapply(wgs[, 7], helper)
wgs[,8] <- sapply(wgs[, 8], helper)
wgs[,9] <- sapply(wgs[, 9], helper)

full.names <- c()
current.names <- colnames(wgs)
for (i in 2:length(current.names)){
    full.names <- c(full.names, paste(current.names[i], "_Gistic", sep = ""))
    full.names <- c(full.names, paste(current.names[i], "_Rearrangments", sep = ""))
}

combined <- matrix(NA, 39, 17)
combined[, 1] <- wgs[, 1]
combined[, 2] <- wgs[, 2]
combined[, 4] <- wgs[, 3]
combined[, 6] <- wgs[, 4]
combined[, 8] <- wgs[, 5]
combined[, 10] <- wgs[, 6]
combined[, 12] <- wgs[, 7]
combined[, 14] <- wgs[, 8]
combined[, 16] <- wgs[, 9]

colnames(combined) <- c("Chromosome", full.names)
combined

## Takes chrom positions, translates them to binary instructions. 0 implies value is max, 1 implies val is a min
## Matrix used to split rearrangment calls between p and q arms using chr position
inputs <- matrix(NA, 39, 3)
rownames(inputs) <- Gistic[, 1]

for (i in 1:12){
    inputs[2 * i - 1, 1] <- mtrx[i, 2]
    inputs[2 * i - 1 ,2] <- 0
    inputs[2 * i - 1 ,3] <- i
    inputs[2 * i, 1] <- mtrx[i, 3]
    inputs[2 * i, 2] <- 1
    inputs[2 * i,3] <- i
}

inputs[25:27, 1:2] <- 1
inputs[25:27, 3] <- c(13, 14, 15)
inputs[38:39, 1:2] <- 1
inputs[38:39, 3] <- c(21, 22)

for (i in 14:18){
    inputs[2 * i, 1] <- mtrx[i + 2, 2]
    inputs[2 * i,2] <- 0
    inputs[2 * i,3] <- i + 2
    inputs[2 * i + 1, 1] <- mtrx[i + 2, 3]
    inputs[2 * i + 1, 2] <- 1
    inputs[2 * i + 1, 3] <- i + 2
}

## Loops through each individual in data set
for (i in 1:length(unique(data$individual))){
    ind <- FilterMaf(genes.all.unique, unique(data$individual)[i],"individual")
    
    ## loops through each chromosome per individual, assigning values to gistic
    for (j in 1:nrow(inputs)){
        chr <- FilterMaf(ind, inputs[j,3], "chr")
        count <- 0
        if (inputs[j, 2] == 0){
            count <- sum(chr$pos < inputs[j, 1])
        }else{
            count <- sum(chr$pos > inputs[j, 1])
        }
        combined[j, 2 * i + 1] <- count
    }
}

write.csv(combined, "C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Plots/master_chart.csv")


## Set up the graph
plot(1, type="n", xlab="Chromsome", ylab="count", xlim=c(0, 39), ylim=c(0, 10), xaxt = "n")
axis(1, seq(1:39), wgs[ ,1], las = 2)

idx <- seq(2, 16, 2)

## Loops through by chromosome

for (i in 1:39){
    ## Identify which indices match each condition, then save corresponding value into vector to get added to graph
    current <- combined[i, ]
    gains <- current[idx] == 1
    gains <- idx[gains]
    vals <- current[(gains + 1)]
    points(rep(i, length(vals)), vals, col = "red", pch = 16)
    
    ## Neutral
    neutral <- current[idx] == 0
    neutral <- idx[neutral]
    vals <- current[(neutral + 1)]
    points(rep(i, length(vals)), vals, col = "grey", pch = 16)
    
    ## Losses
    losses <- current[idx] == -1
    losses <- idx[losses]
    vals <- current[(losses + 1)]
    points(rep(i, length(vals)), vals, col = "blue", pch = 16)
}

## Plots

# takes max and min of possible centromere positions
for (i in 1:24){
    tmp <- unique(centromere.positions$chrom)[i]
    vals <- FilterMaf(centromere.positions, tmp, "chrom")
    mtrx[i,1] <- i
    mtrx[i,2] <- min(vals$chromStart)
    mtrx[i,3] <- max(vals$chromEnd)
}
cents <- data.frame(mtrx)

## check snowman vs dranger calls
disc.dranger <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/dRanger/old/dranger_WGS_disc_calls.txt", stringsAsFactors = F)
disc.dranger <- FilterMaf(disc.dranger, "MEN0049", "individual", F)
ph.dranger <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/dRanger/old/dranger_WGS_ph_calls.txt", stringsAsFactors = F)

disc.good.rearrangements <- good.rearrangements[good.rearrangements$Sample %in% unique(master.table[master.table$Cohort == "Discovery", ]$Pair.Name), ]
ph.good.rearrangements <- good.rearrangements[good.rearrangements$Sample %in% unique(master.table[master.table$Cohort == "PH", ]$Pair.Name), ]

for (i in 1:nrow(disc.good.rearrangements)){
    temp <- disc.good.rearrangements[i, "Sample"]
    temp <- substr(temp, 1, 7)
    disc.good.rearrangements[i, "Sample"] <- temp
}


discovery.comparison <- disc.dranger[, c(1,3,5,6,8)]
discovery.comparison[, 6] <- "dranger"
rownames(discovery.comparison) <- seq(nrow(discovery.comparison))

for (i in 1:nrow(disc.good.rearrangements)){
    matches <- discovery.comparison[discovery.comparison$individual == disc.good.rearrangements$Sample[i] & discovery.comparison$chr1 == disc.good.rearrangements$chr1[i]
                                    & discovery.comparison$pos1 %in% seq(from = disc.good.rearrangements$pos1[i] -3,length.out = 6), ]
    if (nrow(matches) == 0){
        discovery.comparison <- rbind(discovery.comparison, c(disc.good.rearrangements[i, 28], disc.good.rearrangements[i, 1], 
                                                              disc.good.rearrangements[i, 2], disc.good.rearrangements[i, 4], disc.good.rearrangements[i,5], "snowman"))
    }else{
        discovery.comparison[as.numeric(rownames(matches)[1]), 6] <- "both"
    }
}

discovery.comparison <- discovery.comparison[order(discovery.comparison$individual, as.numeric(discovery.comparison$chr1), as.numeric(discovery.comparison$pos1)),]
write.csv(discovery.comparison, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/discovery_cohort_comparison.csv", row.names = F, quote = F)






rearrangements[, 38:39] <- "Intergenic"
colnames(rearrangements)[38:39] <- c("gene1.or.cancer", "gene2.or.cancer")

## Creats columns for genes or cancer genes within 100kb
for (i in 1:nrow(rearrangements)){
    if (rearrangements$gene1[i] == "Intergenic"){
        if (rearrangements$cancer_genes100kb_1[i] != "none"){
            rearrangements$gene1.or.cancer[i] <- rearrangements$cancer_genes100kb_1[i]
        }
    }else{
        rearrangements$gene1.or.cancer[i] <- rearrangements$gene1[i]
    }
    if (rearrangements$gene2[i] == "Intergenic"){
        if (rearrangements$cancer_genes100kb_2[i] != "none"){
            rearrangements$gene2.or.cancer[i] <- rearrangements$cancer_genes100kb_2[i]
        }
    }else{
        rearrangements$gene2.or.cancer[i] <- rearrangements$gene2[i]
    }
}

rearrangements[, 17] <- NA
## Create cancer-related column for easy viewing
for (i in 1:nrow(rearrangements)){
    temp1 <- rearrangements$cancer_gene1[i]
    temp2 <- rearrangements$cancer_gene2[i]
    if (temp1 == "none"){
        if (temp2 == "none"){
            rearrangements[i, 17] <- "none"
        }else{
            rearrangements[i, 17] <- temp2
        }
    }else{
        if (temp2 == "none"){
            rearrangements[i, 17] <- temp1
        }else{
            rearrangements[i, 17] <- paste(temp1, temp2, sep = "_")
        }
    }
}
colnames(rearrangements)[17] <- "cancer.related"








## ------------------------------------------------------------------------------------------------------------------------------------------
## dRanger text




## Gets start and end of each gene
hg.reference <- read.delim("C:/Users/Noah/Onedrive/Work/Coding/R/dbs/all_genes_hg38.txt", stringsAsFactors = F, comment.char = "", header = T)
reference <- hg.reference[, c(3, 5,6,13)]

## Gets centromere position
centromere.positions <- read.delim("C:/Users/Noah/Onedrive/Work/Coding/R/dbs/centromere_positions.txt", stringsAsFactors = F, comment.char = "", header = T)
centromere.positions <- centromere.positions[order(centromere.positions$chrom), ]
mtrx <- matrix(NA, 24, 3)
colnames(mtrx) <- c("chr", "start", "end")

# takes max and min of possible centromere positions
for (i in 1:24){
    tmp <- unique(centromere.positions$chrom)[i]
    vals <- FilterMaf(centromere.positions, tmp, "chrom")
    mtrx[i,1] <- i
    mtrx[i,2] <- min(vals$chromStart)
    mtrx[i,3] <- max(vals$chromEnd)
}
cents <- data.frame(mtrx)

## read in snowman data
snowman.folder.rear <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/rearrangements")
snowman.rearrangements <- NULL
for (i in list.files(snowman.folder.rear)){
    temp <- read.csv(paste(snowman.folder.rear, list.files(snowman.folder.rear)[i], sep = "/"),
                     stringsAsFactors = F)  
    snowman.rearrangements <- rbind(snowman.rearrangements, temp)    
}

## Read in data, remove MEN049 which is not a meningioma
data <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/dRanger.txt", stringsAsFactors = F)
data <- FilterMaf(data, "MEN0049", "individual", F)

## gets all genes
genes1 <-data[,c(1,3,5, 24)]
genes2 <- data[,c(1,6, 8, 26)]
colnames(genes2)[2:4] <- c("chr1", "pos1", "gene1")
genes.all <- rbind(genes1, genes2)
colnames(genes.all) <- c("individual", "chr", "pos", "gene")

## Plot rearrangments per gene
genes.all.4 <- ReccurentMaf(genes.all,"gene", 3)
tbl <- table(genes.all.4$gene)
tbl <- tbl[order(tbl, decreasing = T)]
par(mar = c(8, 4,4,2))
barplot(tbl, las = 2, main = "# of rearrangments per gene, minimum 4")

## Plot unique rearrangments per gene
genes.all.unique <- PerSampleMaf(genes.all, "gene", identifier.column =  "individual")
genes.all.unique.2 <- ReccurentMaf(genes.all.unique, "gene")
tbl <- table(genes.all.unique.2$gene)
tbl <- tbl[order(tbl, decreasing = T)]
par(mar = c(8, 4,4,2))
barplot(tbl, las = 2, main = "# of unique rearrangments per gene, min 2")

## Plot unique rearrangments per reasonable distance
## extract location, and distance from actual gene
data$site1.simple 
x <- sapply(data$site1[1:10], function(x) strsplit(x, " ")[[1]][c(1,2)], USE.NAMES = F)
grep("[0-9]kb", "3kb")
strsplit("3kb", "[0-9]")


## Find rearrangmentes within region bordering genes
regions <- c()
names <- c()
interval <- 2
small <- genes.all[!duplicated(genes.all$gene), ]

## Loops through gene list
for (i in 1:nrow(small)){
    cur <- small$gene[i]
    chr <- small$chr[i]
    squad <- paste(cur, "_posse", sep = "")
    min <- 0
    
    ## gets all matching tx start/end positions for gene
    short.list <- reference[reference$name2 == cur, ]
    if (nrow(short.list) > 0){
        
        ## sets min and max based on values
        if (min(short.list$txStart) > interval){
            min <- min(short.list$txStart) - interval
        }
        max <- max(short.list$txEnd) + interval
        
        squad_size <- sum(genes.all.unique$pos > min & genes.all.unique$pos < max & genes.all.unique$chr == chr)
        regions <- c(regions, squad_size)
        names <- c(names, squad)
    }
}
names(regions) <- names
regions <- regions[order(regions, decreasing = T)]
par(mar = c(10, 4,4,2))
barplot(regions[1:50], las = 2, main = "# of rearrangments per region")

# Find rearrangmentes within given arbitrary segments
regions.agnostic <- c()
names.agnostic <- c()
chr.agnostic <- c()
interval <- 200000

## Find max length of each chromosome
the.buck.stops.here <- c()
chromosome.names <- sort(unique(genes.all.unique$chr))
for (i in 1:length(chromosome.names)){
    temp <- FilterMaf(genes.all.unique, chromosome.names[i], "chr")
    the.buck.stops.here <- c(the.buck.stops.here, max(temp$pos))
}

## Loops through chromosomes, breaking up into even regions
for (i in 1:22){
    universe <- FilterMaf(genes.all.unique, i, "chr")
    
    for (j in 1:(the.buck.stops.here[i] / interval)){
        
        ## sets min and max based on values
        min <- (j - 1) * interval
        max <- min + interval
        squad_size <- sum(genes.all.unique$pos > min & genes.all.unique$pos < max)
        
        if (squad_size > 0){
            regions.agnostic <- c(regions.agnostic, squad_size)
            chr.agnostic <- c(chr.agnostic, i)
        }
    }
    
}

plot(regions.agnostic)

regions.agnostic <- regions.agnostic[order(regions.agnostic, decreasing = T)]
barplot(regions.agnostic[1:50], las = 2, main = "# of rearrangments per region")


## Gets all recurrent rearrangments between pairs
matches <- NULL

for (i in 1:nrow(data)){
    # Gets first gene-gene pair
    current.gene <- data$gene1[i]
    gene.chr <- data$chr1[i]
    gene.pos <- data$pos1[i]
    current.target <- data$gene2[i]
    target.chr <- data$chr2[i]
    target.pos <- data$pos2[i]
    friends <- FilterMaf(data[-i,], current.gene, "gene1")
    
    #checks to see if any other instances of same starter gene have same
    if (nrow(friends) > 0){
        idx <- friends$gene2 %in% current.target
        close.friends <- friends[idx,]
        if (nrow(close.friends) > 0){
            temp <- c(rownames(data)[i], current.gene, gene.chr, gene.pos, current.target, 
                    target.chr, target.pos, data$individual[i])
            matches <- rbind(matches, temp)
            
        }
    
    } 
}

matches.table <- matrix(matches, length(matches)/8, 8)
pairwise.matches <- data.frame(matches.table, stringsAsFactors = FALSE)
colnames(pairwise.matches) <- c("Original Index", "Gene1", "Chr1", "Pos1", "Gene2", "Chr2", "pos2", "sample")
x <- pairwise.matches
## gets rid of duplicates from same sample

## Sets current row info
for (i in 1:nrow(x)){
    current.gene <- x[, 2][i]
    current.target <- x[, 5][i]
    current.individual <- x[, 8][i]
    
    # Finds all samples with same start
    friends <- FilterMaf(x[-i,], current.gene, 2)
    
    if (nrow(friends) > 0){
        
        #keeps those with same target
        idx <- friends[, 5] %in% current.target
        close.friends <- friends[idx,]
        
        # marks duplicates if same individual
        if (nrow(close.friends) > 0){
            for (j in 1:nrow(close.friends)){
                if (close.friends[, 8][j] == current.individual){
                    x[as.numeric(rownames(close.friends)[j]), 9] <- 1
                }
                  
            }
        
        }
    
    
    }

} 


x[18, 9] <- NA

pretty <- x[is.na(x[[9]]), ]
pretty <- pretty[, c(2:8)]
pretty <- pretty[order(pretty$Gene1), ]
write.csv(pretty, file = "C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Plots/pretty.csv", row.names = F)

Gistic <- read.delim("C:/Users/Noah/Onedrive/Work/Meningioma/Analysis/GISTIC_Arm_Level.txt", stringsAsFactors = F, header = T)
helper <- function(input){
    if (input < 0){
        input <- -1
    }else{
        if(input >0){
            input <- 1
        }else{
            input <- 0
        }
    }
}

wgs <- Gistic[, c(1,2,5, 6, 15, 18, 19, 20, 21)]
wgs[,2] <- sapply(wgs[, 2], helper)
wgs[,3] <- sapply(wgs[, 3], helper)
wgs[,4] <- sapply(wgs[, 4], helper)
wgs[,5] <- sapply(wgs[, 5], helper)
wgs[,6] <- sapply(wgs[, 6], helper)
wgs[,7] <- sapply(wgs[, 7], helper)
wgs[,8] <- sapply(wgs[, 8], helper)
wgs[,9] <- sapply(wgs[, 9], helper)

full.names <- c()
current.names <- colnames(wgs)
for (i in 2:length(current.names)){
    full.names <- c(full.names, paste(current.names[i], "_Gistic", sep = ""))
    full.names <- c(full.names, paste(current.names[i], "_Rearrangments", sep = ""))
}

combined <- matrix(NA, 39, 17)
combined[, 1] <- wgs[, 1]
combined[, 2] <- wgs[, 2]
combined[, 4] <- wgs[, 3]
combined[, 6] <- wgs[, 4]
combined[, 8] <- wgs[, 5]
combined[, 10] <- wgs[, 6]
combined[, 12] <- wgs[, 7]
combined[, 14] <- wgs[, 8]
combined[, 16] <- wgs[, 9]

colnames(combined) <- c("Chromosome", full.names)
combined

## Takes chrom positions, translates them to binary instructions. 0 implies value is max, 1 implies val is a min
## Matrix used to split rearrangment calls between p and q arms using chr position
inputs <- matrix(NA, 39, 3)
rownames(inputs) <- Gistic[, 1]

for (i in 1:12){
    inputs[2 * i - 1, 1] <- mtrx[i, 2]
    inputs[2 * i - 1 ,2] <- 0
    inputs[2 * i - 1 ,3] <- i
    inputs[2 * i, 1] <- mtrx[i, 3]
    inputs[2 * i, 2] <- 1
    inputs[2 * i,3] <- i
}

inputs[25:27, 1:2] <- 1
inputs[25:27, 3] <- c(13, 14, 15)
inputs[38:39, 1:2] <- 1
inputs[38:39, 3] <- c(21, 22)

for (i in 14:18){
    inputs[2 * i, 1] <- mtrx[i + 2, 2]
    inputs[2 * i,2] <- 0
    inputs[2 * i,3] <- i + 2
    inputs[2 * i + 1, 1] <- mtrx[i + 2, 3]
    inputs[2 * i + 1, 2] <- 1
    inputs[2 * i + 1, 3] <- i + 2
}

## Loops through each individual in data set
for (i in 1:length(unique(data$individual))){
    ind <- FilterMaf(genes.all.unique, unique(data$individual)[i],"individual")
    
    ## loops through each chromosome per individual, assigning values to gistic
    for (j in 1:nrow(inputs)){
        chr <- FilterMaf(ind, inputs[j,3], "chr")
        count <- 0
        if (inputs[j, 2] == 0){
            count <- sum(chr$pos < inputs[j, 1])
        }else{
            count <- sum(chr$pos > inputs[j, 1])
        }
        combined[j, 2 * i + 1] <- count
    }
}

write.csv(combined, "C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Plots/master_chart.csv")


## Set up the graph
plot(1, type="n", xlab="Chromsome", ylab="count", xlim=c(0, 39), ylim=c(0, 10), xaxt = "n")
axis(1, seq(1:39), wgs[ ,1], las = 2)

idx <- seq(2, 16, 2)

## Loops through by chromosome

for (i in 1:39){
    ## Identify which indices match each condition, then save corresponding value into vector to get added to graph
    current <- combined[i, ]
    gains <- current[idx] == 1
    gains <- idx[gains]
    vals <- current[(gains + 1)]
    points(rep(i, length(vals)), vals, col = "red", pch = 16)
    
    ## Neutral
    neutral <- current[idx] == 0
    neutral <- idx[neutral]
    vals <- current[(neutral + 1)]
    points(rep(i, length(vals)), vals, col = "grey", pch = 16)
    
    ## Losses
    losses <- current[idx] == -1
    losses <- idx[losses]
    vals <- current[(losses + 1)]
    points(rep(i, length(vals)), vals, col = "blue", pch = 16)
}


## Length distribution analysis
hist(data$span, 30, xlab = "rearrangmenet distance", main = "Length Distribution")

## inter vs intra-chromosomal 


