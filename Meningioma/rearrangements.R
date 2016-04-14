source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

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
snowman.folder.rear <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/Discovery/rearrangements/")
snowman.rearrangements <- NULL
for (i in 1:length(list.files(snowman.folder.rear))){
    temp <- read.csv(paste(snowman.folder.rear, list.files(snowman.folder.rear)[i], sep = "/"),
                     stringsAsFactors = F)
    temp[, 28] <-  strsplit(list.files(snowman.folder.rear)[i], ".csv")[[1]]
    colnames(temp)[28] <- "Sample"
    temp <- temp[, -29]
    snowman.rearrangements <- rbind(snowman.rearrangements, temp)    
}

## Keep rearrangements only with likely somatic score

snowman.rearrangements <- snowman.rearrangements[snowman.rearrangements$somatic_lod > 4, ]


##Keep only relevant columns
rearrangements <- snowman.rearrangements[, c(28, 1,2,29,31,33,35, 4,5,30, 32,34,36)]
rearrangements[, 14:15] <- "Intergenic"
colnames(rearrangements)[14:15] <- c("gene1.or.cancer", "gene2.or.cancer")

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


## Create cancer-related column for easy viewing
for (i in 1:nrow(rearrangements)){
    temp1 <- rearrangements$cancer_gene1[i]
    temp2 <- rearrangements$cancer_gene2[i]
    if (temp1 == "none"){
        if (temp2 == "none"){
            rearrangements[i, 16] <- "none"
        }else{
            rearrangements[i, 16] <- temp2
        }
    }else{
        if (temp2 == "none"){
            rearrangements[i, 16] <- temp1
        }else{
            rearrangements[i, 16] <- paste(temp1, temp2, sep = "_")
        }
    }
}
colnames(rearrangements)[16] <- "cancer.related"


## gets all genes
genes1 <- rearrangements[,c(1,2,3,14)]
genes2 <- rearrangements[,c(1,8,9,15)]
colnames(genes2)[2:4] <- colnames(genes1)[2:4]
genes.all <- rbind(genes1, genes2)
colnames(genes.all) <- c("sample", "chr", "pos", "gene")

## Plot rearrangments per gene
genes.all.4 <- ReccurentMaf(genes.all,"gene", 3)
genes.all.4 <- FilterMaf(genes.all.4, "Intergenic", "gene", F)
tbl <- table(genes.all.4$gene)
tbl <- tbl[order(tbl, decreasing = T)]
par(mar = c(8, 4,4,2))
barplot(tbl, las = 2, main = "# of rearrangments per gene, minimum 4")

## Plot unique rearrangments per gene
genes.all.unique <- PerSampleMaf(genes.all, "gene", identifier.column =  "sample")
genes.all.unique <- FilterMaf(genes.all.unique, "Intergenic", "gene", F)
genes.all.unique.2 <- ReccurentMaf(genes.all.unique, "gene")
tbl <- table(genes.all.unique.2$gene)
tbl <- tbl[order(tbl, decreasing = T)]
par(mar = c(8, 4,4,2))
barplot(tbl, las = 2, main = "# of unique rearrangments per gene, min 2")

## Plot recurrence-unique rearrangements per gene
genes.modified <- genes.all
genes.modified$sample <- sapply(genes.all$sample, function(x) substr(x,1,7), USE.NAMES = F)
genes.all.unique <- PerSampleMaf(genes.modified, "gene", identifier.column =  "sample")
genes.all.unique <- FilterMaf(genes.all.unique, "Intergenic", "gene", F)
genes.all.unique.2 <- ReccurentMaf(genes.all.unique, "gene")
tbl <- table(genes.all.unique.2$gene)
tbl <- tbl[order(tbl, decreasing = T)]
par(mar = c(8, 4,4,2))
barplot(tbl, las = 2, main = "# of unique rearrangments per gene, min 2")

##Investigate multiple hits
write.csv(genes.all.unique.2, "C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/Plots/rearrangement_targets.csv", row.names = F)


## Gets all recurrent rearrangments between pairs of genes. Replaces intergenic regions with nearby cancer genes
matches <- NULL
for (i in 1:nrow(gene.rearrangements)){
    if (gene.rearrangements$gene1[i] == "Intergenic"){
        if (gene.rearrangements$cancer_genes100kb_1[i] != "none"){
            gene.rearrangements$gene1[i] <- gene.rearrangements$cancer_genes100kb_1[i]
        }
    }
    if (gene.rearrangements$gene2[i] == "Intergenic"){
        if (gene.rearrangements$cancer_genes100kb_2[i] != "none"){
            gene.rearrangements$gene2[i] <- gene.rearrangements$cancer_genes100kb_2[i]
        }
    }
}

gene.rearrangements <- FilterMaf(rearrangements, "Intergenic", "gene1", F)
gene.rearrangements <- FilterMaf(gene.rearrangements, "Intergenic", "gene2", F)

for (i in 1:nrow(gene.rearrangements)){
    # Gets first gene-gene pair
    current.gene <- gene.rearrangements$gene1[i]
    gene.chr <- gene.rearrangements$chr1[i]
    gene.pos <- gene.rearrangements$pos1[i]
    current.target <- gene.rearrangements$gene2[i]
    target.chr <- gene.rearrangements$chr2[i]
    target.pos <- gene.rearrangements$pos2[i]
    youve.got.cancer <- gene.rearrangements$cancer.related[i]
    friends <- FilterMaf(gene.rearrangements[-i,], current.gene, "gene1")
    
    #checks to see if any other instances of same starter gene have same target
    if (nrow(friends) > 0){
        idx <- friends$gene2 %in% current.target
        close.friends <- friends[idx,]
        if (nrow(close.friends) > 0){
            temp <- c(rownames(gene.rearrangements)[i], current.gene, gene.chr, gene.pos, current.target, 
                      target.chr, target.pos, gene.rearrangements$Sample[i], youve.got.cancer)
            matches <- c(matches, temp)
            
        }
        
    } 
}

matches.table <-matrix(matches, length(matches)/9, 9, byrow = T)
pairwise.matches <- data.frame(matches.table, stringsAsFactors = FALSE)
colnames(pairwise.matches) <- c("Original Index", "Gene1", "Chr1", "Pos1", "Gene2", "Chr2", "pos2", "sample", "cancer")
x <- pairwise.matches

## gets rid of duplicates from same sample
for (i in 1:nrow(x)){
    ## Sets current row info
    current.gene <- x[, 2][i]
    current.target <- x[, 5][i]
    current.individual <- x[, 8][i]
    
    # Finds all samples with same start
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
                    x[as.numeric(rownames(close.friends)[1]), 10] <- 1
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
                    x[as.numeric(rownames(close.friends)[j]), 10] <- 1
                }
                
            }
                    
        }
        
    }
    
} 


pretty <- x[is.na(x[[10]]), ]
pretty <- pretty[, c(2:9)]
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


