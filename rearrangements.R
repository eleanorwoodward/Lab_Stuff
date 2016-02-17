source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

## Gets start and end of each gene
hg.reference <- read.delim("C:/Users/Noah/Onedrive/Work/Coding/R/dbs/all_genes_hg38.txt", stringsAsFactors = F, comment.char = "", header = T)
reference <- hg.reference[, c(3, 5,6,13)]

## Gets centromere position
centromere.positions <- read.delim("C:/Users/Noah/Onedrive/Work/Coding/R/dbs/centromere_positions.txt", stringsAsFactors = F, comment.char = "", header = T)
centromere.positions <- FilterMaf(centromere.positions, c("chrX", "chrY"), "chrom", F)
mtrx <- matrix(NA, 22, 3)
colnames(mtrx) <- c("chr", "start", "end")

# takes max and min of possible centromere positions
for (i in 1:22){
    tmp <- unique(centromere.positions$chrom)[i]
    vals <- FilterMaf(centromere.positions, tmp, "chrom")
    mtrx[i,1] <- i
    mtrx[i,2] <- min(vals$chromStart)
    mtrx[i,3] <- max(vals$chromEnd)
}
cents <- data.frame(mtrx)

## Read in data, remove MEN049 which is not a meningioma
data <- read.delim("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/dRanger.txt", stringsAsFactors = F)
data <- FilterMaf(data, "MEN0049", "individual", F)

## gets all genes
genes1 <-data[,c(1,3,5, 24)]
genes2 <- data[,c(1,6, 8, 26)]
colnames(genes2)[2:4] <- c("chr1", "pos1", "gene1")
genes.all <- rbind(genes1, genes2)

## Plot total rearrangments per gene
genes.all.4 <- ReccurentMaf(genes.all,"gene1", 3)
tbl <- table(genes.all.4$gene1)
tbl <- tbl[order(tbl, decreasing = T)]
barplot(tbl, las = 2, main = "# of rearrangments per gene, minimum 4")

## Plot unique rearrangments per gene
genes.all.unique <- PerSampleMaf(genes.all, "gene1", identifier.column =  "individual")
genes.all.unique.2 <- ReccurentMaf(genes.all.unique, "gene1")
PlotMaf(genes.all.unique.2, "gene1", title = "Genes rearranged in at least 2 unique samples")

## Find rearrangmentes within region bordering genes
regions <- c()
names <- c()
interval <- 10000

## Loops through gene list
for (i in 1:nrow(genes.all.unique)){
    cur <- genes.all.unique$gene1[i]
    chr <- genes.all.unique$chr1[i]
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
        
        squad_size <- sum(genes.all.unique$pos1 > min & genes.all.unique$pos1 < max & genes.all.unique$chr1 == chr)
        regions <- c(regions, squad_size)
        names <- c(names, squad)
    }
}
names(regions) <- names
regions <- regions[order(regions, decreasing = T)]
barplot(regions[1:50], las = 2, main = "# of rearrangments per region")

# Find rearrangmentes within given arbitrary segments
regions.agnostic <- c()
names.agnostic <- c()
chr.agnostic <- c()
interval <- 200000

## Find max length of each chromosome
the.buck.stops.here <- c()
chromosome.names <- sort(unique(genes.all.unique$chr1))
for (i in 1:length(chrosome.names)){
    temp <- FilterMaf(genes.all.unique, chromosome.names[i], "chr1")
    the.buck.stops.here <- c(the.buck.stops.here, max(temp$pos1))
}

## Loops through chromosomes, breaking up into even regions
for (i in 1:22){
    universe <- FilterMaf(genes.all.unique, i, "chr1")
    
    for (j in 1:(the.buck.stops.here[i] / interval)){
        
        ## sets min and max based on values
        min <- (j - 1) * interval
        max <- min + interval
        squad_size <- sum(genes.all.unique$pos1 > min & genes.all.unique$pos1 < max)
        
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
    current.gene <- data$gene1[i]
    gene.chr <- data$chr1[i]
    gene.pos <- data$pos1[i]
    current.target <- data$gene2[i]
    target.chr <- data$chr2[i]
    target.pos <- data$pos2[i]
    friends <- FilterMaf(data[-i,], current.gene, "gene1")
    
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
x <- data.frame(matches.table, stringsAsFactors = FALSE)
colnames(x) <- c("Original Index", "Gene1", "Chr1", "Pos1", "Gene2", "Chr2", "pos2", "sample")
x2 <- x
## gets rid of duplicates from same sample

x <- FilterMaf(x, c("LAMA4", "NF2", "ACO1"), "X2")

for (i in 1:nrow(x)){
  current.gene <- x[, 2][i]
  current.target <- x[, 3][i]
  current.individual <- x[, 4][i]
  friends <- FilterMaf(x[-i,], current.gene, "X2")
  if (nrow(friends) > 0){
    idx <- friends[, 3] %in% current.target
    close.friends <- friends[idx,]
    
    if (nrow(close.friends) > 0){
      for (j in 1:nrow(close.friends)){
          if (close.friends[, 4][j] == current.individual){
              x[as.numeric(rownames(close.friends)[j]), 5] <- 1
          }
              
      }
          
      }
      
      
    }
    
} 

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

counter
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
for (i in 2:length(names)){
    full.names <- c(full.names, paste(names[i], "_Gistic", sep = ""))
    full.names <- c(full.names, paste(names[i], "_Rearrangments", sep = ""))
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

for (i in 1:length(unique(data$individual))){
    ind <- FilterMaf(data, unique(data$individual)[i],"individual")
    
}