## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")


## Generate separate mafs for each sample to be analzyed

methylation.folder <- ("C:/Users/Noah/OneDrive/Work/mimi/Methylation Data/homer.calls/inputco")
methylation.names <- NULL

## read in all data, add unique identifier, keep track of sample names
methylation.calls <- NULL
for (i in 1:length(list.files(methylation.folder))){
    temp <- read.delim(paste(methylation.folder, list.files(methylation.folder)[i], sep = "/"),
                       stringsAsFactors = F, comment.char = "#", header = F)
    
    temp[, 14] <-  strsplit(list.files(methylation.folder)[i], ".txt")[[1]][1]
    methylation.names <- c(methylation.names, temp[1,14])
    methylation.calls <- rbind(methylation.calls, temp)

}
colnames(methylation.calls) <- c("PeakID", "chr", "start", "end", "strand", "NormalizedTagCount", "regionsize"
       ,"findPeaksScore", "TotalTags(normalizedtoControlExperiment)",   "ControlTags", "FoldChangevsControl", "p-valuevsControl",
       "ClonalFoldChange", "identifier")
# removes mitochondrial chromosome information
methylation.calls <- methylation.calls[methylation.calls$chr %in% c(1:24, "X", "Y"), ]
aggregated.calls <- NULL

## loops through list of mafs, annotating master maf with presence or absence for each unique methyl mark
## for each sample maf that is given in input list
for (i in 1:length(methylation.names)){
    print(paste("analyzing sample # ", i))
    temp <- FilterMaf(methylation.calls, methylation.names[i], "identifier")
    if (i == 1){
        ## sets up master list with all calls in first maf, chr, start, end, and presence or absence
        aggregated.calls <- temp[, c(2:4)]
        aggregated.calls[, 4] <- 1
        aggregated.calls <- rbind(aggregated.calls, c("Y", 0, 0, 0))
        colnames(aggregated.calls) <- c("chr", "start", "end", methylation.names[i])
        rownames(aggregated.calls) <- 1:nrow(aggregated.calls)
    }else{
        # populate column with zero as default, to be modified if any of the previously added mutations are found in current sample
        aggregated.calls[, i + 3] <- 0
        colnames(aggregated.calls)[i + 3] <- methylation.names[i]
        ## Loops through sample maf, checking each individual row
        for (j in 1:nrow(temp)){
            if(j %% 100 == 0){
                print(paste("Analyzing row ", j, " sample ", i))
            }
            # restrict master database to those calls from same chromosome
            working.calls <- aggregated.calls[aggregated.calls$chr == temp$chr[j], ]
            
            ## remove all points that have an end position smaller than our current start position - 100
            working.calls <- working.calls[working.calls$end > (temp$start[j] - 100), ]
            
            ## remove all points that have a start position greater than current end position + 100
            working.calls <- working.calls[working.calls$start < (temp$end[j] + 100), ]
                
            ## checks if region is already in list
            if (nrow(working.calls) > 0){
                ## if present, add 1 in column of all overlaps
                original.idx <- as.numeric(rownames(working.calls))
                aggregated.calls[original.idx, i +3] <- 1
            }else{
                ## if not, adds it to master list 
                aggregated.calls <- rbind(aggregated.calls, c(temp[j, 2], temp[j, 3], temp[j,4], rep(0, i -1), 1))
                
            }
            
        }
    }
}

for (k in 1:nrow(aggregated.calls)){
    total <- sum(as.numeric(aggregated.calls[k, 4:15]))
    aggregated.calls[k, 16] <- total
}

write.csv(aggregated.calls, "C:/Users/Noah/OneDrive/Work/Mimi/Methylation Data/pre_post_calls.csv", row.names = F)
aggregated.calls <- read.csv("C:/Users/Noah/OneDrive/Work/Mimi/Methylation Data/pre_post_calls.csv", stringsAsFactors = F)

## Permutations testing. 
total.hits <- c()
for (i in 1:1000){
    permuted.matrix <- aggregated.calls
    for (j in 4:(ncol(permuted.matrix) - 1)){
        permuted.matrix[, j] <- sample(permuted.matrix[, j], nrow(permuted.matrix), 
                                       replace = T)
    }
    
    permuted.matrix[, 16] <- rowSums(permuted.matrix[, 4:15])
    
    pass <- sum(permuted.matrix[, 16] > 5)
    total.hits <- c(total.hits, pass)
}
hist(total.hits, xlab = "# of permutations", main = "Histogram of permutations testing")
table(aggregated.calls$V16)

ggplot(aggregated.calls, aes(x = V16)) + geom_histogram() +scale_y_log10() +scale_x_continuous(limits = c(0, 10), breaks = 1:10)
ggplot(aggregated.calls, aes(x = V16)) + geom_histogram(binwidth = 1, center = 0) +scale_y_log10() +scale_x_continuous(limits = c(0, 10), breaks = 1:10)
ggplot(aggregated.calls, aes(x = V16)) + geom_histogram() + scale_y_log10() +scale_x_continuous(limits = c(0, 10), breaks = 1:10)





