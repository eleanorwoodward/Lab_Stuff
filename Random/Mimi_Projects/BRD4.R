## Analysis for BRD4 paper

## Compare the change in weight of methyl marks for genes which scored in the ORF screen for all three marks
D283.orf <- read.delim("C:/Users/Noah/Dropbox/Work/ORF-results/D283_hits.txt", stringsAsFactors = F)
D283.orf.hits <- D283.orf[D283.orf$AVERAGE.JQ1 > 1 | D283.orf$AVERAGE.IBET > 1, -3]
colnames(D283.orf.hits)[1] <- "Gene.Symbol"

D458.orf <- read.delim("C:/Users/Noah/Dropbox/Work/ORF-results/D458_hits.txt", stringsAsFactors = F)
D458.orf.hits <- D458.orf[D458.orf$AVERAGE.JQ1 > 1 | D458.orf$AVERAGE.IBET > 1, ]

## combine orf hits from both screens together
orf.hits <- rbind(D458.orf.hits, D283.orf.hits)


## read in peak calls
D425.sensitive.BRD4.peaks <- read.delim("C:/Users/Noah/OneDrive/Work/Mimi/BRD4_Project/Peak_Analysis/D425_sensitive_JQ1_BRD4.aligned.merged.calls", 
                        stringsAsFactors = F, skip = 1, header = F, comment.char = "#")

D425J1.resistant.BRD4.peaks <- read.delim("C:/Users/Noah/OneDrive/Work/Mimi/BRD4_Project/Peak_Analysis/D425J1_resistant_BRD4.aligned.merged.calls", 
                                   stringsAsFactors = F, skip = 1, header = F, comment.char = "#")

D458.sensitive.BRD4.peaks <- read.delim("C:/Users/Noah/OneDrive/Work/Mimi/BRD4_Project/Peak_Analysis/D458_sensitive_JQ1_BRD4.aligned.merged.calls", 
                                   stringsAsFactors = F, skip = 1, header = F, comment.char = "#")

D458J1.resistant.BRD4.peaks <- read.delim("C:/Users/Noah/OneDrive/Work/Mimi/BRD4_Project/Peak_Analysis/D458J1_resistant_BRD4.aligned.merged.calls", 
                                     stringsAsFactors = F, skip = 1, header = F, comment.char = "#")

peaks.names <- c("PeakID", "chr", "start", "end" , "strand", "NormalizedTagCount", "regionsize", "findPeaksScore", "TotalTags", "ControlTags(normalizedtoIPExperiment)", 
                 "FoldChangevsControl", "p-valuevscontrol", "ClonalFoldChange")
# rename columns
colnames(D425J1.resistant.BRD4.peaks) <- peaks.names
colnames(D425.sensitive.BRD4.peaks) <- peaks.names
colnames(D458J1.resistant.BRD4.peaks) <- peaks.names
colnames(D458.sensitive.BRD4.peaks) <- peaks.names

## gene locations
hg.reference <- read.delim("C:/Users/Noah/Documents/Big Files/dbs/all_genes_hg38.txt", stringsAsFactors = F, comment.char = "", header = T)

## label peak regions with nearest gene


## takes peak calls from two separate Chips, sums the calls for the same gene within the same sample, then compares levels pre and post
CollapseAndCompare <- function(sensitive.data, resistant.data){
    sensitive.data <- sensitive.data[!is.na(sensitive.data$gene), ]
    resistant.data <- resistant.data[!is.na(resistant.data$gene), ]
    sensitive.gene.list <- unique(sensitive.data$gene)
    resistant.gene.list <- unique(resistant.data$gene)
    
    ## initializes comparison data frame with first gene from first sample
    comparison <- data.frame(sensitive.data$gene[1], sensitive.data$chr[1], sum(sensitive.data[sensitive.data$gene == sensitive.data$gene[1], ]$TotalTags),
                             stringsAsFactors = F)
    colnames(comparison) <- c("gene", "chr", "tags.sensitive")
    
    ## first generates collapsed pre-treatment dataframe
    for (i in 2:length(sensitive.gene.list)){
        hits <- sensitive.data[sensitive.data$gene == sensitive.gene.list[i], ]
        comparison[i, 1:2] <- c(hits$gene[1], hits$chr[1])
        comparison[i, 3] <- sum(hits$TotalTags)
    }
    
    ## next look at all genes in resistant, add genes that scored to same place, else add a new column
    comparison$tags.resistant <- NA
    for (j in 1:length(resistant.gene.list)){
        hits <- resistant.data[resistant.data$gene == resistant.gene.list[j], ]
        idx <- which(comparison$gene == hits$gene[1])
        if (length(idx) == 1){
            comparison$tags.resistant[idx] <- sum(hits$TotalTags)
        }else {
            comparison <- rbind(comparison, c(hits$gene[1], hits$chr[1], NA, sum(hits$TotalTags)))
        }
    }
    return(comparison)
    
}


## maps Chip peak to nearest gene within max.distance
ChromPos2Gene <- function(input.data, max.distance){
    for (i in 1:nrow(input.data)){
        if(i %% 1000 == 0){
            print(paste("merging row" ,i))
        }
        chr <- paste("chr", input.data$chr[i], sep = "")
        start <- input.data$start[i]
        end <- input.data$end[i]
        targets <- hg.reference[hg.reference$chrom == chr, ]
        targets <- targets[targets$txStart + max.distance > start & targets$txEnd - max.distance < end, ]
        if (nrow(targets) > 0){
            distances <- abs(targets$txStart - ((start + end) / 2))
            minimum <- min(distances)
            closest <- which(distances == minimum)
            input.data$gene[i] <- targets[closest[1], "name2"]
            input.data$distances[i] <- minimum
            
        }
    }
    return(input.data)
}

## compare peaks
D425J1.resistant.BRD4.peaks <- D425J1.resistant.BRD4.peaks[D425J1.resistant.BRD4.peaks$chr %in% 1:22, ]
D425.sensitive.BRD4.peaks <- D425.sensitive.BRD4.peaks[D425.sensitive.BRD4.peaks$chr %in% 1:22, ]
D425J1.resistant.BRD4.peaks <- ChromPos2Gene(D425J1.resistant.BRD4.peaks, 500000)
D425.sensitive.BRD4.peaks <- ChromPos2Gene(D425.sensitive.BRD4.peaks, 500000)

D458J1.resistant.BRD4.peaks <- D458J1.resistant.BRD4.peaks[D458J1.resistant.BRD4.peaks$chr %in% 1:22, ]
D458J1.resistant.BRD4.peaks <- ChromPos2Gene(D458J1.resistant.BRD4.peaks, 500000)
D458.sensitive.BRD4.peaks <- D458.sensitive.BRD4.peaks[D458.sensitive.BRD4.peaks$chr %in% 1:22, ]
D458.sensitive.BRD4.peaks <- ChromPos2Gene(D458.sensitive.BRD4.peaks, 500000)


D425.compare <- CollapseAndCompare(D425.sensitive.BRD4.peaks, D425J1.resistant.BRD4.peaks)
D425.compare.hits <- D425.compare[D425.compare$gene %in% D458.orf.hits$Gene.Symbol, ]
D425.compare.hits.all <- D425.compare[D425.compare$gene %in% orf.hits$Gene.Symbol, ]


D458.compare <- CollapseAndCompare(D458.sensitive.BRD4.peaks, D458J1.resistant.BRD4.peaks)
D458.compare.hits <- D458.compare[D458.compare$gene %in% D458.orf.hits$Gene.Symbol, ]
D458.compare.hits.all <- D458.compare[D458.compare$gene %in% orf.hits$Gene.Symbol, ]
