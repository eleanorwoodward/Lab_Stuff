## Simplified Meerkat style rearrangement caller

#rtrack <- import("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/psuedogenes.bed", format = "BED")
## read in regions of genome with tandem repeats
vntr <- readRDS("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/gr.repeatMasker.rds")
## read in regions with transposable elements
te <- readRDS("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/tubio_l1.rds")


folder <-("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/Rearrangements")
events.matrix <- NULL
for (i in 1:length(list.files(folder))){
    temp <- read.csv(paste(folder, list.files(folder)[i], sep = "/"),
                     stringsAsFactors = F)
    temp[, 28] <-  strsplit(list.files(folder)[i], ".csv")[[1]]
    colnames(temp)[28] <- "Sample"
    colnames(temp)[29] <- "vcf.info"
    events.matrix <- rbind(events.matrix, temp)    
}



events.matrix <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Snowman/Rearrangements/men.tsv", 
                            stringsAsFactors = F)
events.matrix$meerdog <- "not classified"
for (i in 1:nrow(events.matrix)){
    
    # sets variables for current rearrangement
    pos1 <- events.matrix[i, "pos1"]
    chr1 <- events.matrix[i, "chr1"]
    pos2 <- events.matrix[i, "pos2"]
    chr2 <- events.matrix[i, "chr2"]
    
    ## First check if either end of rearrangment is located within transposable element
    ## buffer this one. 
    te.candidates <- te[te@ranges@start < pos1 & (te@ranges@width + te@ranges@start) > pos1 & te@seqnames == chr1]
    te.candidates2 <- te[te@ranges@start < pos2 & (te@ranges@width + te@ranges@start) > pos2 & te@seqnames == chr2]
    
    if (length(start(ranges(te.candidates))) > 0){
        events.matrix$meerdog[i] <- "TE"
        next
        
    }else if (length(start(ranges(te.candidates2))) > 0){
        events.matrix$meerdog[i] <- "TE"
        next
    }
    ## Next check if overlaps with Variable Number of Tandem Repeats region, either at one end or both ends
#     vntr.candidates <- vntr[vntr@ranges@start < pos1 & (vntr@ranges@width + vntr@ranges@start) > pos1 & vntr@seqnames == chr1]
#     vntr.candidates2 <- vntr[vntr@ranges@start < pos2 & (vntr@ranges@width + vntr@ranges@start) > pos2 & vntr@seqnames == chr2]
#     
#     if (length(start(ranges(vntr.candidates))) > 0){
#         if (length(start(ranges(vntr.candidates2))) > 0){
#             events.matrix$meerdog[i] <- "VNTR-VNTR"
#         }else{
#         events.matrix$meerdog[i] <- "VNTR"
#         }
#     }else if (length(start(ranges(vntr.candidates2))) > 0){
#         events.matrix$meerdog[i] <- "VNTR"
#         
    ## insertions, but not deletions, get binned as NA if they don't match either of these: can ignore 
    ## for now until we integreate copy number plots with arrow calls
#   }else if (events.matrix[i]$event.type == "insertion"){
#       events.matrix[i]$rearrangement.classification <- "NA"
#     

    ## They next check if deletion breakpoints have insertion at the same location.
    #}else 
        if (!is.na(events.matrix$INSERTION[i]) & nchar(events.matrix$INSERTION[i]) > 0){
        
        ## If insertion is > 10 bp vs < 10 bp, diferent calls
        if (nchar(events.matrix$INSERTION[i]) > 10){
            events.matrix$meerdog[i] <- "MMBIR"
        }else{
            events.matrix$meerdog[i] <- "NHEJ"
        }
    
    ## Checks for broad homology at breakpoint: can't do with short reads  
#     }else if (events.matrix[i]$homology > 100){
#         events.matrix[i]$rearrangement.classification <- "NAHR"
#         
    ## Checks for microhomology at breakpoint    
    }else if (!is.na(events.matrix$HOMSEQ[i]) & nchar(events.matrix$HOMSEQ[i]) > 1){
        events.matrix$meerdog[i] <- "MMEJ"
    
    ## otherwise must be NHEJ
    }else{
        events.matrix$meerdog[i] <- "NHEJ"
    }
    
}