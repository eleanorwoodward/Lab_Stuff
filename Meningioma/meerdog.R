## Simplified Meerkat style rearrangement caller

rtrack <- import("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/psuedogenes.bed", format = "BED")
## read in regions of genome with tandem repeats
vntr <- readRDS("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/gr.repeatMasker.rds")
## read in regions with transposable elements
te <- readRDS("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/tubio_l1.rds")

for (i in 1:nrow(events.matrix)){
    
    # sets variables for current rearrangement
    pos1 <- events.matrix[i, "pos1"]
    chr1 <- events.matrix[i, "chr1"]
    pos2 <- events.matrix[i, "pos2"]
    chr2 <- events.matrix[i, "chr2"]
    
    ## First check if either end of rearrangment is located within transposable element
    ## buffer this one. 
    te.candidates <- rtrack[rtrack@ranges@start < pos1 & (rtrack@ranges@width + rtrack@ranges@start) > pos1 & rtrack@seqnames == chr1]
    te.candidates2 <- rtrack[rtrack@ranges@start < pos2 & (rtrack@ranges@width + rtrack@ranges@start) > pos2 & rtrack@seqnames == chr2]
    
    if (length(start(ranges(te.candidates))) > 0){
        events.matrix[i]$meerdog <- "TE"
        
    }else if (length(start(ranges(te.candidates2))) > 0){
        events.matrix[i]$meerdog <- "TE"
        
    }
    ## Next check if overlaps with Variable Number of Tandem Repeats region   
    vntr.candidates <- rtrack[rtrack@ranges@start < pos1 & (rtrack@ranges@width + rtrack@ranges@start) > pos1 & rtrack@seqnames == chr1]
    vntr.candidates2 <- rtrack[rtrack@ranges@start < pos2 & (rtrack@ranges@width + rtrack@ranges@start) > pos2 & rtrack@seqnames == chr2]
    
    if (length(start(ranges(vntr.candidates))) > 0){
        events.matrix[i]$meerdog <- "VNTR"
        
    }else if (length(start(ranges(vntr.candidates2))) > 0){
        events.matrix[i]$meerdog <- "VNTR"
        
     
    ## After first two options exhausted, remaining insertions get binned as NA
    }else if (events.matrix[i]$event.type == "insertion"){
        events.matrix[i]$rearrangement.classification <- "NA"
    
    ## They next check if deletion breakpoints have insertion at the same location. This is just your insertion category
    }else if (events.matrix[i]$insertion == T){
        
        ## If insertion is > 10 bp vs < 10 bp, diferent calls
        if (events.matrix[i]$breakpoint.length > 10){
            events.matrix[i]$rearrangement.classification <- "MMBIR"
        }else{
            events.matrix[i]$rearrangement.classification <- "NHEJ"
        }
    
    ## Checks for broad homology at breakpoint    
    }else if (events.matrix[i]$homology > 100){
        events.matrix[i]$rearrangement.classification <- "NAHR"
        
    ## Checks for microhomology at breakpoint    
    }else if (events.matrix[i]$homology > 1){
        events.matrix[i]$rearrangement.classification <- "MMEJ"
    
    ## otherwise must be NHEJ
    }else{
        events.matrix[i]$rearrangment.classification <- "NHEJ"
    }
    
}else{
    events.matrix$meerdog[i] <- "not evaluated"

            
}