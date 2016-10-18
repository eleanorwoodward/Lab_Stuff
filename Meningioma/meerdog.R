## Simplified Meerkat style rearrangement caller

#rtrack <- import("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/psuedogenes.bed", format = "BED")
## read in regions of genome with tandem repeats
vntr <- readRDS("C:/Users/Noah/Documents/Big Files/dbs/gr.repeatMasker.rds")
## read in regions with transposable elements
te <- readRDS("C:/Users/Noah/Documents/Big Files/dbs/tubio_l1.rds")


events.matrix <- rearrangements
events.matrix$meerdog <- "not classified"
events.backup <- events.matrix
events.matrix <- events.backup

meerdog <- function(events.matrix){
    for (i in 1:nrow(events.matrix)){
        
        # sets variables for current rearrangement
        pos1 <- events.matrix[i, "pos1"]
        chr1 <- events.matrix[i, "chr1"]
        pos2 <- events.matrix[i, "pos2"]
        chr2 <- events.matrix[i, "chr2"]
        
        ## First check if either end of rearrangment is located within transposable element
        ## Jeremiah suggested we add a buffer here
        te.candidates <- te[te@ranges@start < pos1 & (te@ranges@width + te@ranges@start) > pos1 & te@seqnames == chr1]
        te.candidates2 <- te[te@ranges@start < pos2 & (te@ranges@width + te@ranges@start) > pos2 & te@seqnames == chr2]
        
        if (length(start(ranges(te.candidates))) > 0){
            events.matrix$meerdog[i] <- "TE"
            
        }else if (length(start(ranges(te.candidates2))) > 0){
            events.matrix$meerdog[i] <- "TE"
            
        }else{
            # Next check if overlaps with Variable Number of Tandem Repeats region at both ends
            vntr.candidates <- vntr[vntr@ranges@start < pos1 & (vntr@ranges@width + vntr@ranges@start) > pos1 & vntr@seqnames == chr1]
            vntr.candidates2 <- vntr[vntr@ranges@start < pos2 & (vntr@ranges@width + vntr@ranges@start) > pos2 & vntr@seqnames == chr2]
            vntr.candidates.combined <- vntr[vntr@ranges@start < pos1 & (vntr@ranges@width + vntr@ranges@start) > pos2 & vntr@seqnames == chr2]
        
            if (length(start(ranges(vntr.candidates.combined))) > 0 & chr1 == chr2){
            #if (length(start(ranges(vntr.candidates))) > 0 & length(start(ranges(vntr.candidates2))) > 0 & chr1 == chr2){
                events.matrix$meerdog[i] <- "VNTR"
               
            # }else if (events.matrix[i]$event.type == "insertion"){
            #     # insertions, but not deletions, get binned as NA if they don't match either of these: can ignore
            #     # for now until we integreate copy number plots with arrow calls
            #     events.matrix[i]$rearrangement.classification <- "NA"
            # 
            # 
            
            }else{
                ## They next check if deletion breakpoints have insertion at the same location. 
                if (!is.na(events.matrix$INSERTION[i]) & nchar(events.matrix$INSERTION[i]) > 0){
                
                    if (nchar(events.matrix$INSERTION[i]) > 10){
                        ## If insertion is > 10 bp vs < 10 bp, diferent calls
                        events.matrix$meerdog[i] <- "MMBIR"
                        
                    }else{
                        events.matrix$meerdog[i] <- "NHEJ"
                    }
                
                
                }else if (!is.na(events.matrix$HOMSEQ[i])){
                    # Checks for differing levels of homology at breakpoint
                    
                    if (nchar(events.matrix$HOMSEQ[i]) > 100){
                        ## first checks for broad homology: impossible with short reads
                        events.matrix$meerdog[i] <- "NAHR"
                    
                        
                    }else if (nchar(events.matrix$HOMSEQ[i]) > 19){
                        ## next checks medium levels of homology 
                        events.matrix$meerdog[i] <- "SSR"
                    
                    
                    }else if (nchar(events.matrix$HOMSEQ[i]) > 1){
                        ## checks for microhomology
                        events.matrix$meerdog[i] <- "MMEJ"
                    }else{
                        ## otherwise must be NHEJ
                        events.matrix$meerdog[i] <- "NHEJ"
                    }    
                
                }else{
                    ## otherwise must be NHEJ
                    events.matrix$meerdog[i] <- "NHEJ"
                }
            }
        }
    }
    return(events.matrix)
}

x <- meerdog(rearrangements)
x <- meerdog(events.matrix)
