## Simplified Meerkat style rearrangement caller


## Psuedo code

for (i in 1:nrow(events.matrix)){
    
    ## They first check if rearrangment location is covered by transposable element
    ## would best way to do this be with one of the tracks in GISTIC?
    if (events.matrix[i]$base.pair.position %in% TE.track){
        events.matrix[i]$rearrangment.classification <- "TE"
    
    ## They next check if overlaps with Variable Number of Tandem Repeats region   
    }else if (events.matrix[i]$base.pair.position %in% VNTR.track){
        events.matrix[i]$rearrangement.classification <- "VNTR"
        
    ## After first two options exhausted, remaining insertions get binned as NA
    ## Would SVtype be where i find categorization as insertion or deletion? Or should manual step be added here? Seems like it
    ## would be pretty straight forward to classify simple events as insertion or deletion based on DNA orientation info that plots
    ## are putting out
    }else if (events.matrix[i]$event.type == "insertion"){
        events.matrix[i]$rearrangement.classification <- "NA"
    
    ## They next check if deletion breakpoints have insertion at the same location. This is just your insertion category
    }else if (events.matrix[i]$multiple.breakpoints == T){
        
        ## If deletion breakpoint is > 10 bp vs < 10 bp, diferent calls
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
        
}