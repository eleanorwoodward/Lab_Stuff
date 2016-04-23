## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")


## Generate separate mafs for each sample to be analzyed
men97 <- FilterMaf(disc.snindels, "MEN0097-P", "Tumor_Sample_Barcode")
men971 <- FilterMaf(disc.snindels, "MEN0097-P1", "Tumor_Sample_Barcode")
men972 <- FilterMaf(disc.snindels, "MEN0097-P2", "Tumor_Sample_Barcode")
men973 <- FilterMaf(disc.snindels, "MEN0097-P3", "Tumor_Sample_Barcode")

input.list <- list(men97, men971, men972, men973)

men97.het.comp <- NULL

## loops through list of mafs, annotating master maf with presence or absence for each unique mutation
## for each sample maf that is given in input list
for (i in 1:length(input.list)){
    print("i is", i)
    if (i == 1){
        ## sets up master list with all contents of first maf
        men97.het.comp <- input.list[[i]][, c(1, 5)]
        men97.het.comp[, 3] <- 1
        colnames(men97.het.comp)[3] <- "Men97"
        rownames(men97.het.comp) <- 1:nrow(men97.het.comp)
    }else{
        # populate column with zero as default, to be modified if any of the present mutations are found
        men97.het.comp[, i + 2] <- 0
        colnames(men97.het.comp)[i + 2] <- input.list[[i]][1,3]
        ## Loops through sample maf, checking each individual row
        for (j in 1:nrow(input.list[[i]])){
            tmp <- input.list[[i]]
            matches.idx <- men97.het.comp$Hugo_Symbol == tmp[j, 1]
            ## checks if gene found anywhere
            if (!is.na(matches.idx) & sum(matches.idx) > 0){
                matches <- men97.het.comp[matches.idx, ]
                hits <- which(matches$Start_position %in% tmp[j,5])
                ## checks if gene hits have same start position
                if (length(hits) > 0){
                    original.idx <- as.numeric(rownames(matches)[hits[1]])
                    men97.het.comp[original.idx, i +2] <- 1
                ## if not, adds it to master list    
                }else{
                    men97.het.comp <- rbind(men97.het.comp, c(tmp[j, 1], tmp[j, 5], rep(0, i -1), 1))
                }
                    
            }else{
                men97.het.comp <- rbind(men97.het.comp, c(tmp[j, 1], tmp[j, 5], rep(0, i -1), 1))
            }
        }
    }
}

men97.het.comp <- men97.het.comp[order(men97.het.comp$Men97, men97.het.comp$'MEN0097-P1',men97.het.comp$'MEN0097-P2', 
                                       men97.het.comp$'MEN0097-P3', decreasing = T), ]
