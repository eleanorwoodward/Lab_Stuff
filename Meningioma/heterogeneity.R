## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

copy.number.calls <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/all_hg.7.15/broad_values_by_arm.txt", stringsAsFactors = F)
primary.tumors <- c("tum1")
recurrent.tumors <- c("dd")


## Generate separate mafs for each sample to be analzyed


input.list <- list(c("MEN0030.TumorA", "MEN0030.TumorB"), c("MEN0042.TumorA", "MEN0042.TumorB", "MEN0042.TumorC"), 
                   c("MEN0045.TumorA", "MEN0045.TumorB", "MEN0045.TumorC", "MEN0045.TumorD", "MEN0045.TumorE"), c("MEN0048.TumorA", "MEN0048.TumorB", "MEN0048.TumorC", "MEN0048.TumorD"),
                   c("MEN0093.TumorA", "MEN0093.TumorC", "MEN0093.TumorD", "MEN0093.TumorE"), c("MEN0097.Tumor", "MEN0097.TumorA", "MEN0097.TumorB", "MEN0097.TumorC"),
                   c("MEN0101.Tumor", "MEN0101.TumorB"), c("MEN0118.TumorA", "MEN0118.TumorB"), c("MEN0119.TumorA", "MEN0119.TumorB"), c("MEN0120.Tumor", "MEN0120.TumorB"))

combined.mutations <- NULL

## loops through list of mafs, annotating master maf with presence or absence for each unique mutation
for (i in 1:length(input.list[[5]])){
    ## sets up relevant info for each sample
    sample.name <- input.list[[5]][i]
    pair.name <- master.table[master.table$Tumor.Name == sample.name, ]$Pair.Name
    mutations <- FilterMaf(disc.snindels.duplicates, pair.name,"Tumor_Sample_Barcode")
    if (i == 1){
        ## sets up master list with all contents of first maf
        combined.mutations <- mutations[, c(1, 5)]
        combined.mutations[, 3] <- 1
        colnames(combined.mutations)[3] <- sample.name
        rownames(combined.mutations) <- 1:nrow(combined.mutations)
    }else{
        # populate column with zero as default, to be modified if any of the previously added mutations are found in current sample
        combined.mutations[, i + 2] <- 0
        colnames(combined.mutations)[i + 2] <- sample.name
        ## Loops through sample maf, checking each individual row
        for (j in 1:nrow(mutations)){
            tmp <- mutations
            matches.idx <- combined.mutations$Hugo_Symbol == tmp[j, 1]
            ## checks if gene found anywhere
            if (!is.na(matches.idx) & sum(matches.idx) > 0){
                matches <- combined.mutations[matches.idx, ]
                hits <- which(matches$Start_position %in% tmp[j,5])
                ## checks if gene hits have same start position
                if (length(hits) > 0){
                    original.idx <- as.numeric(rownames(matches)[hits[1]])
                    combined.mutations[original.idx, i +2] <- 1
                ## if not, adds it to master list    
                }else{
                    combined.mutations <- rbind(combined.mutations, c(tmp[j, 1], tmp[j, 5], rep(0, i -1), 1))
                }
                    
            }else{
                combined.mutations <- rbind(combined.mutations, c(tmp[j, 1], tmp[j, 5], rep(0, i -1), 1))
            }
        }
    }
}

combined.mutations <- cbind(combined.mutations, rowSums(data.matrix(combined.mutations[, -(1:2)])))

combined.mutations <- combined.mutations[order(combined.mutations[, 3], combined.mutations[, 4], combined.mutations[, 5], combined.mutations[, 6], decreasing = T), ]

combined.mutations <- combined.mutations[combined.mutations$`rowSums(data.matrix(combined.mutations[, -(1:2)]))` != 1, ]
