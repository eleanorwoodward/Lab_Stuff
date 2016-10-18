## Mutations analysis for paper

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

copy.number.calls <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/all_hg.7.15/broad_values_by_arm_updated.txt", stringsAsFactors = F)
copy.number.binary <- copy.number.calls
for(i in 2:ncol(copy.number.calls)){
    for(j in 1:nrow(copy.number.calls)){
        if(copy.number.calls[j,i] != 0){
            copy.number.binary[j,i] <- 1
        }
    }
}


## Generate separate mafs for each sample to be analzyed


input.list <- list(c("MEN0030.TumorA", "MEN0030.TumorB"), c("MEN0042.TumorA", "MEN0042.TumorB", "MEN0042.TumorC"), 
                   c("MEN0045.TumorA", "MEN0045.TumorB", "MEN0045.TumorC", "MEN0045.TumorD", "MEN0045.TumorE"), c("MEN0048.TumorA", "MEN0048.TumorD", "MEN0048.TumorB", "MEN0048.TumorC"),
                   c("MEN0093.TumorC", "MEN0093.TumorD", "MEN0093.TumorA", "MEN0093.TumorE"), c("MEN0097.Tumor", "MEN0097.TumorA", "MEN0097.TumorB", "MEN0097.TumorC"),
                   c("MEN0101.TumorB", "MEN0101.TumorA"), c("MEN0118.TumorA", "MEN0118.TumorB"), c("MEN0119.TumorA", "MEN0119.TumorB"), c("MEN0120.Tumor", "MEN0120.TumorB"), 
                   c("MEN0121.TumorA", "MEN0121.TumorB", "MEN0121.TumorC", "MEN0121.TumorD"), c("MEN0122_RBM.TumorA", "MEN0122_RBM.TumorB", "MEN0122_RBM.TumorC",
                                                                                                "MEN0122_RBM.TumorD", "MEN0122_RBM.TumorE"))

primary.list <- c(1, 1, 1, 1,2,1,1,1,1,1)

private.muts <- c()
ubiq.muts <- c()
percent.muts <- c()
private.scna <- 0
ubiq.scna <- 0
percent.scna <- c()
auto.gene.lists <- list()
men.mafs <- list()

## loops through list containing reccurrences from each sample
for (h in 1:length(input.list)){
    combined.mutations <- NULL
    
    ## loops through list of mafs, annotating master maf with presence or absence for each unique mutation
    for (i in 1:length(input.list[[h]])){
        ## sets up relevant info for each sample
        sample.name <- input.list[[h]][i]
        pair.name <- master.table[master.table$Tumor.Name == sample.name, ]$Pair.Name
        mutations <- FilterMaf(disc.snindels.duplicates, pair.name,"Tumor_Sample_Barcode")
        if (i == 1){
            ## sets up master list with all contents of first maf
            combined.mutations <- mutations[, c(1,2,4,5)]
            combined.mutations[, 5] <- 1
            colnames(combined.mutations)[5] <- sample.name
            rownames(combined.mutations) <- 1:nrow(combined.mutations)
        }else{
            # populate column with zero as default, to be modified if any of the previously added mutations are found in current sample
            combined.mutations[, i + 4] <- 0
            colnames(combined.mutations)[i + 4] <- sample.name
            ## Loops through sample maf, checking each individual row
            for (j in 1:nrow(mutations)){
                tmp <- mutations
                matches.idx <- combined.mutations$Hugo_Symbol == tmp[j, 1]
                ## checks if gene found anywhere
                if (!is.na(matches.idx[1]) & sum(matches.idx) > 0){
                    matches <- combined.mutations[matches.idx, ]
                    hits <- which(matches$Start_position %in% tmp[j,5])
                    ## checks if gene hits have same start position
                    if (length(hits) > 0){
                        ## updates mutation count an average allelic fraction
                        original.idx <- as.numeric(rownames(matches)[hits[1]])
                        combined.mutations[original.idx, i +4] <- 1
                        prev.samples <- sum(as.numeric((combined.mutations[original.idx, 5:(i+3)])))
                        combined.mutations[original.idx, 2] <- as.numeric(combined.mutations[original.idx, 2]) * (prev.samples / (prev.samples + 1)) + 
                            tmp[j,2] * 1/(prev.samples + 1)              
                    ## if not, adds it to master list    
                    }else{
                        combined.mutations <- rbind(combined.mutations, c(tmp[j, 1], tmp[j, 2], tmp[j, 4], tmp[j, 5], rep(0, i -1), 1))
                    }
                        
                }else{
                    combined.mutations <- rbind(combined.mutations, c(tmp[j, 1], tmp[j, 2], tmp[j, 4], tmp[j, 5], rep(0, i -1), 1))
                }
            }
        }
    }
    men.mafs[[h]] <- combined.mutations
    combined.mutations <- cbind(combined.mutations, rowSums(data.matrix(combined.mutations[, -(1:4)])))
    total.samples <- length(input.list[[h]])
    private.muts <- c(private.muts, combined.mutations[combined.mutations[,total.samples + 5] == 1, ][,2])
    ubiq.muts <- c(ubiq.muts, combined.mutations[combined.mutations[,total.samples + 5] == total.samples, ][,2])
    percent.muts <- c(percent.muts, combined.mutations[, total.samples + 5] / total.samples)
    
    ## gets names of all genes that are ubiquitous for change in allelic fractions plots
    ubiq <- max(combined.mutations[, ncol(combined.mutations)])
    auto.gene.lists[[h]] <- combined.mutations[combined.mutations[[ncol(combined.mutations)]] == ubiq, ]$Hugo_Symbol
    
    # scna.index <- colnames(copy.number.calls) %in% input.list[[h]]
    # scna.matrix <- copy.number.binary[, scna.index]
    # scna.matrix <- cbind(scna.matrix, rowSums(scna.matrix))
    # scna.matrix <- scna.matrix[scna.matrix[, total.samples + 1] != 0, ]
    # private.scna <- private.scna + sum(scna.matrix[, total.samples + 1] == 1)
    # ubiq.scna <- ubiq.scna + sum(scna.matrix[, total.samples + 1] == total.samples)
    # percent.scna <- c(percent.scna, scna.matrix[, total.samples + 1] / total.samples)
    
    combined.mutations <- combined.mutations[order(combined.mutations[, 5], combined.mutations[, 6], decreasing = T), ]

    #combined.mutations <- combined.mutations[combined.mutations$`rowSums(data.matrix(combined.mutations[, -(1:2)]))` != 1, ]
    #write.table(combined.mutations, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Phylogenetic Trees/MEN0048_mutations.tsv", sep = "\t", row.names = F, quote = F)
}
## read in bap1 data, force call, add to existing mafs

bap1.folder <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Bap1/mafs/"
bap1.list <- list.files(bap1.folder)
bap1.indel <- NA
for (i in (c(1,3,5,7))){
    temp <- read.delim(paste(bap1.folder, bap1.list[i], sep=""), stringsAsFactors = F, comment.char = "#")
    bap1.indel <- rbind(bap1.indel, temp)
}

bap1.indel.long <- NA
for (i in (c(9,11,13,15,17))){
    temp <- read.delim(paste(bap1.folder, bap1.list[i], sep=""), stringsAsFactors = F, comment.char = "#")
    bap1.indel.long <- rbind(bap1.indel.long, temp)
}


bap1.snp <- NA
for (i in (c(2,4,6,8))){
    temp <- read.delim(paste(bap1.folder, bap1.list[i], sep=""), stringsAsFactors = F, comment.char = "#")
    bap1.snp <- rbind(bap1.snp, temp)
}

bap1.snp.long <- NA
for (i in (c(10,12,14,16,18))){
    temp <- read.delim(paste(bap1.folder, bap1.list[i], sep=""), stringsAsFactors = F, comment.char = "#")
    bap1.snp.long <- rbind(bap1.snp.long, temp)
}
bap1.1 <- rbind(bap1.snp[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "i_tumor_f", "Chromosome", "Start_position")],
                bap1.indel[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "i_tumor_f", "Chromosome", "Start_position")])
    
bap1.2 <- rbind(bap1.snp.long[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "i_tumor_f", "Chromosome", "Start_position")],
               bap1.indel.long[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "i_tumor_f", "Chromosome", "Start_position")])

bap1.1 <- bap1.1[bap1.1$Variant_Classification %in% c(snp.variants, indel.variants, "silent"), ]
bap1.2 <- bap1.2[bap1.2$Variant_Classification %in% c(snp.variants, indel.variants, "silent"), ]


bap1.1.forced <- ForceCallMutations(bap1.1[-1, ], "Tumor_Sample_Barcode", "Hugo_Symbol", "Start_position")
bap1.1.forced <- cbind(bap1.1.forced[, 1:2], bap1.1.forced[, -7])

bap1.2.forced <- ForceCallMutations(bap1.2[-1, ], "Tumor_Sample_Barcode", "Hugo_Symbol", "Start_position")
bap1.2.forced <- cbind(bap1.2.forced[, 1:2], bap1.2.forced[,-8])
men.mafs[[11]] <- bap1.1.forced
men.mafs[[12]] <- bap1.2.forced

## compute sd of each list element mutations
all.stdv <- c()
for (i in 1:length(pair.mut.counts)){
    stdv <- sd(pair.mut.counts[[i]])
    all.stdv <- c(all.stdv, stdv / mean(pair.mut.counts[[i]]))
}
all.stdv

mean(all.stdv)
sd(unlist(pair.mut.counts)) / mean(unlist(pair.mut.counts))

## same for copy number
all.stdv <- c()
pair.scna.counts <- pair.scna.counts[-6]
for (i in 1:length(pair.scna.counts)){
    stdv <- sd(pair.scna.counts[[i]])
    all.stdv <- c(all.stdv, stdv / mean(pair.mut.counts[[i]]))
}
all.stdv

mean(all.stdv)
sd(unlist(pair.scna.counts)) / mean(unlist(pair.scna.counts))

t.test(primary.counts, recurrent.counts)
t.test(primary.scna, recurrent.scna)




## compare average number of mutations shared between any two biopsies from same patient from maf that has already been forcecalled
PatientAverage <- function(maf, filler){
    ## of columns before data starts
    percent.shared <- c()
    if (ncol(maf) <= filler + 1){
        #do nothing, only one sample
    }else{
        for (i in (filler + 1):(ncol(maf) - 1)){
            for (j in (i + 1):ncol(maf)){
                shared <- sum(maf[, i] > 0 & maf[, j] > 0)
                total <- sum(maf[, i] > 0 | maf[, j] > 0)
                percent.shared <- c(percent.shared, shared / total)
            }
        }
    return(percent.shared)
    }
}

## calculates the percent of private mutations that are in the sample, as well as maximum difference
PatientPrimaryRecurrent <- function(maf, filler){
    ## of columns before data starts
    percent.private.recurrence <- c()
    percent.first.last <- 0
    change <- c()
    
    if (ncol(maf) <= filler + 1){
        #do nothing, only one sample
    }else{
        for (i in (filler + 1):(ncol(maf) - 1)){
            private.primary <- sum(maf[, i] > 0 & maf[, i + 1] == 0)
            private.recurrence <- sum(maf[, i + 1] > 0 & maf[, i] == 0)
            
            percent.private.recurrence <- c(percent.private.recurrence, private.recurrence / (private.primary + private.recurrence))
            
            change <- c(change, sum(maf[, i + 1] > 0) / sum(maf[, i] > 0))

        }
        private.first <- sum(maf[, filler + 1] > 0)
        private.last <- sum(maf[, ncol(maf)] > 0)
        percent.first.last <- private.last / (private.last + private.first)
        
        return(list(percent.private.recurrence, percent.first.last, change))
    }
}

ForceCallMutations <- function(maf, sample.identifier, identifier.1, identifier.2 = identifier.1) {
    ## Creates maf for each sample
    samples <- unique(maf[[sample.identifier]])
    combined.mutations <- NA
    identifiers <- 1
    if (!is.na(identifier.2)){
        identifiers <- 2
    }
    for (i in 1:length(samples)){
        ## sets up relevant info for each sample
        sample.name <- samples[i]
        mutations <- FilterMaf(maf, sample.name, sample.identifier)
        if (i == 1){
            ## sets up master list with all contents of first maf
            combined.mutations <- mutations[, c(identifier.1, identifier.2, sample.identifier)]
            combined.mutations[, 3] <- 1
            colnames(combined.mutations)[3] <- sample.name
            rownames(combined.mutations) <- 1:nrow(combined.mutations)
        }else{
            # populate column with zero as default, to be modified if any of the previously added mutations are found in current sample
            combined.mutations[, i + 2] <- 0
            colnames(combined.mutations)[i + 2] <- sample.name
            ## Loops through sample maf, checking each individual row
            for (j in 1:nrow(mutations)){
                ## checks first identifier
                matches.idx <- combined.mutations[[identifier.1]] == mutations[j, identifier.1]
                ## checks if any matches
                if (!is.na(matches.idx[1]) & sum(matches.idx) > 0){
                    matches <- combined.mutations[matches.idx, ]
                    hits <- which(matches[[identifier.2]] %in% mutations[j,identifier.2])
                    ## checks if any matches from second identifier
                    if (length(hits) > 0){
                        ## updates sample if found
                        original.idx <- as.numeric(rownames(matches)[hits[1]])
                        combined.mutations[original.idx, i +2] <- 1
                        ## if not, adds it to master list    
                    }else{
                        combined.mutations <- rbind(combined.mutations, c(mutations[j, identifier.1], mutations[j, identifier.2], rep(0, i -1), 1))
                    }
                    
                }else{
                    combined.mutations <- rbind(combined.mutations, c(mutations[j, identifier.1], mutations[j, identifier.2], rep(0, i - 1), 1))
                }
            }
        }
    }
    return(combined.mutations)
}

test.drive <- disc.snindels.duplicates[disc.snindels.duplicates$Tumor_Sample_Barcode %in%
                                           c("MEN0048-P0", "MEN0048-P2", "MEN0048-P3", "MEN0048G-P1"), ]
test.drive <- test.drive[test.drive$Hugo_Symbol %in% c("ACSL6", "ADAMTS18", "BBS2", "ASB11"), ]
test.1 <- ForceCallMutations(test.drive, "Tumor_Sample_Barcode", "Hugo_Symbol", "Start_position")

test.1[1:335, 1] == combined.mutations[1:335, 1]

combined.mutations[, 1] %in% test.1[, 1]

sum(combined.mutations[, 8] == 1)
sum(test.1[,6] == 1)



## read in and preprocess all maf files for comparison
## start with those that haven't been force called
study7 <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations/Study7_per_patient_comparison_all_pts.txt", 
                     stringsAsFactors = F)
colnames(study7)[2] <- "biopsy"
study7.mafs <- list()
study7.pats <- unique(study7$ID)
for (i in 1:length(study7.pats)){
    study7.mafs[[i]] <- ForceCallMutations(study7[study7$ID == study7.pats[i], ], "biopsy", "chrom", "start")
}

study8 <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations/Study8_per_patient_comparison_all_pts.txt", 
                     stringsAsFactors = F)
study8.mafs <- list()
study8.pats <- unique(study8$Patient)
for (i in 1:length(study8.pats)){
    study8.mafs[[i]] <- ForceCallMutations(study8[study8$Patient == study8.pats[i], ], "Sample", "chrom", "start")
}

study9 <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations/Study9_per_patient_comparison_all_pts.txt", 
                     stringsAsFactors = F, skip = 1)
study9.mafs <- list()
study9.pats <- unique(study9$patient)
for (i in 1:length(study9.pats)){
    study9.mafs[[i]] <- ForceCallMutations(study9[study9$patient == study9.pats[i], ], 
                                           "Tumor_Sample_Barcode", "Chromosome", "Start_position")
}

## then do those that have been force called already
study2 <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations/Study2_per_patient_comparison_all_pts.txt", 
                     stringsAsFactors = F)
study2.pats <- unique(study2$Sample)
study2.percentages <- c()
study2.numbers <- c()
for (i in 1:length(study2.pats)){
    temp <- study2[study2$Sample == study2.pats[i], ]
    shared <- sum(rowSums(temp[, 5:6]) == 2)
    total <- nrow(temp)
    study2.percentages <- c(study2.percentages, shared/total)
    study2.numbers <- c(study2.numbers, total)
}

study5 <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations/Study5_per_patient_comparison_all_pts.txt", 
                     stringsAsFactors = F)
study5.pats <- unique(study5$Patient.ID)
study5.percentages <- c()
study5.numbers <- c()
for (i in 1:length(study5.pats)){
    temp <- study5[study5$Patient.ID == study5.pats[i], ]
    shared <- sum(temp[, 3] == 1)
    total <- nrow(temp)
    study5.percentages <- c(study5.percentages, shared/total)
    study5.numbers <- c(study5.numbers, total)
}
file.directory <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations"
file.list <- list.files(file.directory)

## Study 1
study1.files <- file.list[3:8]
study1.mafs <- list()
for (i in 1:length(study1.files)){
    study1.mafs[[i]] <- read.delim(paste(file.directory,study1.files[i], sep = "/"), stringsAsFactors = F)
}

## Study 3
study3.file.directory <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/ID_protected/"
study3.files <- list.files(study3.file.directory)
study3.files <- study3.files[-(1:12)]
study3.mafs <- list()
temp <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/ID_protected/EC-001-ECM1-AB.maf",
                   stringsAsFactors = F)
current.pat <- "000"
pat.num <- 0
for (i in 1:length(study3.files)){
    temp <- read.delim(paste(study3.file.directory, study3.files[i], sep = ""), stringsAsFactors = F)
    next.pat <- strsplit(study3.files[i], "-")[[1]]
    if (next.pat[2] == current.pat){
        previous <- study3.mafs[[pat.num]]
        current <- temp[, 328]
        combined <- cbind(previous, current)
        colnames(combined)[ncol(combined)] <- next.pat[3]
        study3.mafs[[pat.num]] <- combined
    }else{
        current.pat <- next.pat[2]
        pat.num <- pat.num + 1
        combined <- temp[, c(1,5,328)]
        colnames(combined)[3] <- next.pat[3]
        study3.mafs[[pat.num]] <- combined
    }
}

## Study 4
study4.files <- file.list[24:28]
study4.mafs <- list()
for (i in 1:length(study4.files)){
    temp <- read.delim(paste(file.directory,study4.files[i], sep = "/"), stringsAsFactors = F)
    colnames(temp)[c(4,8,9,10)] <- c("start", "s1", "s2", "s3")
    study4.mafs[[i]] <- temp[, c(3:4, 8:10)]
}

## Study 6
study6.files <- file.list[33:42]
study6.mafs <- list()
for (i in 1:length(study6.files)){
    study6.mafs[[i]] <- read.delim(paste(file.directory,study6.files[i], sep = "/"), stringsAsFactors = F)
}

## Study 10
study10.files <- file.list[11:20]
study10.mafs <- list()
for (i in 1:length(study10.files)){
    temp <- read.delim(paste(file.directory,study10.files[i], sep = "/"), stringsAsFactors = F)
    study10.mafs[[i]] <- temp[, -c(2:4)]
}

## Study 11
study11.files <- file.list[22:28]
study11.mafs <- list()
for (i in 1:length(study11.files)){
    temp <- read.delim(paste(file.directory,study11.files[i], sep = "/"), stringsAsFactors = F)
    study11.mafs[[i]] <- temp[, -c(1:2)]
}

## compute pairwise comparisons
study1.percentages <- c()
study1.averages <- c()
study1.numbers <- c()
study1.first.last <- c()
study1.primary.recurrence <- c()
for (i in 1:length(study1.mafs)){
    temp <- study1.mafs[[i]]
    percents <- PatientAverage(temp, 0)
    primary.vs.recurrent <- PatientPrimaryRecurrent(temp, 0)
    study1.primary.recurrence <- c(study1.primary.recurrence,primary.vs.recurrent[[1]])
    study1.first.last <- c(study1.first.last, primary.vs.recurrent[[2]])
    study1.percentages <- c(study1.percentages, percents)
    study1.averages <- c(study1.averages, mean(percents))
    study1.numbers <- c(study1.numbers, nrow(temp))
    
}

study3.percentages <- c()
study3.averages <- c()
study3.numbers <- c()
for (i in 1:length(study3.mafs)){
    temp <- study3.mafs[[i]]
    percents <- PatientAverage(temp, 2)
    study3.percentages <- c(study3.percentages, percents)
    study3.averages <- c(study3.averages, mean(percents))
    study3.numbers <- c(study3.numbers, nrow(temp))
}

study4.mafs[[5]] <- study4.mafs[[5]][-6440, ]
study4.percentages <- c()
study4.averages <- c()
study4.numbers <- c()
for (i in 1:length(study4.mafs)){
    temp <- study4.mafs[[i]]
    percents <- PatientAverage(temp, 2)
    study4.percentages <- c(study4.percentages, percents)
    study4.averages <- c(study4.averages, mean(percents))
    study4.numbers <- c(study4.numbers, nrow(temp))
}

temp <- study6.mafs[[9]]

study6.percentages <- c()
study6.averages <- c()
study6.numbers <- c()
for (i in 1:length(study6.mafs)){
    temp <- study6.mafs[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 1)
    study6.percentages <- c(study6.percentages, percents)
    study6.averages <- c(study6.averages, mean(percents))
    study6.numbers <- c(study6.numbers, nrow(temp))
}

temp <- study7.mafs[[5]]

study7.percentages <- c()
study7.averages <- c()
study7.numbers <- c()
for (i in 1:length(study7.mafs)){
    temp <- study7.mafs[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 2)
    study7.percentages <- c(study7.percentages, percents)
    study7.averages <- c(study7.averages, mean(percents))
    study7.numbers <- c(study7.numbers, nrow(temp))
}

study8.percentages <- c()
study8.numbers <- c()
for (i in 1:length(study8.mafs)){
    temp <- study8.mafs[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 2)
    study8.percentages <- c(study8.percentages, percents)
    study8.numbers <- c(study8.numbers, nrow(temp))
}



study9.percentages <- c()
study9.averages <- c()
study9.numbers <- c()
for (i in 1:length(study9.mafs)){
    temp <- study9.mafs[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 2)
    study9.percentages <- c(study9.percentages, percents)
    study9.averages <- c(study9.averages, mean(percents))
    study9.numbers <- c(study9.numbers, nrow(temp))
}


study10.percentages <- c()
study10.averages <- c()
study10.numbers <- c()
for (i in 1:length(study10.mafs)){
    temp <- study10.mafs[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 1)
    study10.percentages <- c(study10.percentages, percents)
    study10.averages <- c(study10.averages, mean(percents))
    study10.numbers <- c(study10.numbers, nrow(temp))
}


study11.percentages <- c()
study11.averages <- c()
study11.numbers <- c()
for (i in 1:length(study11.mafs)){
    temp <- study11.mafs[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 0)
    study11.percentages <- c(study11.percentages, percents)
    study11.averages <- c(study11.averages, mean(percents))
    study11.numbers <- c(study11.numbers, nrow(temp))
}

men.mafs.temp <- men.mafs[-c(3,6)]

## percent shared betweeen all possible pairwise combinations per patient
studyMen.percentages <- c()

## average percent shared between all possible combinations within each patient
studyMen.averages <- c()

## total number of called mutations from all samples from given patient
studyMen.numbers <- c()

## difference in private mutationsr between first and last recurrence
studyMen.first.last <- c()

## difference between any two sequential biopsies
studyMen.primary.recurrence <- c()
studyMen.primary.recurrence.average <- c()
studyMen.primary.recurrence.change <- c()

## change in total number of mutations from recurrence to recurrence
studyMen.primary.recurrence.change.absolute <- c()


mutation.deviation <- list()
for (i in 1:length(men.mafs)){
    temp <- men.mafs.temp[[i]]
    temp <- temp[, !is.na(temp[2, ])]
    temp <- temp[!is.na(temp[, 2]), ]
    percents <- PatientAverage(temp, 4)
    primary.vs.recurrent <- PatientPrimaryRecurrent(temp, 4)
    primary.vs.recurrent.vals <- primary.vs.recurrent[[1]]
    studyMen.primary.recurrence.change.absolute <- c(studyMen.primary.recurrence.change.absolute, primary.vs.recurrent[[3]])
    studyMen.primary.recurrence <- c(studyMen.primary.recurrence,primary.vs.recurrent.vals)
    studyMen.primary.recurrence.average <- c(studyMen.primary.recurrence.average,mean(primary.vs.recurrent.vals))
    studyMen.first.last <- c(studyMen.first.last, primary.vs.recurrent[[2]])
    studyMen.percentages <- c(studyMen.percentages, percents)
    studyMen.averages <- c(studyMen.averages, mean(percents))
    studyMen.numbers <- c(studyMen.numbers, nrow(temp))
    for (j in 2:(length(primary.vs.recurrent.vals))){
        percent <- (primary.vs.recurrent.vals[j] / primary.vs.recurrent.vals[j - 1])
        studyMen.primary.recurrence.change <- c(studyMen.primary.recurrence.change, percent)
    }

    # gets # of mutations
    num.muts <- c()
    for (j in 5:ncol(temp)){
        num.muts <- c(num.muts, sum(temp[, j] > 0))
    }
    mutation.deviation[[i]] <- num.muts
}

mean(studyMen.primary.recurrence.change.absolute)


## regenerate percent variance across samples
mutation.deviation.vals <- c()
for (i in 1:length(mutation.deviation)){
    mutation.deviation.vals <- c(mutation.deviation.vals, sd(mutation.deviation[[i]]) / mean(mutation.deviation[[i]]))
}

mean(mutation.deviation.vals)


## look at heterogeneity in copy number alterations
focal.cna <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/all_hg.7.15/all_lesions.conf_99.txt", 
                        stringsAsFactors = F)
focal.cna <- focal.cna[1:61, ]

focal.mafs <- list(focal.cna[, 11:12], focal.cna[, 14:16], focal.cna[, 17:21], focal.cna[, 22:25], 
                   focal.cna[, c(33, 35:37)], focal.cna[, 47:48], focal.cna[, 65:66], focal.cna[, 67:68], focal.cna[, 69:70])

focal.percents <- c()
focal.percents.mean <- c()
for (i in 1:length(focal.mafs)){
    temp <- focal.mafs[[i]]
    patient.percents <- PatientAverage(temp, 0)
    focal.percents <- c(focal.percents, patient.percents)
    focal.percents.mean <- c(focal.percents.mean, mean(patient.percents))
}

broad.cna <- cbind(focal.cna[1:39, 1:8], copy.number.binary)
broad.mafs <- list(broad.cna[, 11:12], broad.cna[, 14:16], broad.cna[, 17:21], broad.cna[, 22:25], 
                   broad.cna[, c(33, 35:37)], broad.cna[, 47:48], broad.cna[, 65:66], broad.cna[, 67:68], broad.cna[, 69:70])

broad.percents <- c()
broad.percents.mean <- c()
for (i in 1:length(broad.mafs)){
    temp <- broad.mafs[[i]]
    patient.percents <- PatientAverage(temp, 0)
    broad.percents <- c(broad.percents, patient.percents)
    broad.percents.mean <- c(broad.percents.mean, mean(patient.percents))
}


all.cna <- rbind(focal.cna[-9], broad.cna[, -9])
all.cna <- cbind(all.cna[,1], all.cna)

all.mafs <- list(all.cna[, 11:12], all.cna[, 14:16], all.cna[, 17:21], all.cna[, 22:25], 
                   all.cna[, c(33, 35:37)], all.cna[, 47:48], all.cna[, 65:66], all.cna[, 67:68], all.cna[, 69:70])

all.percents <- c()
all.percents.mean <- c()
for (i in 1:length(all.mafs)){
    temp <- all.mafs[[i]]
    patient.percents <- PatientAverage(temp, 0)
    all.percents <- c(all.percents, patient.percents)
    all.percents.mean <- c(all.percents.mean, mean(patient.percents))
}






## export percents for pan-can comparison
export.list.1 <- list(study1.percentages, study2.percentages, study3.percentages, study4.percentages, study5.percentages,
                           study6.percentages, study7.percentages, study8.percentages, study9.percentages, study10.percentages, 
                           studyMen.percentages)
export.list.2 <- list(study1.averages, study2.percentages, study3.averages, study4.averages, study5.percentages, study6.averages, 
                     study7.averages, study8.percentages, study9.averages, study10.averages, studyMen.averages)
max.length <- 0
for (i in 1:length(export.list.2)){
    if (length(export.list.2[[i]]) > max.length){
        max.length <- length(export.list.2[[i]])
    }
}

for (i in 1:length(export.list.2)){
    export.list.2[[i]] <- c(export.list.2[[i]], rep(NA, max.length - length(export.list.2[[i]])))
}

export.2 <- data.frame(export.list.2)
write.csv(export.2, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heterogeneity/shared.comparison.fixed_calculation_average.csv", 
          row.names = F)

all.numbers <- c(study1.numbers, study2.numbers, study3.numbers, study4.numbers, study5.numbers,
                                      study6.numbers, study7.numbers, study8.numbers, study9.numbers, study10.numbers, 
                                      studyMen.numbers)

all.averages <- c(study1.averages, study2.percentages, study3.averages, study4.averages, study5.percentages,
                 study6.averages, study7.averages, study8.percentages, study9.averages, study10.averages, 
                 studyMen.averages)

remove <- all.numbers > 4000

all.numbers <- all.numbers[!remove]
all.averages <- all.averages[!remove]
plot(all.averages, all.numbers)



## check and see about samples that don't have reconstructable phylogeny
men0048.combined.calls <- men.mafs[[4]]
men0048.combined.calls <- men0048.combined.calls[order(men0048.combined.calls[, 5], men0048.combined.calls[, 6], men0048.combined.calls[, 7], men0048.combined.calls[, 8], 
                                                       decreasing = T), ]

men0093.combined.calls <- men.mafs[[5]]
men0093.combined.calls <- men0093.combined.calls[order(men0093.combined.calls[, 5], men0093.combined.calls[, 6], men0093.combined.calls[, 7], men0093.combined.calls[, 8], 
                                                       decreasing = T), ]
write.csv(men0093.combined.calls, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Phylogenetic Trees/men0093.mutations.csv", row.names = F)

write.csv(disc.snindels, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/Heterogeneity/data/meningioma_allelic_fraction.csv", row.names = F)
## get info to make accurate trees
## MEn0030
tree <- men.mafs[[1]]
sum(tree[, 5] > 0 & tree[,6] > 0)
sum(tree[, 5] > 0)
sum(tree[, 6] > 0)

##MEN042
tree <- men.mafs[[2]]
tree <- tree[order(tree[, 5], tree[,6], tree[,7], decreasing = T), ]
sum(tree[, 5] > 0 & tree[,6] > 0 & tree[,7] > 0)
sum(tree[, 5] > 0 & tree[,7] > 0)
sum(tree[, 5] > 0)
sum(tree[, 6] > 0)
sum(tree[, 7] > 0)





## Figs
## create stacked bar plot for each patient with totally private mutations, totally shared mutations, and anything in between for each patient. 

## creats a db with clonality counts

sample.vector <- c()
clonal.vector <- c()
clonal.muts <- list()

for (i in 1:length(men.mafs)){
    maf <- men.mafs[[i]]
    maf.matrix <- as.matrix(maf[, -c(1:4)])
    maf.matrix <- apply(maf.matrix,2,as.numeric)
    clonal.vals <- rowSums(maf.matrix)
    samples <- ncol(maf.matrix)
    current.clonal <- rep("private", length(clonal.vals))
    current.clonal[clonal.vals == samples] <- "clonal"
    current.clonal[clonal.vals < samples & clonal.vals > 1] <- "subclonal"
    clonal.vector <- c(clonal.vector, current.clonal)
    sample.vector <- c(sample.vector, rep(colnames(maf)[5], nrow(maf)))
    clonal.muts[[i]] <- maf$Hugo_Symbol[current.clonal == "clonal"]
    
}

clonal.df <- data.frame(sample.vector, clonal.vector)
clonal.df$sample.vector <- factor(clonal.df$sample.vector, levels=names(sort(table(clonal.df$sample.vector), decreasing = T)))
clonal.df$clonal.vector <- factor(clonal.df$clonal.vector, levels = c("clonal", "subclonal", "private"))

ggplot(data=clonal.df, aes(x=sample.vector, fill=clonal.vector)) + geom_bar() + scale_fill_grey()

## load ccf data
absolute.folder <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/ABSOLUTE/mafs"
absolute.ccfs <- NULL
for (i in 1:length(list.files(absolute.folder))){
    temp <- readRDS(paste(absolute.folder, list.files(absolute.folder)[i], sep = "/"))
    absolute.ccfs <- rbind(absolute.ccfs, temp[, c("pair_id", "Hugo_Symbol", "Chromosome", "Start_position", "Variant_Classification", "i_tumor_f", "ccf_hat")])
}

absolute.ccfs <- absolute.ccfs[absolute.ccfs$Variant_Classification %in% c(snp.variants, indel.variants, "Silent"), ]
absolute.ccfs <- absolute.ccfs[!is.na(absolute.ccfs$i_tumor_f), ]


## plot evolution of allelic fractions over time
## takes a list of genes, and returns their allelic fraction for all samples in given sublist
input.list.temp <- input.list[-c(1,3,6)]
clonal.muts.temp <- clonal.muts[-c(1,3,6)]
time.series.list <- list()

## calculates change in ccf over time, both for average of all and for average of clonal mutations
ccf.lfc <- c()
ccf.clonal.lfc <- c()
ccf.first.last <- c()

for (i in 1:length(input.list.temp)){
    samples <- input.list.temp[[i]]
    fractions <- c()
    recurrence <- c()
    gene <- c()
    of.interest <- c()
    prior.ccf <- c()
    first.ccf <- 0
    for (j in 1:length(samples)){
        pair.name <- master.table[master.table$Tumor.Name == samples[j], ]$Pair.Name[1]
        mutations <- FilterMaf(absolute.ccfs, pair.name,"pair_id")
        mutations <- PerSampleMaf(mutations, "Hugo_Symbol", identifier.column = "pair_id")
        average <- mean(mutations$ccf_hat)
        mutations <- mutations[mutations$Hugo_Symbol %in% clonal.muts.temp[[i]], ]
        fractions <- c(fractions, mutations$ccf_hat, average)
        recurrence <- c(recurrence, rep(pair.name, nrow(mutations) + 1))
        of.interest <- c(of.interest, rep("genes", nrow(mutations)), "average")
        
        gene <- c(gene, mutations$Hugo_Symbol, "Average")
        if ("NF2" %in% mutations$Hugo_Symbol){
            idx <- max(which(gene %in% "NF2"))
            of.interest[idx] <- "NF2"
        }
        if (j == 1){
            first.ccf <- average
        }
        if (j > 1){
            ccf.lfc <- c(ccf.lfc, average / prior.ccf[1])
            ccf.clonal.lfc <- c(ccf.clonal.lfc, mean(mutations$ccf_hat / prior.ccf[2]))
        }
        
        if (j == length(samples)){
            ccf.first.last <- c(ccf.first.last, average / first.ccf)
        }
        prior.ccf <- c(average, mean(mutations$ccf_hat))
    }
    
    
    time.series <- data.frame(fractions, recurrence, gene, of.interest)
    # Map sex to color
    x_levels <- unique(recurrence)
    time.series$recurrence <- factor(time.series$recurrence, levels = x_levels)
    time.series.list[[i]] <- time.series
}

ggplot(data=time.series.list[[2]], aes(x=recurrence, y=fractions, group=gene, colour=of.interest)) +
    geom_line(size=1.5, aes(linetype=of.interest)) +
    geom_point(size=2.5) +
    xlab("Recurrence") +
    ylab("Allelic Fraction") +
    scale_linetype_manual(values=c("dotted", "solid", "twodash")) +
    scale_color_grey()


## look at 42 specifically, break up by direction of change
time.series1 <- time.series[0, ]
time.series2 <- time.series1
for (i in 1:length(unique(time.series$gene))){
    temp <- time.series[time.series$gene == unique(time.series$gene)[i], ]
    if (temp$fractions[2] > temp$fractions[1]){
        time.series1 <- rbind(time.series1, temp)
    }else{
        time.series2 <- rbind(time.series2, temp)
    }
}
time.series2 <- rbind(time.series2, time.series[time.series$gene == "Average", ])
time.series1$recurrence <- factor(time.series1$recurrence, levels = unique(time.series1$recurrence))
time.series2$recurrence <- factor(time.series2$recurrence, levels = unique(time.series2$recurrence))





