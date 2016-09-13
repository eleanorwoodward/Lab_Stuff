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
                   c("MEN0101.Tumor", "MEN0101.TumorA"), c("MEN0118.TumorA", "MEN0118.TumorB"), c("MEN0119.TumorA", "MEN0119.TumorB"), c("MEN0120.Tumor", "MEN0120.TumorB"))

primary.list <- c(1, 1, 1, 1,2,1,1,1,1,1)

private.muts <- c()
ubiq.muts <- c()
percent.muts <- c()
private.scna <- 0
ubiq.scna <- 0
percent.scna <- c()
auto.gene.lists <- list()
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

## calculate percentages of mutation types
length(private.muts) / length(percent.muts)
length(ubiq.muts) / length(percent.muts)
mean(percent.muts)

private.scna /length(percent.scna)
ubiq.scna / length(percent.scna)
mean(percent.scna)

# calculate difference in allelic fractions
t.test(as.numeric(private.muts), as.numeric(ubiq.muts))

pair.mut.counts <- vector("list", length(input.list))
pair.scna.counts <- vector("list", length(input.list))
primary.counts <- c()
primary.scna <- c()
recurrent.counts <- c()
recurrent.scna <- c()

## summary calculations for section
for (i in 1:length(input.list)){
    sample.names <- input.list[[i]]
    current.counts <- c()
    current.scna.counts <- c()
    for (j in 1:length(sample.names)){
        current.name <- master.table[master.table$Tumor.Name == sample.names[j], ]$Pair.Name
        muts <- nrow(FilterMaf(disc.snindels.duplicates, current.name,"Tumor_Sample_Barcode"))
        current.counts <- c(current.counts, muts)
        scnas <- sum(copy.number.calls[[sample.names[j]]] != 0)
        current.scna.counts <- c(current.scna.counts, scnas)
        if (primary.list[i] == j){
            primary.counts <- c(primary.counts, muts)
            primary.scna <- c(primary.scna, scnas)
        }else{
            recurrent.counts <- c(recurrent.counts, muts)
            recurrent.scna <- c(recurrent.scna, scnas)
            
        }
    }
    pair.mut.counts[[i]] <- current.counts
    pair.scna.counts[[i]] <- current.scna.counts
}

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
for (i in 1:length(pair.scna.counts)){
    stdv <- sd(pair.scna.counts[[i]])
    all.stdv <- c(all.stdv, stdv / mean(pair.mut.counts[[i]]))
}
all.stdv

mean(all.stdv)
sd(unlist(pair.scna.counts)) / mean(unlist(pair.scna.counts))

t.test(primary.counts, recurrent.counts)
t.test(primary.scna, recurrent.scna)



## look at tumors where recurrences were of diffrent grades


## plot evolution of allelic fractions over time
# Basic line graph with points: stolen from Rcookbook
ggplot(data=dat1, aes(x=time, y=total_bill, group=sex, shape=sex)) + 
    geom_line(size=1.5) + 
    geom_point(size=3, fill="white") +
    scale_shape_manual(values=c(22,21))

## takes a list of genes, and returns their allelic fraction for all samples in given sublist
samples <- input.list[[2]]
# gene.lists <- list(c("SRP13", "CRIPAK", "LRP1B"))
# auto.gene.lists see above
fractions <- c()
recurrence <- c()
gene <- c()
of.interest <- c()
for (i in 1:length(samples)){
    pair.name <- master.table[master.table$Tumor.Name == samples[i], ]$Pair.Name[1]
    mutations <- FilterMaf(disc.snindels.duplicates, pair.name,"Tumor_Sample_Barcode")
    mutations <- mutations[mutations$Variant_Classification != "Silent", ]
    mutations <- PerSampleMaf(mutations, "Hugo_Symbol")
    average <- mean(mutations$i_tumor_f)
    mutations <- mutations[mutations$Hugo_Symbol %in% auto.gene.lists[[2]], ]
    fractions <- c(fractions, mutations$i_tumor_f, average)
    recurrence <- c(recurrence, rep(pair.name, nrow(mutations) + 1))
    of.interest <- c(of.interest, rep("genes", nrow(mutations)), "average")
    
    gene <- c(gene, mutations$Hugo_Symbol, "Average")
    if ("NF2" %in% mutations$Hugo_Symbol){
        idx <- max(which(gene %in% "NF2"))
        of.interest[idx] <- "NF2"
    }
}
time.series <- data.frame(fractions, recurrence, gene, of.interest)
# Map sex to color
x_levels <- unique(recurrence)
time.series$recurrence <- factor(time.series$recurrence, levels = x_levels)
ggplot(data=time.series, aes(x=recurrence, y=fractions, group=gene, colour=gene)) +
    geom_line(size=1.5, aes(linetype=of.interest)) +
    geom_point(size=2.5) +
    xlab("Recurrence") +
    ylab("Allelic Fraction") +
    scale_linetype_manual(values=c("dotted", "solid", "longdash"))


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



## compare average number of mutations shared between any two biopsies from same patient from maf that has already been forcecalled
PatientAverage <- function(maf, filler){
    ## of columns before data starts
    percent.shared <- c()
    for (i in (filler + 1):(ncol(maf) - 1)){
        for (j in (i + 1):ncol(maf)){
            shared <- rowSums(maf[, c(i,j)]) == 2
            not.shared <- rowSums(maf[, c(i,j)]) == 1
            percent.shared <- c(percent.shared, shared / (shared + not.shared))
        }
    }
    return(percent.shared)
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
for (i in 1:length(study2.pats)){
    temp <- study2[study2$Sample == study2.pats[i], ]
    shared <- sum(rowSums(temp[, 5:6]) == 2)
    total <- nrow(temp)
    study2.percentages <- c(study2.percentages, shared/total)
}

study5 <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Heterogeneity/Previously published/mutations/Study5_per_patient_comparison_all_pts.txt", 
                     stringsAsFactors = F)
study5.pats <- unique(study5$Patient.ID)
study5.percentages <- c()
for (i in 1:length(study5.pats)){
    temp <- study5[study5$Patient.ID == study5.pats[i], ]
    shared <- sum(temp[, 3] == 1)
    total <- nrow(temp)
    study5.percentages <- c(study5.percentages, shared/total)
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

study.export <- data.frame(study2.percentages, study5.percentages, )



