## Analysis of validation cohort

## compute probability of seeing given number of muts based on priors of discovery rate

## which mutations would we expected to see at given rate based on incidence in discovery cohort
disc.count <- 3
disc.size <- 37
val.count <-11
val.size <- 69
prob <- 1 - pbinom(val.count - 1, val.size, disc.count / disc.size)


## figures out which validation genes to get rid of. 
gene.list <- genelist.generator[genelist.generator$i_tumor_f > .094, ]


validation.list <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/validation_list.txt", stringsAsFactors = F)
meningioma.list <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/men_list.txt", stringsAsFactors = F)
schwan.list <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/schwan.list.txt"
                          , stringsAsFactors = F)

for (i in 1:nrow(gene.list)){
    cur <- gene.list$Tumor_Sample_Barcode[[i]]
    temp <- strsplit(cur, "G-")[[1]][1]
    temp <- strsplit(temp, "-")[[1]][1]
    gene.list$Tumor_Sample_Barcode[[i]] <- temp
}

gene.list <- PerSampleMaf(gene.list, "Hugo_Symbol", "Tumor_Sample_Barcode")
gene.list.2 <- ReccurentMaf(gene.list, "Hugo_Symbol")

validation.list <- validation.list[order(validation.list[,1]), ]
combined <- validation.list
combined <- cbind(combined, combined$Gene %in% meningioma.list[,1], combined$Gene %in% schwan.list[,1], 
                  combined$Gene %in% gene.list$Hugo_Symbol, combined$Gene %in% gene.list.2$Hugo_Symbol, combined$Gene %in% val.snindels$Hugo_Symbol, 
                  combined$Gene %in% val.snindels2$Hugo_Symbol)
colnames(combined) <- c("Gene", "Notes", "In.Meningioma.List", "In.Schwan.List", "Mutated.In.Disc", "Mutated.2.In.Discovery", 
                        "Mutated.in.val", "Mutated.2.in.val")

combined <- combined[, c(1:3, 6,8,5,7,4)]

combined <- combined[order(combined[,3], combined[,4], combined[,5], combined[,6], combined[, 7]), ]

combined <- cbind(combined, combined[,3] == F & combined[,4] == F & combined[,5] == F & combined[,6] == F & combined[, 7] == F)
colnames(combined)[9] <- "drop"
 
write.csv(combined, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/decision.time.csv", row.names = F)




sort(val.snindels[val.snindels$Hugo_Symbol == "NF2", ]$Tumor_Sample_Barcode)
val.snindels[val.snindels$Hugo_Symbol %in% gsea.validation.pi3k, ]

val.snindels[val.snindels$Tumor_Sample_Barcode %in% val.snindels[val.snindels$Hugo_Symbol %in% gsea.validation.pi3k, ]$Tumor_Sample_Barcode, ]

PlotMaf(val.snindels2, "Hugo_Symbol", 40, title = "Genes Mutated at least ")


val.snindels.clean <- val.snindels[val.snindels$i_tumor_f > .099, ]
val.snindels.clean <- ReccurentMaf(val.snindels.clean, "Hugo_Symbol")

write.csv(val.snindels.clean, "C:/Users/Noah/OneDrive/Work/Meningioma/GSEA/validation.mutated.2.csv")


## permutation testing for enrichment in validation sequenced genes
## define space of all possible genes
all.genes <- validation.list[,1]

## PIK3CA enrichment genes defined by all genes in pathway
cpdb.hits <- c("AKT1", "ERBB3", "EGFR", "NRG1", "ERBB2", "CDKN1A", "TSC2", "MTOR", 
               "FGFR3", "RICTOR", "PIK3CA")
path.hits <- sum(cpdb.hits %in% val.snindels.clean$Hugo_Symbol)
all.hits <- sum(all.genes %in% val.snindels.clean$Hugo_Symbol) - path.hits

## set up matrix
permuted <- matrix(NA, length(all.genes), 2)
colnames(permuted) <- c("pathway", "score")

## defines categories
permuted[,1] <- rep(c(1,0), c(length(cpdb.hits), length(all.genes) - length(cpdb.hits)))

## defines baseline categories
permuted[,2] <- rep(c(1,0,1,0), c(path.hits, length(cpdb.hits) - path.hits, all.hits, length(all.genes) - length(cpdb.hits) - all.hits))

all.means <- c()

for (i in 1:10000){
    vals <- sample(permuted[,2], nrow(permuted), replace = T)
    mean.diff <- mean(vals[permuted[,1] == 1]) - mean(vals[permuted[,1] == 0])
    all.means <- c(all.means, mean.diff)
}

hist(all.means)
