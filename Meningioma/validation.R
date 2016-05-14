## Validation Analysis

## compute probability of seeing given number of muts based on priors of discovery rate

drop.list <- c("M1130-Tumor", "M1560-Tumor", "M46-Tumor")

## M460 > M46

disc.count <- 3
disc.size <- 37
val.count <-11
val.size <- 69
prob <- 1 - pbinom(val.count - 1, val.size, disc.count / disc.size)




## figures out which validation genes to get rid of
gene.list <- genelist.generator[genelist.generator$i_tumor_f > .094, ]
val.snindels <- val.snindels[val.snindels$i_tumor_f > .09, ]


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
val.snindels <- PerSampleMaf(val.snindels, "Hugo_Symbol")
val.snindels.2 <- ReccurentMaf(val.snindels, "Hugo_Symbol")

validation.list <- validation.list[order(validation.list[,1]), ]
combined <- validation.list
combined <- cbind(combined, combined$Gene %in% meningioma.list[,1], combined$Gene %in% schwan.list[,1], 
                  combined$Gene %in% gene.list$Hugo_Symbol, combined$Gene %in% gene.list.2$Hugo_Symbol, combined$Gene %in% val.snindels$Hugo_Symbol, 
                  combined$Gene %in% val.snindels.2$Hugo_Symbol)
colnames(combined) <- c("Gene", "Notes", "In.Meningioma.List", "In.Schwan.List", "Mutated.In.Disc", "Mutated.2.In.Discovery", 
                        "Mutated.in.val", "Mutated.2.in.val")

combined <- combined[, c(1:3, 6,8,5,7,4)]

combined <- combined[order(combined[,3], combined[,4], combined[,5], combined[,6], combined[, 7]), ]

combined <- cbind(combined, combined[,3] == F & combined[,4] == F & combined[,5] == F & combined[,6] == F & combined[, 7] == F)
colnames(combined)[9] <- "drop"
 
write.csv(combined, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Mutations/decision.time.csv", row.names = F)



