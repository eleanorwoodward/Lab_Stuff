## Pituitary WGS analysis
setwd("C:/Users/Noah/Syncplicity Folders/Pit (Linda Bi)/WGS project/")

## read in SVs, mutations, and indels
files <- list.files("Snowman/SVs/", "*.csv")

for (i in 1:length(files)){
    if (i == 1){
        pit.svs <- read.delim(paste("Snowman/SVs/", files[i], sep = ""), stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
    }else{
    temp <- read.delim(paste("Snowman/SVs/", files[i], sep = ""), stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
    pit.svs <- rbind(pit.svs, temp)
    }
}

files <- list.files("Snowman/Indels/", "*.annotated")

for (i in 1:length(files)){
    if (i == 1){
        pit.indels <- read.delim(paste("Snowman/Indels/", files[i], sep = ""), stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
        pit.indels$Tumor_Sample_Barcode <- files[i]
    }else{
        temp <- read.delim(paste("Snowman/Indels/", files[i], sep = ""), stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
        temp$Tumor_Sample_Barcode <- files[i]
        pit.indels <- rbind(pit.indels, temp)
    }
}


files <- list.files("MuTect/", "*.annotated")

for (i in 1:length(files)){
    if (i == 1){
        pit.snvs <- read.delim(paste("MuTect/", files[i], sep = ""), stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
    }else{
        temp <- read.delim(paste("MuTect/", files[i], sep = ""), stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
        pit.snvs <- rbind(pit.snvs, temp)
    }
}

genes <- pit.svs[, c("gene1", "sample")]
colnames(genes)[1] <- "gene2"
genes <- rbind(genes, pit.svs[, c("gene2", "sample")])

temp <- genes
temp$gene2 <- factor(temp$gene2, levels = names(sort(table(temp$gene2),decreasing = TRUE)[1:20]))
temp <- temp[!is.na(temp$gene2), ]
temp <- temp[temp$gene2 != "", ]
ggplot(data = temp, aes(x = gene2, fill = sample)) + geom_bar()

write.table(pit.svs, "Snowman/SVs/combined_SV_calls.txt", sep = "\t", row.names = FALSE)

## 1D analysis for mutations

## too many mutations

distance.1d <- 5000
pit.snvs.short <- pit.snvs[, c("Chromosome", "Tumor_Sample_Barcode", "Start_position")]
for (i in 1:nrow(pit.snvs.short)){
    pit.snvs.short$matches_1d[i] <- nrow(pit.snvs.short[pit.snvs.short$Tumor_Sample_Barcode != pit.snvs.short$Tumor_Sample_Barcode[i] & pit.snvs.short$Chromosome == pit.snvs.short$Chromosome[i] &
                                                   pit.snvs.short$Start_position > pit.snvs.short$Start_position[i] - distance.1d &
                                                       pit.snvs.short$Start_position < pit.snvs.short$Start_position[i] + distance.1d, ])
}


## restrict to smaller window
short.matches.idx <- which(pit.snvs.short$matches_1d != 0)
distance.1d.short <- 500
pit.snvs.short$matches_1d_short <- 0
for (i in short.matches.idx){
    pit.snvs.short$matches_1d_short[i] <- nrow(pit.snvs.short[pit.snvs.short$Tumor_Sample_Barcode != pit.snvs.short$Tumor_Sample_Barcode[i] & pit.snvs.short$Chromosome == pit.snvs.short$Chromosome[i] &
                                                      pit.snvs.short$Start_position > pit.snvs.short$Start_position[i] - distance.1d.short &
                                                          pit.snvs.short$Start_position < pit.snvs.short$Start_position[i] + distance.1d.short, ])
}

## look for hotspots
hot.matches.idx <- which(pit.snvs.short$matches_1d_short != 0)
pit.snvs.short$matches_1d_hot <- 0
for (i in hot.matches.idx){
    pit.snvs.short$matches_1d_hot[i] <- nrow(pit.snvs.short[pit.snvs.short$Tumor_Sample_Barcode != pit.snvs.short$Tumor_Sample_Barcode[i] & pit.snvs.short$Chromosome == pit.snvs.short$Chromosome[i] &
                                                            pit.snvs.short$Start_position == pit.snvs.short$Start_position[i], ])
}


## 1D analysis for indels

distance.1d <- 5000
pit.indels.short <- pit.indels[, c("Chromosome", "Tumor_Sample_Barcode", "Start_position", "Hugo_Symbol", "End_position", "Variant_Type", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele2")]
pit.indels.short <- pit.indels.short[!duplicated(pit.indels.short[, 2:4]), ]
for (i in 1:nrow(pit.indels.short)){
    pit.indels.short$matches_1d[i] <- nrow(pit.indels.short[pit.indels.short$Tumor_Sample_Barcode != pit.indels.short$Tumor_Sample_Barcode[i] & pit.indels.short$Chromosome == pit.indels.short$Chromosome[i] &
                                                      pit.indels.short$Start_position > pit.indels.short$Start_position[i] - distance.1d &
                                                          pit.indels.short$Start_position < pit.indels.short$Start_position[i] + distance.1d, ])
}


## restrict to smaller window
short.matches.idx <- which(pit.indels.short$matches_1d != 0)
distance.1d.short <- 500
pit.indels.short$matches_1d_short <- 0
for (i in short.matches.idx){
    pit.indels.short$matches_1d_short[i] <- nrow(pit.indels.short[pit.indels.short$Tumor_Sample_Barcode != pit.indels.short$Tumor_Sample_Barcode[i] &
                                                                      pit.indels.short$Chromosome == pit.indels.short$Chromosome[i] &
                                                                      pit.indels.short$Start_position > pit.indels.short$Start_position[i] - distance.1d.short & 
                                                                      pit.indels.short$Start_position < pit.indels.short$Start_position[i] + distance.1d.short, ])
}

## look for hotspots
hot.matches.idx <- which(pit.indels.short$matches_1d_short != 0)
pit.indels.short$matches_1d_hot <- 0
for (i in hot.matches.idx){
    pit.indels.short$matches_1d_hot[i] <- nrow(pit.indels.short[pit.indels.short$Tumor_Sample_Barcode != pit.indels.short$Tumor_Sample_Barcode[i] & 
                                                                    pit.indels.short$Chromosome == pit.indels.short$Chromosome[i] &
                                                          pit.indels.short$Start_position == pit.indels.short$Start_position[i], ])
}


## group indels for those that co-occured
annotate.idx <- which(pit.indels.short$matches_1d_hot != 0)
for (i in annotate.idx){
    pit.indels.short$matches_1d_hot[pit.indels.short$Start_position == pit.indels.short$Start_position[i]] <- paste("group", i, sep = "")
}

## group indels for those that were within 1 kb
annotate.idx <- which(pit.indels.short$matches_1d_short != 0)
for (i in annotate.idx){
    pit.indels.short$matches_1d_short[pit.indels.short$Start_position > pit.indels.short$Start_position[i] - distance.1d.short & 
                                          pit.indels.short$Start_position < pit.indels.short$Start_position[i] + distance.1d.short]  <- paste("group", i, sep = "")
}

## group indels for those that were within 10 kb
annotate.idx <- which(pit.indels.short$matches_1d != 0)
for (i in annotate.idx){
    pit.indels.short$matches_1d[pit.indels.short$Start_position > pit.indels.short$Start_position[i] - distance.1d & 
                                          pit.indels.short$Start_position < pit.indels.short$Start_position[i] + distance.1d]  <- paste("group", i, sep = "")
}



## 1D analysis for SVs

pit.svs.short <- pit.svs[, c("sample", "chr1", "pos1", "gene1", "gene1_100kb", "chr2", "pos2", "gene2", "gene2_100kb")]
temp <- pit.svs.short[, c(1, 6:9)]
pit.svs.short <- pit.svs.short[, -c(6:9)]
colnames(temp) <- colnames(pit.svs.short)
pit.svs.short <- rbind(pit.svs.short, temp)


distance.1d <- 5000
for (i in 1:nrow(pit.svs.short)){
    pit.svs.short$matches_1d[i] <- nrow(pit.svs.short[pit.svs.short$sample != pit.svs.short$sample[i] & pit.svs.short$chr1 == pit.svs.short$chr1[i] &
                                                                pit.svs.short$pos1 > pit.svs.short$pos1[i] - distance.1d &
                                                                pit.svs.short$pos1 < pit.svs.short$pos1[i] + distance.1d, ])
}


## group svs for those that were within 10 kb
annotate.idx <- which(pit.svs.short$matches_1d != 0)
for (i in annotate.idx){
    pit.svs.short$matches_1d[pit.svs.short$pos1 > pit.svs.short$pos1[i] - distance.1d & 
                                    pit.svs.short$pos1 < pit.svs.short$pos1[i] + distance.1d]  <- paste("group", i, sep = "")
}











