## Pituitary Tumors Analysis


## Iniatialize launch sequence
pancan_names <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "MEN", "OV", "PIT", "READ", "UCEC")
broad.list <- c()
focal.list <- c()
pit.focal.list <- c()
pit.broad.list <- c()
pit.disrupted <- c(1, 4, 9, 16, 22, 23, 26, 28, 32, 34, 37, 42)
pit.disrupted.list <- c("PIT001.Tumor", "PIT007.Tumor", "PIT013.Tumor", "PIT1002.Tumor","PIT1009.Tumor","PIT1010.Tumor","PIT1013.Tumor","PIT1015.Tumor",
                        "PIT1022.Tumor","PIT202.Tumor", "PIT311.Tumor", "PIT504.Tumor")

## Loop through pancan data set, calculation ratio of focal to broad disruption across samples
for (i in 1:length(pancan_names)){
broad.names <- list.files('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/broad/')
focal.names <- list.files('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/')
  
broad <- read.delim(paste('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/broad/', 
                          broad.names[i], sep = ""), stringsAsFactors = FALSE)
focal <- read.delim(paste('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/', 
                          focal.names[i], sep = ""), stringsAsFactors = FALSE)

ratio <- focal$focal_disruption_from_median / broad$broad_disruption_from_median

file.name <- paste('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/ratio/', 
                   pancan_names[i], '.csv', sep = "")
write.csv(ratio, file.name, row.names = FALSE)
}

## Loop through pancan data, calculate mean focal disruption for each tumor type
list <- c()
for (i in 1:length(pancan_names)){
    focal.names <- list.files('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/')
    
    focal <- read.delim(paste('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/', 
                              focal.names[i], sep = ""), stringsAsFactors = FALSE)
    
    list <- c(list, mean(focal$focal_disruption_from_median))
}

list

## Calculate ratio for pit samples

broad.pit <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/broad/pit_disruptionbroad_per_sample.151114.txt'
                      , stringsAsFactors = FALSE)
focal.pit <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/pit_disruptionfocal_per_sample.pt27.151114.txt' 
                          ,stringsAsFactors = FALSE)
ratio.pit <- focal.pit$focal_disruption_from_median / broad.pit$broad_disruption_from_median


# Broad alterations in non-disrupted subset of tumors
broad.pit.quiet <- FilterMaf(broad.pit, pit.disrupted.list, "sample_id", F)

focal.pit.disrupted <- FilterMaf(focal.pit, pit.disrupted.list, "sample_id", T)

## Focal alterations in disrupted subtype
mean(focal.pit.disrupted$focal_disruption_from_median) / mean(list)



mean(broad.pit.quiet$broad_disruption_from_mean)


## Allelic fraction information per tumor
setwd("C:/Users/Noah/OneDrive/Work/Coding/R/output/")
x <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/Mutation-Indels/Pit20140506.final_analysis_set.maf.txt', 
                stringsAsFactors = FALSE)
x <- FilterMaf(x, "Silent", "Variant_Classification", FALSE)

sample.list <- unique(x$Tumor_Sample_Barcode)
for (i in 1:length(sample.list)){
  fig.name <- sample.list[i]
    pdf(paste(fig.name,"_allelic_fractions", ".pdf", sep = ""))
    temp <- FilterMaf(x, fig.name, "Tumor_Sample_Barcode")
    genelist <- unique(temp$Hugo_Symbol)
    avgs <- AverageMaf(temp, genelist, "Hugo_Symbol","i_tumor_f")
    barplot(avgs, names.arg = genelist, main = paste("Average Allelic Fraction of Mutations in sample ", fig.name
            , sep = ""), las = 2)
    abline(h = mean(avgs))
    dev.off()
    
}

disrupted.maf <- FilterMaf(x, pit.disrupted.list, "Tumor_Sample_Barcode")

quiet.maf <- FilterMaf(x, pit.disrupted.list, "Tumor_Sample_Barcode", FALSE)

disrupted.filter <- disrupted.maf[disrupted.maf$Tumor_Sample_Barcode != "PIT504-Tumor", ]
t.test(table(disrupted.filter$Tumor_Sample_Barcode),  table(quiet.maf$Tumor_Sample_Barcode))

## Mutation Analysis
x.recurrent <- ReccurentMaf(x, "Hugo_Symbol")
PlotMaf(x.recurrent, "Hugo_Symbol")

## Plot relationship between mutation rate and power to detect
power <- c()
power.3 <- c()
power.med <- c()
power.med.3 <- c()
mut.rate <- seq(.0075, .20, .01)

for(i in 1:length(mut.rate)){
power <- c(power, pbinom(2, 41, mut.rate[i], FALSE))
power.3 <- c(power.3, pbinom(3, 42, mut.rate[i], FALSE))
power.med <- c(power.med, pbinom(2, 32, mut.rate[i], FALSE))
power.med.3 <- c(power.med.3, pbinom(3, 32, mut.rate[i], FALSE))
}

plot(mut.rate, power, xlab = "Mutation Rate", ylab = "Power to Detect Mutations", 
     main = "Figure 5: Power Calculation", pch = 16)
points(mut.rate, power.3, pch = 15)
points(mut.rate, power.med, pch = 1)
points(mut.rate, power.med.3, pch = 16)

  
legend("topleft", c("Power for at least 3 mutations", "Power for at least 4 mutations", 
                    "At least 3 mutations, 10 bad samples", 
                    "At least 4 mutations, 10 bad samples"), pch =c(0, 15, 1, 16))


## disrupted vs non-disrupted in validation

validation.data <- read.delim("C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/data/R_input.txt", stringsAsFactors = F)
validation.data <- validation.data[validation.data$Oncopanel. == 1, ]
validation.data <- validation.data[validation.data$Pathology == "Pituitary Adenoma", ]
pituitary.comut <- validation.data[, c("pathology.anatomic","disruption", "Chr.1.Loss")]
write.csv(pituitary.comut, "C:/Users/Noah/Dropbox/Work/Pits/ClinCanRes resubmission/Figs/fig1.validation.csv", row.names = F)
pituitary.comut <- pituitary.comut[order(pituitary.comut$pathology.anatomic, pituitary.comut$disruption, pituitary.comut$Chr.1.Loss), ]

table(validation.data$pathology.anatomic == "Null", validation.data$disruption)
table(validation.data$Atypical. == 1, validation.data$disruption)



## Reviews
## RNA expression

rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GDS4275_full.soft", skip = 100)
rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GDS4859_full.soft", skip = 69, nrows = 200)
rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GSE51618_family.soft", skip = 90, nrows = 200)
rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GSE22812_family.soft", skip = 108, nrows = 200)
rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GSE22812_series_matrix.txt", skip = 71, nrows = 200)
rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GSE22812_family.xml.tar", skip = 1, nrows = 200)
rna.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/RNA/GSE46311_family.soft", skip = 116, nrows = 200)


roseta <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE22812/downloaded_key.txt", stringsAsFactors = F)
upload.file <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE22812/CodeLink_Human_Whole_Genome.chip.txt", stringsAsFactors = F)


for (i in 1:nrow(upload.file)){
    oldval <- upload.file$Probe.Set.ID[i]
    upload.file[i, 1] <- roseta[roseta$PROBE_NAME == oldval, 1]
}

write.csv(upload.file, "C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE22812/fixed.key.csv")

## remove unexpressed genes: 20% or more below threshold value. 
precollapse.sample3 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE22812/GSE22812null.gct",
                      stringsAsFactors = F, skip = 2)

precollapse.sample3 <- precollapse.sample3[!is.na(rowSums(precollapse.sample3[, -(1:2)])), ]

bye.bye <- rep(F, nrow(precollapse.sample3))

for (i in 1:nrow(precollapse.sample3)){
    low <- (precollapse.sample3[i, -(1:2)] < 5)
    total <- sum(low)
    if (total > 2){
        bye.bye[i] <- T
    }else{
        idx <- which(low)
        precollapse.sample3[i, idx + 2] <- 5
    }
}

precollapse.sample3 <- precollapse.sample3[bye.bye, ]


write.table(precollapse.sample3, "C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE22812/GSE22812null.fixed.gct", 
            row.names = F, sep = "\t", quote = F)

## load samples
sample1 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/GDS4275null.collapsed.gct",
                      stringsAsFactors = F, skip = 2)

#sample2 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/GDS4859null.collapsed.gct", stringsAsFactors = F, skip = 2)

sample3 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/GSE22812null.collapsed.gct",
                      stringsAsFactors = F, skip = 2)

sample.expressed.3 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE22812/GSE22812null.fixed.collapsed.gct",
                      stringsAsFactors = F, skip = 2)

sample4 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/GSE51618null.collapsed.gct",
                      stringsAsFactors = F, skip = 2)

sample5 <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/Clustering/GSE46311_expression_for_CIN_analysis.txt", stringsAsFactors = F)
sample5 <- sample5[, -c(1,3)]
colnames(sample5)[1] <- "Name"

## remove samples that have all NAs
sample3 <- sample3[!is.na(rowSums(sample3[, -(1:2)])), ]
sample4 <- sample4[!(rowSums(sample4[, -(1:2)], na.rm = T) == 0), ]
sample1 <- sample1[!(rowSums(sample1[, -(1:2)], na.rm = T) == 0), ]
sample5 <- sample5[!(rowSums(sample5[, -(1:2)], na.rm = T) == 0), ]
sample.list <- list(sample4, sample3, sample1, sample5)
sample.list <- list(sample4, sample3, sample1)
combat.ready <- NA

for (i in 1:length(sample.list)){
    print(paste("analyzing sample # ", i))
    temp <- sample.list[[i]]
    if (i == 1){
        ## sets up master list with all genes found in othe samples
        combat.ready <- temp
        combat.ready <- combat.ready[combat.ready$Name %in% sample3$Name, ]
        combat.ready <- combat.ready[combat.ready$Name %in% sample1$Name, ]
        #combat.ready <- combat.ready[combat.ready$Name %in% sample5$Name, ]
        
    }else{
        # populate column with zero as default, to be modified with value of expression in subsequent arrays
        temp <- temp[temp$Name %in% sample4$Name, ]
        combat.ready[, (ncol(combat.ready) + 1):(ncol(combat.ready) + ncol(temp) - 2)] <- NA
        current.empties <- (ncol(combat.ready) - (ncol(temp) - 3)):ncol(combat.ready)
        
        ## Loops through sample maf, checking each individual row
        bad.boys <- c()
        for (j in 1:nrow(temp)){
            if(j %% 100 == 0){
                print(paste("Analyzing row ", j, " sample ", i))
            }
           
            ## checks if region is already in list
            hits <- combat.ready$Name == temp$Name[j] 
            if (sum(hits) > 0){
                ## if present, add 1 in column of all overlaps
                row <- which(hits)
                combat.ready[row, current.empties] <- temp[j, -(1:2)]
            }else{
                ## if not, do nothing, as combat will fail
            }
            
        }
        colnames(combat.ready)[current.empties] <- colnames(temp)[-(1:2)]
    }
}

## redo combat values
combat.ready <- combat.ready[rowSums(combat.ready[, -(1:2)]) > - 2, ]
combat.ready[, 1] <- 1:nrow(combat.ready)


for (i in 1:46){
    if (i < 10){
        testing[, i + 2] <- rnorm(5000, 10, 1.5)
    }else if (i < 20){
        testing[, i + 2] <- rnorm(5000, 7, 2.5)
    }else if (i < 30){
        testing[, i + 2] <- rnorm(5000, 12, 0.5)
    }else if (i < 40){
        testing[, i + 2] <- rnorm(5000, 2, 4.5)
    }else{
        testing[, i + 2] <- rnorm(5000, 6, 6)
    }
}
write.table(combat.ready, "C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/combat.ready.small.tsv", row.names = F,
            quote = F, sep = "\t")

expression.combat <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/combat.ready.combat.gct", stringsAsFactors = F, 
                                skip = 2, header = T)

combat.disrupt <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/ComBat/sample_name_key.txt", stringsAsFactors = F)
## Copy number data

x <- read.delim("C:/Users/Noah/OneDrive/Work/Coding/R/dbs/all_genes_hg38.txt")
bands <- read.delim("C:/Users/Noah/Documents/Big Files/dbs/c1.all.v5.1.symbols.gmt", header = FALSE, stringsAsFactors = F)
bands <- cbind(bands[, 1:2], bands)
for (i in 1:nrow(bands)){
    vals <- strsplit(bands[i,1], split = "p")
    if (length(vals[[1]]) == 2){
        bands[i, 2] <- vals[[1]][1]
        bands[i, 3] <- "p"
        bands[i, 4] <- vals[[1]][2]
    }else{
        vals <- strsplit(bands[i,1], split = "q")
        bands[i, 2] <- vals[[1]][1]
        bands[i, 3] <- "q"
        bands[i, 4] <- vals[[1]][2]
    }
}

bands <- bands[order(bands[, 2]), ]


## decide whether using expressed or all genes here
sample3 <- sample3[rowSums(sample3[, -(1:2)]) > - 2, ]
sample.expressed.3 <- sample.expressed.3[rowSums(sample.expressed.3[, -c(1:2)]) > -2, ]
sample5 <- sample5[rowSums(sample5[, -(1:2)]) > - 2, ]
expression <- sample.expressed.3[rowSums(sample.expressed.3[, -(1:2)]) > - 2, ]
expression <- sample3[rowSums(sample3[, -(1:2)]) > - 2, ]

colnames(expression) <- c("name", "description", "sample1-d", "sample6-q", "sample3-d", "sample11-q", "sample8-q", "sample4-d", "sample7-q",
                          "sample2-d",  "sample12-q", "sample13-d", "sample10-d", "sample9-d", "sample5-d" )

expression <- expression[, order(colnames(expression))]
expression <- cbind(expression[, 1:3], expression[, 8:15], expression[, 4:7])
expression[, 1] <- 1:8383
write.table(expression, "C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/Clustering/GSE22182.cleaned.txt", 
            row.names = F, sep = "\t", quote = F)

##5, 5, 5p
# 8, 8q, 8q,8
# 1q, 1q, 1q
# 9, 9 , 9, 9
# 7 7 7p
# 3 3 3p

## losses
# 11p 11p 11 11 11
# 15q

# 14q, 14q, #19p, 19p


## 11p loss: samples 1-5
bands.11p <- bands[bands[, 2] == "chr11" & bands[, 3] == "p", ]
genes.11p <- unlist(bands.11p[, -(1:4)])
expression.11p <- expression[expression$name %in% genes.11p, ]

## 3p gain: 1, 3 10
bands.3p <- bands[bands[, 2] == "chr3" & bands[, 3] == "p", ]
genes.3p <- unlist(bands.3p[, -(1:4)])
expression.3p <- expression[expression$name %in% genes.3p, ]

## 5p gain:
bands.5p <- bands[bands[, 2] == "chr5" & bands[, 3] == "p", ]
genes.5p <- unlist(bands.5p[, -(1:4)])
expression.5p <- expression[expression$name %in% genes.5p, ]

## 8q gain: 
bands.8q <- bands[bands[, 2] == "chr8" & bands[, 3] == "q", ]
genes.8q <- unlist(bands.8q[, -(1:4)])
expression.8q <- expression[expression$name %in% genes.8q, ]

## 1q gain: 
bands.1q <- bands[bands[, 2] == "chr1" & bands[, 3] == "q", ]
genes.1q <- unlist(bands.1q[, -(1:4)])
expression.1q <- expression[expression$name %in% genes.1q, ]

## 9p gain: 
bands.9p <- bands[bands[, 2] == "chr9" & bands[, 3] == "p", ]
genes.9p <- unlist(bands.9p[, -(1:4)])
expression.9p <- expression[expression$name %in% genes.9p, ]


## 9q gain: 
bands.9q <- bands[bands[, 2] == "chr9" & bands[, 3] == "q", ]
genes.9q <- unlist(bands.9q[, -(1:4)])
expression.9q <- expression[expression$name %in% genes.9q, ]

## 7p gain: 
bands.7p <- bands[bands[, 2] == "chr7" & bands[, 3] == "p", ]
genes.7p <- unlist(bands.7p[, -(1:4)])
expression.7p <- expression[expression$name %in% genes.7p, ]

## 15q loss: 
bands.15q <- bands[bands[, 2] == "chr15" & bands[, 3] == "q", ]
genes.15q <- unlist(bands.15q[, -(1:4)])
expression.15q <- expression[expression$name %in% genes.15q, ]

## log fold change 11p loss
fcs <- c()
for (i in 1:nrow(expression.11p)){
    fc <- mean(as.numeric(expression.11p[i, c(3:7)])) / mean(as.numeric(expression.11p[i, c(8:15)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)


## log fold change 3p gain
fcs <- c()
for (i in 1:nrow(expression.3p)){
    fc <- mean(as.numeric(expression.3p[i, c(3,5,12)])) / mean(as.numeric(expression.3p[i, c(4,6:11,13:15)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)


## log fold change 5p gain
fcs <- c()
for (i in 1:nrow(expression.5p)){
    fc <- mean(as.numeric(expression.5p[i, (c(1,4,10) + 2)])) / mean(as.numeric(expression.5p[i, -(c(-1, 0,1,4,10) + 2)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)

## log fold change 8q gain
fcs <- c()
for (i in 1:nrow(expression.8q)){
    fc <- mean(as.numeric(expression.8q[i, (c(1,3,5,11) + 2)])) / mean(as.numeric(expression.8q[i, -(c(-1, 0,1,3,5,11) + 2)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)

## log fold change 1q gain
fcs <- c()
for (i in 1:nrow(expression.1q)){
    fc <- mean(as.numeric(expression.1q[i, (c(3,4,5) + 2)])) / mean(as.numeric(expression.1q[i, -(c(-1, 0,3,4,5) + 2)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)


## log fold change 9p gain
fcs <- c()
for (i in 1:nrow(expression.9p)){
    fc <- mean(as.numeric(expression.9p[i, (c(3,9,10,12) + 2)])) / mean(as.numeric(expression.9p[i, -(c(-1,0, 3,9,10,12) + 2)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)

## log fold change 9q gain
fcs <- c()
for (i in 1:nrow(expression.9q)){
    fc <- mean(as.numeric(expression.9q[i, (c(3,9,10,12) + 2)])) / mean(as.numeric(expression.9q[i, -(c(-1,0, 3,9,10,12) + 2)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)

## log fold change 7p gain
fcs <- c()
for (i in 1:nrow(expression.7p)){
    fc <- mean(as.numeric(expression.7p[i, (c(9,10,13) + 2)])) / mean(as.numeric(expression.7p[i, -(c(-1,0,9,10,13) + 2)]))
    fcs <- c(fcs, log2(fc))
}
hist(fcs)
t.test(fcs)


## log fold change 15q loss
fcs <- c()
for (i in 1:nrow(expression.15q)){
    fc <- mean(as.numeric(expression.15q[i, (c(5,8,10) + 2)])) / mean(as.numeric(expression.15q[i, -(c(-1,0,5,8,10) + 2)]))
    fcs <- c(fcs, log(fc))
}
hist(fcs)
t.test(fcs)



## Cin70 genes
CIN70 <- c("TPX2", "PRC1", "FOXM1", "CDC2", "C20orf24", "TGIF2", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C", 
           "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", 
           "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "CKAP5", "NUP205", "CDC20", 
           "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", 
           "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "BRRN1", "MCM10", "LSM4", "LUZP5", "ASF1B", "ZWINT", 
           "TOPK", "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A", "GPI", "PDCD2L", 
           "SFRS2", "HDGF", "NXT1", "NEK2", "DHCR7", "STK6", "NDUFAB1", "KIAA0286", "KIF4A")
1x
2x
3x
4x
5x
6x
7x
8x
9x
10x
11x
12x
13
14
15x
16x
17x
18
19x
20x
21
22x

CIN70[CIN70 %in% expression.combat[, 1]]

expression.zed <- expression

expression.uri <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/GSE46311_expression_for_CIN_analysis.txt", stringsAsFactors = F)

expression.uri.zed <- expression.uri

for (i in 1:nrow(expression.zed)){
    avg <- rowSums(expression.zed[i, -(1:2)]) / 13
    expression.zed[i, -(1:2)] <- expression.zed[i, -(1:2)] - avg
    stdv <- sd(as.numeric(expression.zed[i, -(1:2)]))
    expression.zed[i, -(1:2)] <- expression.zed[i, -(1:2)] / stdv
}

for (i in 1:nrow(expression.uri.zed)){
    avg <- rowSums(expression.uri.zed[i, -(1:4)]) / 16
    expression.uri.zed[i, -(1:4)] <- expression.uri.zed[i, -(1:4)] - avg
    stdv <- sd(as.numeric(expression.uri.zed[i, -(1:4)]))
    expression.uri.zed[i, -(1:4)] <- expression.uri.zed[i, -(1:4)] / stdv
}

disrupted <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample9", "sample10", "sample13")
disrupted.uri <- colnames(expression.uri)[-(1:4)][c(F,T,F,F,T,F,T,T,T,T,T,F,T,T,F,F)]
nondisrupted.uri <- colnames(expression.uri)[-(1:4)][!c(F,T,F,F,T,F,T,T,T,T,T,F,T,T,F,F)]


cin.score <- colSums(expression.zed[expression.zed$name %in% CIN70, -(1:2)])
cin.score.uri <- colSums(expression.uri.zed[expression.uri.zed$Gene.symbol %in% CIN70, -(1:4)])
cin.score.uri[nondisrupted.uri]

## original dataset
cin.score[(names(cin.score) %in% disrupted)]
t.test(as.numeric(cin.score[(names(cin.score) %in% disrupted)]), as.numeric(cin.score[!(names(cin.score) %in% disrupted)]))
wilcox.test(as.numeric(cin.score[(names(cin.score) %in% disrupted)]), as.numeric(cin.score[!(names(cin.score) %in% disrupted)]))

## new dataset
t.test(as.numeric(cin.score.uri[disrupted.uri]), as.numeric(cin.score.uri[nondisrupted.uri]))
wilcox.test(as.numeric(cin.score.uri[disrupted.uri]), as.numeric(cin.score.uri[nondisrupted.uri]))

total.disrupted <- c(as.numeric(cin.score.uri[disrupted.uri]), as.numeric(cin.score[(names(cin.score) %in% disrupted)]))
total.quiet <- c(as.numeric(cin.score.uri[nondisrupted.uri]), as.numeric(cin.score[!(names(cin.score) %in% disrupted)]))


## combined
t.test(c(as.numeric(cin.score.uri[disrupted.uri]), as.numeric(cin.score[(names(cin.score) %in% disrupted)])),
      c(as.numeric(cin.score.uri[nondisrupted.uri]), as.numeric(cin.score[!(names(cin.score) %in% disrupted)])))
wilcox.test(c(as.numeric(cin.score.uri[disrupted.uri]), as.numeric(cin.score[(names(cin.score) %in% disrupted)])),
            c(as.numeric(cin.score.uri[nondisrupted.uri]), as.numeric(cin.score[!(names(cin.score) %in% disrupted)])))

ggdf <- as.data.frame(matrix(c(total.quiet, total.disrupted, rep(0, length(total.quiet)), rep(1, length(total.disrupted))), 29, 2))
colnames(ggdf) <- c("Score", "Functional")

ggplot(ggdf, aes(x = Functional, y = Score)) + geom_point()

## copy number data
copy.number.matrix <- matrix(0, 44, 13)
colnames(copy.number.matrix) <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8",
                                                          "sample9", "sample10", "sample11", "sample12", "sample13" )

rownames(copy.number.matrix) <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q", "7p", "7q",
                                  "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q", "13p", "13q", "14p",
                                  "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p",
                                  "20q", "21p", "21q", "22p", "22q")



copy.number.matrix[c("3p", "5p", "5q", "8p", "8q", "14q", "19p"), "sample1"] <- 1
copy.number.matrix[c("11p"), "sample1"] <- -1
copy.number.matrix[c("4q"), "sample2"] <- 1
copy.number.matrix[c("1p", "11p"), "sample2"] <- -1
copy.number.matrix[c("3p", "8q", "9p", "9q", "14q", "19p"), "sample3"] <- 1
copy.number.matrix[c("11p", "11q"), "sample3"] <- -1
copy.number.matrix[c("1q", "5p", "5q", "15q", "19q"), "sample4"] <- 1
copy.number.matrix[c("11p", "11q", "17p"), "sample4"] <- -1
copy.number.matrix[c("8q"), "sample5"] <- 1
copy.number.matrix[c("1p", "4p", "4q", "11p", "11q", "13q", "15q"), "sample5"] <- -1
copy.number.matrix[c("15q", "2p"), "sample8"] <- -1
copy.number.matrix[c("7p", "7q", "9p", "9q"), "sample9"] <- 1
copy.number.matrix[c("3p", "3q","5p","7p", "7q","9q", "9q"), "sample10"] <- 1
copy.number.matrix[c("15q"), "sample10"] <- -1
copy.number.matrix[c("8p", "8q"), "sample11"] <- 1
copy.number.matrix[c("9p", "9q"), "sample12"] <- 1
copy.number.matrix[c("7p", "20p", "20q"), "sample13"] <- 1
copy.number.matrix[c("13q"), "sample13"] <- -1

colnames(copy.number.matrix) <- c("sample1-d", "sample2-d", "sample3-d", "sample4-d", "sample5-d", "sample6-q", "sample7-q",
                                  "sample8-m", "sample9-m", "sample10-d", "sample11-q", "sample12-q", "sample13-d" )
colnames(expression.zed)[-(1:2)] <- c("sample1-d", "sample2-d", "sample3-d", "sample4-d", "sample5-d", "sample6-q", "sample7-q",
                                  "sample8-m", "sample9-m", "sample10-d", "sample11-q", "sample12-q", "sample13-d" )
copy.number.matrix <- t(copy.number.matrix)
signed.d <- dist(copy.number.matrix)
signed.clust <- hclust(signed.d)
plot(signed.clust, main = "copy number clustering-signed")

unsigned.d <- dist(abs(copy.number.matrix))
unsigned.clust <- hclust(unsigned.d)
plot(unsigned.clust, main = "copy number clustering-unsigned")

signed.d.gene <- dist(t(expression.zed[, -(1:2)]))
signed.clust.gene <- hclust(signed.d.gene)
plot(signed.clust.gene, main = "gene exprssion clustering- signed")

unsigned.d.gene <- dist(abs(t(expression.zed[, -(1:2)])))
unsigned.clust.gene <- hclust(unsigned.d.gene)
plot(unsigned.clust.gene, main = "gene expression clustering-unsigned")

all.samples.d <- dist(t(expression.combat[, -(1:2)]))
all.samples.clust <- hclust(all.samples.d)
plot(all.samples.clust, main = "gene expression clustering-unsigned")

## create gene expression matrix
gene.expression.matrix <- matrix(0, 44, 13)
colnames(gene.expression.matrix) <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8",
                                      "sample9", "sample10", "sample11", "sample12", "sample13" )

rownames(gene.expression.matrix) <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q", "7p", "7q",
                                      "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q", "13p", "13q", "14p",
                                      "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p",
                                      "20q", "21p", "21q", "22p", "22q")

colnames(expression.zed)[-(1:2)] <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8",
                                      "sample9", "sample10", "sample11", "sample12", "sample13" )

adios <- is.na(bands[, 4])
bands <- bands[!adios, ]

new.names <- unique(bands[, 2])[-(23:24)]

for (i in 1:length(new.names)){
    band.segments <- bands[bands[, 2] == new.names[i] & bands[, 3] == "p", ]
    if (nrow(band.segments) != 0){
        genes <- unlist(band.segments[, -(1:4)])
        expression.vals <- expression.zed[expression.zed$name %in% genes, ]
        gene.expression.matrix[i, ] <- colSums(expression.vals[-(1:2)]) / nrow(expression.vals)
        rownames(gene.expression.matrix)[i] <- paste(new.names[i], "p", sep = "")
    }
    
    band.segments <- bands[bands[, 2] == new.names[i] & bands[, 3] == "q", ]
    if (nrow(band.segments) != 0){
        genes <- unlist(band.segments[, -(1:4)])
        expression.vals <- expression.zed[expression.zed$name %in% genes, ]
        gene.expression.matrix[i, ] <- colSums(expression.vals[-(1:2)]) / nrow(expression.vals)
        rownames(gene.expression.matrix)[i] <- paste(new.names[i], "q", sep = "")
        
    }
    
}

colnames(gene.expression.matrix) <- c("sample1-d", "sample2-d", "sample3-d", "sample4-d", "sample5-d", "sample6-q", "sample7-q",
                                  "sample8-m", "sample9-m", "sample10-d", "sample11-q", "sample12-q", "sample13-d" )

## calculate distance matrix using signed (+,-) values
signed.d.gene.matrix <- dist(t(gene.expression.matrix[, -(1:2)]))
signed.clust.gene.matrix <- hclust(signed.d.gene.matrix)
plot(signed.clust.gene.matrix, main = "gene exprssion clustering- signed")

unsigned.d.gene.matrix <- dist(abs(t(gene.expression.matrix[, -(1:2)])))
unsigned.clust.gene.matrix <- hclust(unsigned.d.gene.matrix)
plot(unsigned.clust.gene.matrix, main = "gene expression clustering-unsigned")


## combat them, zscoring by row
combat.zed <- expression.combat
for (i in 1:nrow(combat.zed)){
    avg <- rowSums(combat.zed[i, -(1:2)]) / 46
    combat.zed[i, -(1:2)] <- combat.zed[i, -(1:2)] - avg
    stdv <- sd(as.numeric(combat.zed[i, -(1:2)]))
    combat.zed[i, -(1:2)] <- combat.zed[i, -(1:2)] / stdv
}

new.names <- read.delim("C:/Users/Noah/hubiC/Pituitary Tumor Paper/Revision/Clustering/combat_key.txt", header = F, stringsAsFactors = F)
new.names <- new.names[order(new.names[, 1], decreasing = F), ]
combat.zed <- combat.zed[, order(colnames(combat.zed), decreasing = F)]
cbind(new.names[-1, 1], colnames(combat.zed)[-c(1,48)])
colnames(combat.zed)[-c(1,48)] <- new.names[-1, 2]

unsigned.combat <- dist(abs(t(combat.zed[, -c(1,48)])))
unsigned.clust.combat <- hclust(unsigned.combat)
plot(unsigned.clust.combat, main = "gene expression clustering-unsigned")

cin.score <- colSums(combat.zed[combat.zed$Name %in% CIN70, -c(1,48)])



## calculating r squared
## average genes together for affected vs unaffected group, 

# consructs gene matrices for z-scored version
expression.zed.11p <- expression.zed[expression.zed$name %in% genes.11p, ]

expression.zed.3p <- expression.zed[expression.zed$name %in% genes.3p, ]

expression.zed.5p <- expression.zed[expression.zed$name %in% genes.5p, ]

expression.zed.8q <- expression.zed[expression.zed$name %in% genes.8q, ]

expression.zed.1q <- expression.zed[expression.zed$name %in% genes.1q, ]

expression.zed.9p <- expression.zed[expression.zed$name %in% genes.9p, ]

expression.zed.9q <- expression.zed[expression.zed$name %in% genes.9q, ]

expression.zed.7p <- expression.zed[expression.zed$name %in% genes.7p, ]

expression.zed.15q <- expression.zed[expression.zed$name %in% genes.15q, ]

## creats x and y value for regression; codes disrupted as whichever is greater in value

expression.values <- c()
disrupted.status <- c()

expression.values <- c(expression.values, rowSums(expression.zed.11p[, c(3:7)]), rowSums(expression.zed.11p[, c(8:15)]))
disrupted.status <- c(disrupted.status, rep(0, nrow(expression.zed.11p)), rep(1, nrow(expression.zed.11p)))

expression.values <- c(expression.values, rowSums(expression.zed.3p[, c(3,5,12)]), rowSums(expression.zed.3p[, c(4,6:11,13:15)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.3p)), rep(0, nrow(expression.zed.3p)))

expression.values <- c(expression.values, rowSums(expression.zed.5p[, (c(1,4,10) + 2)]), rowSums(expression.zed.5p[, -(c(-1, 0,1,4,10) + 2)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.5p)), rep(0, nrow(expression.zed.5p)))

expression.values <- c(expression.values, rowSums(expression.zed.8q[, (c(1,3,5,11) + 2)]), rowSums(expression.zed.8q[, -(c(-1, 0,1,3,5,11) + 2)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.8q)), rep(0, nrow(expression.zed.8q)))

expression.values <- c(expression.values, rowSums(expression.zed.1q[, (c(3,4,5) + 2)]), rowSums(expression.zed.1q[, -(c(-1, 0,3,4,5) + 2)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.1q)), rep(0, nrow(expression.zed.1q)))

expression.values <- c(expression.values, rowSums(expression.zed.9p[, (c(3,9,10,12) + 2)]), rowSums(expression.zed.9p[, -(c(-1,0, 3,9,10,12) + 2)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.9p)), rep(0, nrow(expression.zed.9p)))

expression.values <- c(expression.values, rowSums(expression.zed.9q[, (c(3,9,10,12) + 2)]), rowSums(expression.zed.9q[, -(c(-1,0, 3,9,10,12) + 2)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.9q)), rep(0, nrow(expression.zed.9q)))

expression.values <- c(expression.values, rowSums(expression.zed.7p[, (c(9,10,13) + 2)]), rowSums(expression.zed.7p[, -(c(-1,0,9,10,13) + 2)]))
disrupted.status <- c(disrupted.status, rep(1, nrow(expression.zed.7p)), rep(0, nrow(expression.zed.7p)))


## fit linear model on expression values with disrupted or not as dependent variable

trial.model <- lm(expression.values ~ disrupted.status)
summary(trial.model)

# 
# ## NMF clustering
# 
# example.matrix <- NULL
# 
# increasing.means <- rbind(matrix(rnorm(1000, 20,5), 100,10), matrix(rnorm(1000, 35, 5), 100,10), matrix(rnorm(100, 55, 5), 100,10))
# decreasing.means <- rbind(matrix(rnorm(1000, 45,5), 100,10), matrix(rnorm(1000, 26, 5), 100,10), matrix(rnorm(100, 20, 5), 100,10))
# neither.means <- rbind(matrix(rnorm(1000, 30,5), 100,10), matrix(rnorm(1000, 55, 5), 100,10), matrix(rnorm(100, 20, 5), 100,10))
# 
# example.matrix <- cbind(increasing.means, decreasing.means, neither.means)
# colnames(example.matrix) <- c(paste("inc", 1:10, sep = ""), paste("dec", 1:10, sep = ""), paste("neith", 1:10, sep = ""))
# 
# data(esGolub)
# esGolub
# esGolub <- esGolub[1:200, ]
# # remove the uneeded variable ✬Sample✬ from the phenotypic
# # data
# esGolub$Sample <- NULL
# # default NMF algorithm
# res <- nmf(esGolub, 3)
# res
# fit(res)
# V.hat <- fitted(res)
# dim(V.hat)
# summary(res)
# summary(res, target = esGolub)
# summary(res, class = esGolub$Cell)
# # get matrix W
# w <- basis(res)
# dim(w)
# # get matrix H
# h <- coef(res)
# dim(h)
# 
# res <- nmf(example.matrix, 3)
# 
# nmf.dist <- dist(t(coef(res)))
# nmf.clust <- hclust(nmf.dist)
# plot(nmf.clust, main = "gene expression clustering-unsigned")
# ## decided to use genepattern instead



# ##SOM analysis for clustering
# library(kohonen)
# 
# # Create a training data set (rows are samples, columns are variables
# # Here I am selecting a subset of my variables available in "data"
# data_train <- data[, c(2,4,5,8)]
# 
# # Change the data frame with training data to a matrix
# # Also center and scale all variables to give them equal importance during
# # the SOM training process. 
# data_train_matrix <- scale(t(example.matrix))
# data_train_matrix <- matrix(rnorm(5000, 20,5), 100,50)
# 
# # Create the SOM Grid - you generally have to specify the size of the 
# # training grid prior to training the SOM. Hexagonal and Circular 
# # topologies are possible
# som_grid <- somgrid(xdim = 10, ydim=10, topo="rectangular")
# 
# # Finally, train the SOM, options for the number of iterations,
# # the learning rates, and the neighbourhood are available
# 
# som_model <- som(data_train_matrix, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = T, n.hood = "circular")
# plot(som_model, type="changes")
# plot(som_model, type="count")
# 
# 
# 
# data(wines)
# set.seed(7)
# 
# training <- sample(nrow(wines), 120)
# Xtraining <- scale(wines[training, ])
# Xtest <- scale(wines[-training, ],
#                center = attr(Xtraining, "scaled:center"),
#                scale = attr(Xtraining, "scaled:scale"))
# 
# som.wines <- som(Xtraining, grid = somgrid(5, 5, "hexagonal"))
# 
# som.prediction <- predict(som.wines, newdata = Xtest,
#                           trainX = Xtraining,
#                           trainY = factor(wine.classes[training]))
# table(wine.classes[-training], som.prediction$prediction)
# str(wines)
# plot(som.wines, type = "count")
# 
# kohmap <- xyf(scale(wines), classvec2classmat(wine.classes),
#               grid = somgrid(5, 5, "hexagonal"), rlen=100)
# plot(kohmap, type="changes")
# plot(kohmap, type="codes", main = c("Codes X", "Codes Y"))
# plot(kohmap, type="counts")
# 
# ## palette suggested by Leo Lopes
# coolBlueHotRed <- function(n, alpha = 1) {
#     rainbow(n, end=4/6, alpha=alpha)[n:1]
# }
# plot(kohmap, type="quality", palette.name = coolBlueHotRed)
# plot(kohmap, type="mapping", 
#      labels = wine.classes, col = wine.classes+1,
#      main = "mapping plot")
# 
# ## add background colors to units according to their predicted class labels
# xyfpredictions <- classmat2classvec(predict(kohmap)$unit.predictions)
# bgcols <- c("gray", "pink", "lightgreen")
# plot(kohmap, type="mapping", col = wine.classes+1,
#      pchs = wine.classes, bgcol = bgcols[as.integer(xyfpredictions)], 
#      main = "another mapping plot")
# 
# ## Show 'component planes'
# set.seed(7)
# sommap <- som(scale(wines), grid = somgrid(6, 4, "hexagonal"))
# plot(sommap, type = "property", property = sommap$codes[,1],
#      main = colnames(sommap$codes)[1])
# 
# ## Another way to show clustering information
# plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")
# ## use hierarchical clustering to cluster the codebook vectors
# som.hc <- cutree(hclust(dist(sommap$codes)), 5)
# add.cluster.boundaries(sommap, som.hc)
# 
# ## and the same for rectangular maps
# set.seed(7)
# sommap <- som(scale(wines),grid = somgrid(6, 4, "rectangular"))
# plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")
# ## use hierarchical clustering to cluster the codebook vectors
# som.hc <- cutree(hclust(dist(sommap$codes)), 5)
# add.cluster.boundaries(sommap, som.hc)

sample3.zed <- sample3
## zscore the data
for (i in 1:nrow(sample3.zed)){
    avg <- rowSums(sample3.zed[i, -(1:2)]) / 13
    sample3.zed[i, -(1:2)] <- sample3.zed[i, -(1:2)] - avg
    stdv <- sd(as.numeric(sample3.zed[i, -(1:2)]))
    sample3.zed[i, -(1:2)] <- sample3.zed[i, -(1:2)] / stdv
}

sample3.zed.expressed <- sample.expressed.3
## zscore the data
for (i in 1:nrow(sample3.zed.expressed)){
    avg <- rowSums(sample3.zed.expressed[i, -(1:2)]) / 13
    sample3.zed.expressed[i, -(1:2)] <- sample3.zed.expressed[i, -(1:2)] - avg
    stdv <- sd(as.numeric(sample3.zed.expressed[i, -(1:2)]))
    sample3.zed.expressed[i, -(1:2)] <- sample3.zed.expressed[i, -(1:2)] / stdv
}

sample5.zed <- sample5
## zscore the data
for (i in 1:nrow(sample5.zed)){
    avg <- rowSums(sample5.zed[i, -(1:2)]) / 16
    sample5.zed[i, -(1:2)] <- sample5.zed[i, -(1:2)] - avg
    stdv <- sd(as.numeric(sample5.zed[i, -(1:2)]))
    sample5.zed[i, -(1:2)] <- sample5.zed[i, -(1:2)] / stdv
}




## euclidian distance
sample3.disruption <- combat.disrupt$disrupted[-(1:16)]
sample3.disruption[11] <- "Yes"
sample3.disruption.index <- sample3.disruption == "Yes"
sample5.disruption <- combat.disrupt$disrupted[1:16]
sample5.disruption.index <- sample5.disruption == "Yes"

colnames(sample3)[c(F,F,sample3.disruption.index)]

## comparing all samples to one another from same cohorts
sample3.disruption.distances <- dist(t(sample3.zed[, c(F, F, sample3.disruption.index)]), diag = F)
sample3.quiet.distances <- dist(t(sample3.zed[, c(F, F, !sample3.disruption.index)]), diag = F)

wilcox.test(sample3.disruption.distances, sample3.quiet.distances)


sample3.comparison.distances <- c()
# compare all samples from different classes
for (i in 1:sum(sample3.disruption.index)){
    idx <- which(sample3.disruption.index)[i]
    temp <- cbind(sample3.zed[, idx + 2], sample3.zed[, c(F, F, !sample3.disruption.index)])
    temp.d <- dist(t(temp), diag = F)
    sample3.comparison.distances <- c(sample3.comparison.distances, temp.d[1:sum(!sample3.disruption.index)])
}

wilcox.test(c(sample3.disruption.distances, sample3.quiet.distances), sample3.comparison.distances)



## same thing for expressed only versions of sample 3

## comparing all samples to one another from same cohorts
sample3.disruption.distances.exp <- dist(t(sample3.zed.expressed[, c(F, F, sample3.disruption.index)]), diag = F)
sample3.quiet.distances.exp <- dist(t(sample3.zed.expressed[, c(F, F, !sample3.disruption.index)]), diag = F)

wilcox.test(sample3.disruption.distances.exp, sample3.quiet.distances.exp)


sample3.comparison.distances.exp <- c()
# compare all samples from different classes
for (i in 1:sum(sample3.disruption.index)){
    idx <- which(sample3.disruption.index)[i]
    temp <- cbind(sample3.zed.expressed[, idx + 2], sample3.zed.expressed[, c(F, F, !sample3.disruption.index)])
    temp.d <- dist(t(temp), diag = F)
    sample3.comparison.distances.exp <- c(sample3.comparison.distances.exp, temp.d[1:sum(!sample3.disruption.index)])
}

wilcox.test(c(sample3.disruption.distances, sample3.quiet.distances), sample3.comparison.distances)


## same thing for sample5
sample5.disruption.distances <- dist(t(sample5.zed[, c(F, F, sample5.disruption.index)]), diag = F)
sample5.quiet.distances <- dist(t(sample5.zed[, c(F, F, !sample5.disruption.index)]), diag = F)

wilcox.test(sample5.disruption.distances, sample5.quiet.distances)


sample5.comparison.distances <- c()
# compare all samples from different classes
for (i in 1:sum(sample5.disruption.index)){
    idx <- which(sample5.disruption.index)[i]
    temp <- cbind(sample5.zed[, idx + 2], sample5.zed[, c(F, F, !sample5.disruption.index)])
    temp.d <- dist(t(temp), diag = F)
    sample5.comparison.distances <- c(sample5.comparison.distances, temp.d[1:sum(!sample5.disruption.index)])
}

t.test(c(sample5.disruption.distances, sample5.quiet.distances), sample5.comparison.distances)

distance.matrix.5 <- as.matrix(dist(t(sample5[, -(1:2)]), diag = F, upper = T))
t.test(unlist(distance.matrix.5[, sample5.disruption.index]), unlist(distance.matrix.5[, !sample5.disruption.index]))

distance.matrix.3 <- as.matrix(dist(t(sample3.zed[, -(1:2)]), diag = F, upper = T))
t.test(unlist(distance.matrix.3[, sample3.disruption.index]), unlist(distance.matrix.3[, !sample3.disruption.index]))

distance.matrix.3.exp <- as.matrix(dist(t(sample3.zed.expressed[, -(1:2)]), diag = F, upper = T))
t.test(as.vector(distance.matrix.3.exp[, sample3.disruption.index]), unlist(distance.matrix.3.exp[, !sample3.disruption.index]))

wilcox.test(c(as.vector(distance.matrix.3[, sample3.disruption.index]), as.vector(distance.matrix.5[, sample5.disruption.index])),
c(as.vector(distance.matrix.3[, !sample3.disruption.index]), as.vector(distance.matrix.5[, !sample5.disruption.index])))


# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name


d <- dist(t(sample5[, -(1:2)])) # euclidean distances between the rows
d <- dist(t(sample3.zed[, -(1:2)])) # euclidean distances between the rows

fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
#plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
 #    main="Metric	MDS",	type="n")
#text(x, y, labels = sample5.disruption, cex=.7)

disrupted <- sample5.disruption
qplot(x,y, color = disrupted, geom = "point", size = I(4))

ggframe1 <- data.frame(c(sample5.disruption.distances, sample5.quiet.distances), 
                      c(rep("disrupted", length(sample5.disruption.distances)), rep("quiet", length(sample5.quiet.distances))))
ggframe2 <- data.frame(c(sample3.disruption.distances, sample3.quiet.distances), 
                      c(rep("disrupted", length(sample3.disruption.distances)), rep("quiet", length(sample3.quiet.distances))))

ggframe1[,1] <- scale(ggframe1[,1])
ggframe2[,1] <- scale(ggframe2[,1])
colnames(ggframe1) <- c("Value", "Status")
colnames(ggframe2) <- c("Value", "Status")


ggframe <- rbind(ggframe1, ggframe2)

ggplot(ggframe, aes(x = Status, y = Value)) + geom_point()
t.test(ggframe[ggframe$Status == "disrupted", ]$Value, ggframe[ggframe$Status == "quiet", ]$Value)



## collapse data sets

table(sample5$Name %in% sample3$Name)

small.sample3 <- sample3[sample3$Name %in% sample5$Name, ]
small.sample5 <- sample5[sample5$Name %in% sample3$Name, ]
small.sample5 <- small.sample5[!duplicated(small.sample5$Name), ]
small.sample3 <- small.sample3[order(small.sample3$Name), ]
small.sample5 <- small.sample5[order(small.sample5$Name), ]
small.combined <- cbind(small.sample3[, -(1:2)], small.sample5[, -(1:2)])

d <- dist(t(small.combined)) # euclidean distances between the rows

fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

disrupted <- c(sample3.disruption, sample5.disruption)
qplot(x,y, color = disrupted, geom = "point", size = I(4))


## collapse data sets: z score version

table(sample5$Name %in% sample3$Name)

small.sample3 <- sample3.zed.expressed[sample3.zed.expressed$Name %in% sample5$Name, ]
small.sample5 <- sample5.zed[sample5$Name %in% sample3.zed.expressed$Name, ]
small.sample5 <- small.sample5[!duplicated(small.sample5$Name), ]
small.sample3 <- small.sample3[order(small.sample3$Name), ]
small.sample5 <- small.sample5[order(small.sample5$Name), ]
small.combined <- cbind(small.sample3[, -(1:2)], small.sample5[, -(1:2)])

d <- dist(t(small.combined)) # euclidean distances between the rows


fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
z <- fit$points[,3]

disrupted <- c(sample3.disruption, sample5.disruption)
disrupted.index <- disrupted == "Yes"
qplot(x,y, color = disrupted, geom = "point", size = I(4))
qplot(x,z, color = disrupted, geom = "point", size = I(4))
qplot(y,z, color = disrupted, geom = "point", size = I(4))



## check to make sure still significant
d.matrix <- as.matrix(d)
wilcox.test(d.matrix[, disrupted == "Yes"], d.matrix[, disrupted == "No"])


## redo dot plot with data that is from same distance calculation
small.combined.disruption.distances <- dist(t(small.combined[, disrupted.index]), diag = F)
small.combined.quiet.distances <- dist(t(small.combined[, !disrupted.index]), diag = F)


small.combined.comparison.distances <- c()
# compare all samples from different classes
for (i in 1:sum(disrupted.index)){
    idx <- which(disrupted.index)[i]
    temp <- cbind(small.combined[, idx], small.combined[ !disrupted.index])
    temp.d <- dist(t(temp), diag = F)
    small.combined.comparison.distances <- c(small.combined.comparison.distances, temp.d[1:sum(!disrupted.index)])
}

t.test(small.combined.disruption.distances, small.combined.quiet.distances)
t.test(small.combined.quiet.distances, small.combined.comparison.distances)

write.csv(ggframe3, "C:/Users/Noah/Dropbox/Work/Pits/ClinCanRes resubmission/Figs/drafts/figure_4_new.vals.csv")

ggframe3 <- data.frame(c(small.combined.disruption.distances,small.combined.quiet.distances, small.combined.comparison.distances), 
                       c(rep("disrupted", length(small.combined.disruption.distances)), rep("quiet", length(small.combined.quiet.distances)), 
                         rep("comparison", length(small.combined.comparison.distances))))

colnames(ggframe3) <- c("Value", "Status")


ggplot(ggframe3, aes(x = Status, y = Value)) + geom_point()
t.test(ggframe[ggframe$Status == "disrupted", ]$Value, ggframe[ggframe$Status == "quiet", ]$Value)


ggframe4 <- ggframe3
ggframe4[, 2] <- "All"
ggframe4 <- rbind(ggframe3, ggframe4)
ggplot(ggframe4, aes(x = Status, y = Value)) + geom_point()



## disrupted with disrupted, quiet with quiet, disrupted with quiet
