## Cals for Oncopanel Analysis
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

## Gets CSV data, formats and cleans it
onc.data <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/R_input.txt', stringsAsFactors = FALSE, comment.char = "", 
                   header = T, strip.white = T)
onc.data <- onc.data[1:sum(onc.data$Last != ""), ]

onc.data <- onc.data[onc.data$Oncopanel. == 1, ]


by.gene <- NULL

## loops through databse, extracts gene level information

for (i in 1:nrow(onc.data)){
    if (onc.data$Mutation...DNA.Variants[i] == "" | onc.data$Mutation...DNA.Variants[i] == "no results"){
        # do nothing
    }else{
        lst <- strsplit(onc.data$Mutation...DNA.Variants[i], ";")[[1]]
        for (j in 1:length(lst)){
            split1 <- strsplit(lst[j], "c.")[[1]]
            name = trimws(split1[1])
            split2 <- strsplit(split1[2], " ")
            aa <- split2[[1]][1]
            protien <- split2[[1]][2]
            remainder <- strsplit(split1[2], ",")[[1]][2]
            temp <- onc.data[i, c(1:3, 23, 33,35:38)]
            temp1 <- data.frame(name, aa, protien, remainder, stringsAsFactors = F)
            temp <- cbind(temp, temp1)
            by.gene <- rbind(by.gene, temp)
        }
    }
}

colnames(by.gene)[10:13] <- c("gene", "amino.acid", "protein", "rest")

## gets total number of occurences for each gene
for (i in 1:nrow(by.gene)){
    count = sum(by.gene[,10 ] == by.gene[i, 10])
    by.gene[i, 14] <- count
}
write.table(by.gene, "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/by_gene.txt", row.names = F, quote = F, sep = "\t")
## last modified 9/23/16


## Analysis

# Check to see if disrupted or chr loss correlates with subtype

adenomas <- onc.data[onc.data$Pathology == "Pituitary Adenoma", ]
by.gene.patient <- PerSampleMaf(by.gene, "gene", "BWH_MRN")
by.gene.adenomas <- FilterMaf(by.gene.patient, c("Functional", "Null", "Gonadotrope"), "pathology.anatomic")


# correlation of 1, 11, and disruption with subtype
yvals <- c("Chr.1.Loss", "Chr.11.loss", "disruption")
mtrx1 <- matrix(NA, 1, 3)
for (i in 1:3){
    mini <- CleanYourRoom(adenomas[, c("pathology.anatomic", yvals[i])])
    tbl <- table(mini[, 1], mini[, 2])
    math.party <- fisher.test(tbl)
    mtrx1[1, i] <- signif(math.party$p.value, 2)
}
mtrx1

## correlation between disruption and other clinical features
table(onc.data$Atypical., onc.data$disruption)
table(onc.data$Recurrent.Tumor., onc.data$disruption)

## correlation of mutational load with clinical characteristics
colums <-c("BWH_MRN", "DOB", "Pathology.Subtype", "Sex", "Recurrent.Tumor.", "disruption", "MIB....", "pathology.clinical")
patient.stats <- NULL
mrns <- unique(by.gene.adenomas$BWH_MRN)
for (i in 1:length(mrns)){
    if (i == 1){
        patient.stats <- onc.data[onc.data$BWH_MRN == mrns[i], colums]
        patient.stats$muts <- sum(by.gene.adenomas$BWH_MRN == mrns[1])
    }else{
    muts <- sum(by.gene.adenomas$BWH_MRN == mrns[i])
    patient.stats <- rbind(patient.stats, cbind(onc.data[onc.data$BWH_MRN == mrns[i], colums], muts))
    }
}

### convert to testable format
numeric <- as.numeric(patient.stats$MIB....)
patient.stats$MIB_Low <- is.na(numeric) | numeric < 3
patient.stats$DOB <- mdy(patient.stats$DOB)

today <- mdy("09-01-2016")
x <- as.period(today- patient.stats$DOB)
patient.stats$Age <- as.numeric(x, "years")

p.values <- c()

p.values <- c(p.values, t.test(patient.stats[patient.stats$Recurrent.Tumor. == 1, ]$muts, patient.stats[patient.stats$Recurrent.Tumor. == 0, ]$muts )$p.value)

p.values <- c(p.values, t.test(patient.stats[patient.stats$Sex == "F" | patient.stats$Sex == "Female", ]$muts, patient.stats[patient.stats$Sex == "M" 
                                                                                                                             | patient.stats$Sex == "Male", ]$muts )$p.value)

p.values <- c(p.values, t.test(patient.stats[patient.stats$disruption == 1, ]$muts, patient.stats[patient.stats$disruption == 0, ]$muts )$p.value)

p.values <- c(p.values, t.test(patient.stats[patient.stats$MIB_Low == T, ]$muts, patient.stats[patient.stats$MIB_Low == F, ]$muts )$p.value)

p.values <- c(p.values, t.test(patient.stats[patient.stats$pathology.clinical == "Nonfunctional", ]$muts, patient.stats[patient.stats$pathology.clinical == "Functional", ]$muts )$p.value)



p.values <- c(p.values, t.test(patient.stats[patient.stats$Age == T, ]$muts, patient.stats[patient.stats$Age == F, ]$muts )$p.value)


## Chr loss freq

up.data <- adenomas[adenomas$Up != "", ]
down.data <- adenomas[adenomas$Down != "", ]
up.chr <- up.data$Up
down.chr <- down.data$Down
down.list <- strsplit(down.chr, ",")
up.list <- strsplit(up.chr, ",")

up <- c()
for (i in 1:length(up.list)){
    up <- c(up, up.list[[i]])
}

down <- c()
for (i in 1:length(down.list)){
    down <- c(down, down.list[[i]])
}

up <- as.numeric(up)
down <- as.numeric(down)

down <- down[!is.na(down)]


hist(up, breaks = 1:22)
hist(down, breaks = 1:23)

## mutations


## cohort statistics
## generate maf for functional and nonfunctional adenomas

func.adenomas <- by.gene.adenomas[by.gene.patient$pathology.clinical == "Functional", ]
nonfunc.adenomas <- by.gene.adenomas[by.gene.patient$pathology.clinical == "Nonfunctional", ]

# test difference in mutation rates
func.table <- table(func.adenomas$gene)
nonfunc.table <- table(nonfunc.adenomas$gene)

wilcox.test(mutations.nonfunc, mutations.func)

## Gene level analysis

## generate plot of mutations in all tumors
setwd("C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/figures + tables/")
pdf("mutation barplot all tumors.pdf", width = 10, height = 7)
PlotMaf(by.gene.patient, "gene", 18, title = "Genes Mutated in at least 4 pituitary tumor patients")
dev.off()

## generate plot of adenoma mutations
pdf("mutation barplot adenomas.pdf", width = 10, height = 7)

PlotMaf(by.gene.adenomas, "gene", 25, title = "Genes Mutated in at least 4 adenoma patients")
dev.off()


nrow(by.gene.patient) / 106

length(unique(by.gene.patient$gene))



## plot mutations in functional vs nonfunctional
## Comparison of functional vs nonfunctional tumors
pdf("mutations functional vs nonfunctional.pdf", width = 15, height = 7)
combined.table <- EqualizeTable(table(func.adenomas$gene),table(nonfunc.adenomas$gene))
barplot(combined.table, main = "Comparison of functional and nonfunctional adenomas", 
        legend.text = c("Functional adenomas", "Nonfunctional Adenomas"), las = 2, beside = TRUE, args.legend = list(x = "topleft"))
dev.off()
adenoma.table <- table(by.gene.adenomas$gene)
table(by.gene.adenomas$Pathology.Subtype)
gh.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "GH", ]$gene)
gh.df <- data.frame(gh.table, rep("GH", length(gh.table)))
fsh.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "FSH", ]$gene)
fsh.df <- data.frame(fsh.table, rep("FSH", length(fsh.table)))
acth.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "ACTH", ]$gene)
acth.df <- data.frame(acth.table, rep("ACTH", length(acth.table)))
prl.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "PRL", ]$gene)
prl.df <- data.frame(prl.table, rep("PRL", length(prl.table)))
null.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype %in% c("Null", "null"), ]$gene)
null.df <- data.frame(null.table, rep("Null", length(null.table)))
colnames(gh.df) <- c("Gene", "Freq", "Subtype")
colnames(fsh.df) <- c("Gene", "Freq", "Subtype")
colnames(acth.df) <- c("Gene", "Freq", "Subtype")
colnames(prl.df) <- c("Gene", "Freq", "Subtype")
colnames(null.df) <- c("Gene", "Freq", "Subtype")
subtype.df <- rbind(gh.df, fsh.df, acth.df, prl.df, null.df)
## creates data frame with per subtype rate
## 
## takes all genes with at least 3 mutations
by.gene.adenomas <- by.gene.adenomas[by.gene.adenomas$Pathology.Subtype != "hyperplasia", ]
inclusion.list <- dimnames(table(by.gene.adenomas$gene)[table(by.gene.adenomas$gene) > 3])[[1]]
overall.counts <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% inclusion.list,]$gene), decreasing = T)
subtype.df <- subtype.df[subtype.df$Gene %in% inclusion.list, ]
subtype.df$Gene <- factor(subtype.df$Gene, levels = names(overall.counts))
ggplot(subtype.df, aes(x = Gene, y = Freq, fill = Subtype)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks=(seq(2, 14, 2))) +
    scale_fill_grey() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


## reorder by pathway members
list1 <- c("BRCA1", "BRCA2", "PRKDC")
list2 <- c("ARID1B", "ARID1A", "ASXL1", "BRD4", "CUX1")
list3 <- c("ATM", "DEPDC", "GLI1", "GLI2", "GLI3", "GNAS", "NOTCH1", "NOTCH2", "NTRK1", "PTPRD", "PIK3CA", "TCF3", "TSC2")
others <- inclusion.list[!(inclusion.list %in% c(list1, list2, list3))]
counts1 <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% list1,]$gene), decreasing = T)
counts2 <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% list2,]$gene), decreasing = T)
counts3 <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% list3,]$gene), decreasing = T)
counts4 <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% list4,]$gene), decreasing = T)
counts5 <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% list5,]$gene), decreasing = T)
counts.other <- sort(table(by.gene.adenomas[by.gene.adenomas$gene %in% others,]$gene), decreasing = T)
counts <- c(counts1, counts2, counts3, counts.other)
subtype.df$Gene_subtype <- factor(subtype.df$Gene, levels = names(counts))
ggplot(subtype.df, aes(x = Gene_subtype, y = Freq, fill = Subtype)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks=(seq(2, 14, 2))) +
    scale_fill_grey() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))



## Hotspot mutations
hot <- ReccurentMaf(by.gene.patient, "amino.acid")
hot.table <- table(hot$gene)
hot.list <- names(hot.table[hot.table > 1])

PlotMaf(by.gene.adenomas[by.gene.adenomas$gene %in% hot.list, ], "gene", "Genes with hotspot mutations")
## GSEA





## genets
write.csv(adenoma.table[adenoma.table > 4], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/all_adenoma_5ormore.csv")
write.csv(adenoma.table[adenoma.table > 3], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/all_adenoma_4ormore.csv")
write.csv(adenoma.table[adenoma.table > 2], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/all_adenoma_3ormore.csv")
write.csv(adenoma.table[adenoma.table > 1], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/all_adenoma_2ormore.csv")

write.csv(func.table[func.table > 1], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/func_adenoma_2ormore.csv")
write.csv(func.table[func.table > 2], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/func_adenoma_3ormore.csv")
write.csv(func.table[func.table > 3], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/func_adenoma_4ormore.csv")

write.csv(nonfunc.table[nonfunc.table > 1], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/nonfunc_adenoma_2ormore.csv")
write.csv(nonfunc.table[nonfunc.table > 2], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/nonfunc_adenoma_3ormore.csv")
write.csv(nonfunc.table[nonfunc.table > 3], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/nonfunc_adenoma_4ormore.csv")

write.csv(gh.table[gh.table > 1], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/gh_adenoma_2ormore.csv")
write.csv(gh.table[gh.table > 2], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/gh_adenoma_3ormore.csv")
write.csv(gh.table[gh.table > 3], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/gh_adenoma_4ormore.csv")

write.csv(acth.table[acth.table > 1], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/acth_adenoma_2ormore.csv")
write.csv(acth.table[acth.table > 2], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/acth_adenoma_3ormore.csv")
write.csv(acth.table[acth.table > 3], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/acth_adenoma_4ormore.csv")

write.csv(prl.table[prl.table > 1], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/prl_adenoma_2ormore.csv")
write.csv(prl.table[prl.table > 2], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/prl_adenoma_3ormore.csv")
write.csv(prl.table[prl.table > 3], "C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Paper/GeNets/prl_adenoma_4ormore.csv")

