## Cals for Oncopanel Analysis
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

## Gets CSV data, formats and cleans it
onc.data <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/R_input.txt', stringsAsFactors = FALSE, comment.char = "", 
                   header = T, strip.white = T)
onc.data <- onc.data[1:sum(onc.data$Last != ""), ]

by.gene <- NULL

## loops through databse, extracts gene level information

cases <- onc.data

for (i in 1:nrow(cases)){
    if (cases$Mutation...DNA.Variants[i] == "" | cases$Mutation...DNA.Variants[i] == "no results"){
        # do nothing
    }else{
        lst <- strsplit(cases$Mutation...DNA.Variants[i], ";")[[1]]
        for (j in 1:length(lst)){
            split1 <- strsplit(lst[j], "c.")[[1]]
            name = trimws(split1[1])
            split2 <- strsplit(split1[2], " ")
            aa <- split2[[1]][1]
            protien <- split2[[1]][2]
            remainder <- strsplit(split1[2], ",")[[1]][2]
            temp <- cases[i, c(1:3, 23, 33,35:38)]
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
## last modified 9/8/16
by.gene <- read.delim("C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/by_gene.txt", stringsAsFactors = F)
by.gene.copy <- by.gene
for (i in 1:20){
    temp <- strsplit(by.gene$gene[i], " ")[[1]]
    if (length(temp) > 1){
        by.gene$gene[i] <- temp[2]
    }else{
        by.gene$gene[i] <- temp[1]
    }
    
}

temp <- table(by.gene$gene)

## Analysis

# Check to see if disrupted or chr loss correlates with subtype
onc.data <- onc.data[onc.data$Oncopanel. == 1, ]

adenomas <- onc.data[onc.data$Pathology == "Pituitary Adenoma", ]

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
by.gene.patient <- PerSampleMaf(by.gene, "gene", "BWH_MRN")
func.adenomas <- by.gene.patient[by.gene.patient$pathology.clinical == "Functional", ]
nonfunc.adenomas <- by.gene.patient[by.gene.patient$pathology.clinical == "Nonfunctional", ]

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
by.gene.adenomas <- FilterMaf(by.gene.patient, c("Functional", "Null"), "pathology.anatomic")
PlotMaf(by.gene.adenomas, "gene", 25, title = "Genes Mutated in at least 4 adenoma patients")
dev.off()


nrow(by.gene.patient) / 106

length(unique(by.gene.patient$gene))

## Hotspot mutations
hot <- ReccurentMaf(by.gene.patient, "amino.acid")
hot.table <- table(hot$gene)
hot.list <- names(hot.table[hot.table > 1])

PlotMaf(by.gene.adenomas[by.gene.adenomas$gene %in% hot.list, ], "gene", "Genes with hotspot mutations")
## GSEA


## plot mutations in functional vs nonfunctional
## Comparison of functional vs nonfunctional tumors
pdf("mutations functional vs nonfunctional.pdf", width = 15, height = 7)
combined.table <- EqualizeTable(table(func.adenomas$gene),table(nonfunc.adenomas$gene))
barplot(combined.table, main = "Comparison of functional and nonfunctional adenomas", 
        legend.text = c("Functional adenomas", "Nonfunctional Adenomas"), las = 2, beside = TRUE, args.legend = list(x = "topleft"))
dev.off()
adenoma.table <- table(by.gene.adenomas$gene)
gh.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "HGH", ]$gene)
gh.df <- data.frame(gh.table, rep("GH", length(gh.table)))
acth.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "ACTH", ]$gene)
acth.df <- data.frame(acth.table, rep("ACTH", length(acth.table)))
prl.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype == "Prolactin", ]$gene)
prl.df <- data.frame(prl.table, rep("PRL", length(prl.table)))
null.table <- table(by.gene.adenomas[by.gene.adenomas$Pathology.Subtype %in% c("Null", "null"), ]$gene)
null.df <- data.frame(null.table, rep("Null", length(null.table)))
colnames(gh.df) <- c("Gene", "Freq", "Subtype")
colnames(acth.df) <- c("Gene", "Freq", "Subtype")
colnames(prl.df) <- c("Gene", "Freq", "Subtype")
colnames(null.df) <- c("Gene", "Freq", "Subtype")
subtype.df <- rbind(gh.df, acth.df, prl.df, null.df)
## creates data frame with per subtype rate
## 

inclusion.list <- dimnames(table(by.gene.adenomas$gene)[table(by.gene.adenomas$gene) > 3])[[1]]
subtype.df <- subtype.df[subtype.df$Gene %in% inclusion.list, ]
ggplot(subtype.df, aes(x = Gene, y = Freq, fill = Subtype)) + 
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


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

