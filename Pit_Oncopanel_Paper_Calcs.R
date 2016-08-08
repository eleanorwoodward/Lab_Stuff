## Cals for Oncopanel Analysis
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

## Gets CSV data, formats and cleans it
onc.data <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/R_input.txt', stringsAsFactors = FALSE, comment.char = "", 
                   header = T, strip.white = T)
onc.data <- onc.data[1:322, ]

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

by.gene <- read.delim("C:/Users/Noah/Syncplicity Folders/Pituitary Oncopanel Project/Data/by_gene.txt", stringsAsFactors = F)

## Analysis

# Check to see if disrupted or chr loss correlates with subtype
onc.data <- onc.data[onc.data$Oncopanel. == 1, ]

adenomas <- onc.data[onc.data$Pathology == "Pituitary Adenoma", ]

# correlation of 1, 11, and disruption with subtype
yvals <- c("Chr.1.Loss", "Chr.11.loss", "disruption")
mtrx1 <- matrix(NA, 1, 3)
for (i in 1:3){
    mini <- CleanYourRoom(adenomas[, c("pathology.simple", yvals[i])])
    tbl <- table(mini[, 1], mini[, 2])
    math.party <- fisher.test(tbl)
    mtrx1[1, i] <- signif(math.party$p.value, 2)
}
mtrx1


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

## cohort statistics
func.adenomas <- by.gene[by.gene$pathology.simple == "Functional", ]

nonfunc.adenomas <- by.gene[by.gene$pathology.simple == "Null", ]

mutations.func <- table(func.adenomas$BWH_MRN)
mutations.nonfunc <- table(nonfunc.adenomas$BWH_MRN)

t.test(mutations.nonfunc, mutations.func)

## Gene level analysis

PlotMaf(by.gene, "gene", 30, title = "Genes Mutated in at least 5 patients")

by.gene.adenomas <- FilterMaf(by.gene, c("Functional", "Null"), "pathology.simple")

nrow(by.gene) / 106

length(unique(by.gene$gene))

## Hotspot mutations
hot <- ReccurentMaf(by.gene, "amino.acid")




