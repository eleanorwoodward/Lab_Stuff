## analysis of mutation data from Oncopanel calls

source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

setwd("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/OncDRS data/")


## read in raw data
all.mutations <- read.csv("REQ_ID08_65337_ONCOPANEL_MUTATION_RESULTS.csv", stringsAsFactors = F)
all.mutations.tier1.2 <- all.mutations[all.mutations$TIER_ID < 3, ]

all.cancer.diag <- read.csv("REQ_ID08_65337_CANCER_DIAGNOSIS_CAREG.csv", stringsAsFactors = F)



## define groupings for analysis
gliomas <- c("Anaplastic Astrocytoma", "Astrocytoma", "Ganglioglioma", "Anaplastic Oligoastrocytoma", "Diffuse Glioma", "Glioblastoma", "Low-Grade Glioma, NOS", 
             "Oligoastrocytoma", "Paraganglioma", "Anaplastic Ganglioglioma", "Anaplastic Oligodendroglioma", "Glioblastoma Multiforme", "Small Cell Glioblastoma")

meningiomas <- c("Anaplastic Meningioma", "Atypical Meningioma", "Meningioma", "Rhabdoid Meningioma")






## gets all histologies in dataset
write.table(table(all.mutations[!(duplicated(all.mutations$PATIENT_ID)), ]$PRIMARY_CANCER_DIAGNOSIS), "../Analysis/all_tumor_types.txt", sep = "\t", quote = F, 
            row.names = F)

## all histologies other than major ones
other.mutations <- all.mutations[!(all.mutations$PRIMARY_CANCER_DIAGNOSIS %in% c(gliomas, meningiomas)), ]
table(other.mutations[!(duplicated(other.mutations$PATIENT_ID)), ]$PRIMARY_CANCER_DIAGNOSIS)

table (all.cancer.diag$HISTOLOGY_DESCR)

length(unique(all.cancer.diag$PATIENT_ID))

length(unique(all.mutations[all.mutations$BIOPSY_SITE_TYPE != "Metastatic Recurrence", ]$PATIENT_ID))



## exploratory analysis plots

## most frequently mutated genes

## all tier 1 & 2 mutations in all samples
pdf("../Analysis/all_tier1and2_mutations.pdf", width = 14)
PlotMaf(all.mutations.tier1.2, "BEST_EFF_GENE", title = "All Tier 1 and 2 mutations in all cancers")
dev.off()

## all tier 1 & 2 mutations in all gliomas
pdf("../Analysis/all_tier1and2_mutations_gliomas.pdf", width = 14)
PlotMaf(all.mutations.tier1.2[all.mutations.tier1.2$PRIMARY_CANCER_DIAGNOSIS %in% gliomas, ], "BEST_EFF_GENE", title = "All Tier 1 and 2 mutations in all gliomas")
dev.off()

## all tier 1 & 2 mutations in all meningiomas
pdf("../Analysis/all_tier1and2_mutations_meningiomas.pdf", width = 14)
PlotMaf(all.mutations.tier1.2[all.mutations.tier1.2$PRIMARY_CANCER_DIAGNOSIS %in% meningiomas, ], "BEST_EFF_GENE", title = "All Tier 1 and 2 mutations in all meningiomas")
dev.off()


## all tier 1 & 2 mutations in non-meningioma, non-gliomass
pdf("../Analysis/all_tier1and2_mutations_meningiomas.pdf", width = 14)
PlotMaf(all.mutations.tier1.2[!(all.mutations.tier1.2$PRIMARY_CANCER_DIAGNOSIS %in% c(meningiomas, gliomas)), ], "BEST_EFF_GENE", 
        title = "All tier 1 and 2 mutations in all non-meningiomas, non gliomas")
dev.off()



## recurrent hotspots across the cohort
## Hotspot finder
snindels11 <- ReccurentMaf(snindels, "Hugo_Symbol", 4)
gene.list <- sort(unique(snindels11$Hugo_Symbol))
mtrx <- matrix(0, 40, length(gene.list))
for(i in 1:length(gene.list)){
    temp <- FilterMaf(snindels11, gene.list[i], "Hugo_Symbol")
    x <- table(temp$Start_position)
    x <- sort(unname(x), decreasing = TRUE)
    if (length(x) > 40){
        x <- x[1:40]
    }
    
    if (length(x) < 40){
        temp <- 40 - length(x)
        x <- c(x, rep(0, temp))
    }
    mtrx[,i] <- x
}

