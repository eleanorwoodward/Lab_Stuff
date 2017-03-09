## calculations for revision
validation.snps[validation.snps$Tumor_Sample_Barcode == "MG-40-tumor" & validation.snps$Hugo_Symbol == "NF2", ][, c("Hugo_Symbol", "Variant_Classification", 
                                                                                                                     "Protein_Change", "Codon_Change", "Start_position")]

validation.indels[validation.indels$Tumor_Sample_Barcode == "MG-193-tumor" &ver validation.indels$Hugo_Symbol == "NF2", ][, c("Hugo_Symbol", "Variant_Classification", 
                                                                                                                        "Protein_Change", "Codon_Change")]

discovery.coding.snps[discovery.coding.snps$Tumor_Sample_Barcode == "MEN0095-P" & discovery.coding.snps$Hugo_Symbol == "NF2",
                    ][, c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "Codon_Change")]

discovery.coding.indels[discovery.coding.indels$Tumor_Sample_Barcode == "MEN0025-P" & discovery.coding.indels$Hugo_Symbol == "NF2",
                      ][, c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "Codon_Change")]

ph.indel[ph.indel$Tumor_Sample_Barcode == "MEN_PH_LG_43-pair", ][, c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "Codon_Change")]
                                                                                      
ph.indel[ph.indel$Tumor_Sample_Barcode == "MEN0019-pair" & ph.indel$Hugo_Symbol == "NF2", ][, c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "Codon_Change")]



## expression analysis
## taken from tutorial at http://homer.salk.edu/homer/basicTutorial/affymetrix.html
file.folder <- "C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Expression/files/"
source("http://bioconductor.org/biocLite.R")

v?# Each of these commands tells Bioconductor to download and install each package
biocLite("affy")
biocLite("oligo")
nbiocLite("limma")

biocLite("hugene10stv1cdf")


# load the oligo library
library(oligo)

# Read in the CEL files in the directory
celFiles <- list.files(file.folder)
setwd(file.folder)
biocLite("hugene10st.db")
affyRaw <- read.celfiles(celFiles)

# You might need to install and load a package for the specific array you are using (this example is mouse gene 2.0 ST)
# It may try to load it automatically, but may fail.  Install & load the library manually if this happens.
library(pd.mogene.2.0.st) 
eset <- rma(affyRaw)

# Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Expression/files/normalized.txt")

expression <- read.delim("C:/Users/noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Expression/files/1445672.expression.gct", stringsAsFactors = F, header = F)
expression.cleaned <- expression
colnames(expression.cleaned) <- expression.cleaned[3, ]
expression.cleaned <- expression.cleaned[-(1:3), ]
table(expression.cleaned$Description == "NA // NA // NA")
expression.gene <- expression.cleaned[expression.cleaned$Description != "NA // NA // NA", ]
expression.gene$gene1 <- "NA"
expression.gene$gene2 <- "NA"
expression.gene$gene3 <- "NA"
expression.gene.named <- expression.gene

for (i in 1:nrow(expression.gene)){
    expression.gene.named[i, 17:19] <- strsplit(expression.gene$Description[i], " // ")[[1]]
    
}

expression.gene.named$mean <- rowSums(data.matrix(expression.gene.named[, 3:16]))
fivenum(expression.gene.named$mean)
boxplot.stats(expression.gene.named$mean)
quantile(expression.gene.named$mean, probs = seq(0,1, .1))
hist(expression.gene.named$mean)
expression.high <- expression.gene.named[expression.gene.named$mean > 73, ]
hist(expression.gene.named$mean)

481 / (1361 + 481)





## check and see which genes show up in oncopanel
samples <- read.delim("C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Oncopanel/Sample_names.txt", stringsAsFactors = F, strip.white = T)
idx <- unique(samples$MG.149) %in% ccgd.snindels$Tumor_Sample_Barcode
unique(samples$MG.149.tumor)[idx]
unique(samples$MG.149.tumor)[!idx]

mutations <- ccgd.snindels[ccgd.snindels$Tumor_Sample_Barcode %in% unique(samples$MG.149.tumor), ]

mutations <- mutations[mutations$Hugo_Symbol %in% genes$Genes, ]
genes <- read.delim("C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Oncopanel/Gene_list.txt", stringsAsFactors = F, strip.white = T)

View(validation.coding.snps[validation.coding.snps$Tumor_Sample_Barcode == "MG-229-tumor", ])
table(validation.coding.snps$Tumor_Sample_Barcode == "MG-59-tumor")

write.csv(mutations, "C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Oncopanel/potential_genes.csv", row.names = F)


genes.val <- read.delim("C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Oncopanel/validation_gene_list.txt", stringsAsFactors = F, strip.white = T)
genes.val <- genes.val[-157, ]

samples.val <- disc.snindels.all[disc.snindels.all$Tumor_Sample_Barcode %in% c("MEN0100-P", "MEN0105-P", "MEN0092-P", "MEN0103-P", "MEN0108-P", 
                                                                               "MEN0119-P1", "MEN0102-P", "MEN0099-P", "MEN0025G-P", "MEN0104-P", "MEN0095-P"), ]
samples.val <- samples.val[samples.val$Hugo_Symbol %in% genes.val, ]

write.csv(samples.val, "C:/Users/Noah/Dropbox/Work/Meningioma/Genomic Medicine Revision/Oncopanel/targeted_validation_overlap.csv", row.names = F)

View(ccgd.snindels[ccgd.snindels$Tumor_Sample_Barcode == "MG-69-tumor", ])
