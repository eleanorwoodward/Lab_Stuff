## Pituitary Tumors Analysis

pancan_names <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "MEN", "OV", "PIT", "READ", "UCEC")
data = matrix(0, 3, length(pancan_names))
colnames(data) <- pancan_names
broad.list <- c()
focal.list <- c()
pit.focal.list <- c()
pit.broad.list <- c()
pit.disrupted <- c(1, 4, 9, 16, 22, 23, 26, 28, 32, 34, 37, 42)
pit.disrupted.list <- c("PIT001-Tumor", "PIT007-Tumor", "PIT013-Tumor", "PIT1002-Tumor","PIT1009-Tumor","PIT1010-Tumor","PIT1013-Tumor","PIT1015-Tumor",
                        "PIT1022-Tumor","PIT202-Tumor", "PIT311-Tumor", "PIT504-Tumor")

for (i in 1:length(pancan_names)){
  
 
broad.names <- list.files('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/broad/')
focal.names <- list.files('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/')
  
broad <- read.delim(paste('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/broad/', 
                          broad.names[i], sep = ""), stringsAsFactors = FALSE)
focal <- read.delim(paste('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/genomedisruption/output_files_ng_edited/focal/', 
                          focal.names[i], sep = ""), stringsAsFactors = FALSE)

if (pancan_names[i] == "PIT"){
  pit.focal.list <- focal$focal_disruption_from_median
  pit.broad.list <- broad$broad_disruption_from_median
  
}else{
avg_focal = mean(focal$focal_disruption_from_median)
avg_broad = mean(broad$broad_disruption_from_median)
ratio <- avg_focal / avg_broad
data[1,i] <- avg_focal
data[2, i] <- avg_broad
data[3, i] <- ratio
broad.list <- c(broad.list, broad$broad_disruption_from_median)
focal.list <- c(focal.list, focal$focal_disruption_from_median)
}
}
rownames(data) <- c("average focal", "average broad", "ratio")
data
t.test(focal.list, pit.focal.list[pit.disrupted])
mean(focal.list)
mean(pit.focal.list)
barplot(data[3, ], main = "Comparison of Ratio of Focal to Broad in various tumor types")


## Allelic fraction information per tumor
setwd("C:/Users/Noah/OneDrive/Work/R/output/")
x <- read.delim('C:/Users/Noah/Syncplicity Folders/Pituitary (Linda Bi)/Mutation-Indels/Pit20140506.final_analysis_set.maf.txt', 
                stringsAsFactors = FALSE)
sample.list <- unique(x$Tumor_Sample_Barcode)
for (i in 1:length(sample.list)){
  fig.name <- sample.list[2]
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


