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
