# SCNA analysis

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")

distances <- dist(gistic.calls)
model <- hclust(distances)
plot(model)
sample.info <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]
sample.info <- sample.info[model$order, ]
sample.info <- t(sample.info)
write.csv(sample.info, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Figs/cluster.label.csv")
sample.names <- sample.info[3, ]

write.table(sample.info[3,], "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/figs/clustering/cluster_order.txt", sep = "\t", row.names = F, quote = F)

