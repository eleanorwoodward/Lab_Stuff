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

master.table[master.table$Grade == "I" & master.table$Cohort != "onc", ]$percent.disruption
master.table[master.table$Grade != "I" & master.table$Cohort != "onc", ]$percent.disruption
master.table[master.table$Grade == "III" & master.table$Cohort != "onc", ]$percent.disruption
rep(master.table[master.table$Subtype == "Rhabdoid", ]$percent.disruption, 3)

copy.number <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/GISTIC/total.gistic.10.14/broad_values_by_arm.txt", stringsAsFactors = F)
copy.number <- t(copy.number)
colnames(copy.number) <- copy.number[1,]
copy.number <- copy.number[-1, ]
copy.number <- cbind(copy.number, 0)
copy.number[,40] <- sapply(1:57, function(i) sum(as.numeric(copy.number[i,]) < 0))
copy.number <- cbind(copy.number, 0)
copy.number[,41] <- sapply(1:57, function(i) sum(as.numeric(copy.number[i,]) > 0))

rhaboid.tumors <- c(12,22,41:47, 51:52)
hg.no.rhaboid <- copy.number[c(11, 13:21, 23:40, 48:50), 40:41]
hg.rhaboid <- copy.number[rhaboid.tumors, 40:41]

# 
# lg <- rbind(copy.number[1:10, 40:41], copy.number[53:57, 40:41])
# hg <- copy.number[11:52, 40:41]
# loss.counts <- c(sum(as.numeric(lg[,1])), sum(as.numeric(hg[,1])))
# gain.counts <- c(sum(as.numeric(lg[,2])), sum(as.numeric(hg[,2])))
# prop.test(cbind(loss.counts, gain.counts))
# 
# loss.counts <- c(sum(as.numeric(hg.no.rhaboid[,1])), sum(as.numeric(hg.rhaboid[,1])))
# gain.counts <- c(sum(as.numeric(hg.no.rhaboid[,2])), sum(as.numeric(hg.rhaboid[,2])))
# prop.test(cbind(loss.counts, gain.counts))


## Gene Pattern Formating
gp.table <- master.table[master.table$Analsysis.Set. == 1 | master.table$Cohort == "PH", ]
gp.table <- gp.table[, c("Pair.Name", "Grade", "Subtype", "NF2.snp.indel.rearrangement", "chr22.loss", "Chr1.loss", "X3p.loss", "X6p.loss", 
                             "X6q.loss", "X7p.loss", "X14q.loss", "X17q.gain", "X18p.loss", "X18q.loss", "X20q.gain")]
gp.flip <- t(gp.table)
write.csv(gp.flip, "C:/Users/Noah/OneDrive/Work/Meningioma/GenePattern/master.table.csv", col.names = F)

## not significant
nof2.vals <- master.table[master.table$Analsysis.Set. == 1 & master.table$NF2.snp.indel.rearrangement == 0  & master.table$chr22.loss == 0, ]$percent.disruption
nf2.vals <- master.table[master.table$Analsysis.Set. == 1 & (master.table$NF2.snp.indel.rearrangement == 1  | master.table$chr22.loss == 1), ]$percent.disruption
nf2.vals <- nf2.vals[!is.na(nf2.vals)]
t.test(nf2.vals, nof2.vals)
