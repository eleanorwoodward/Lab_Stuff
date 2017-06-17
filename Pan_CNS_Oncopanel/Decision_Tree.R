## Decision Tree / Random forest classification for data interprtation
## http://www.saedsayad.com/decision_tree.htm
source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/MafFunctions.R")

## generate fake data for testing
random <- sapply(1:20, function(x){rbinom(1200,1, .5)}, USE.NAMES = FALSE)
colnames(random) <- paste("random", 1:ncol(random))
non.random <- sapply(c(.1, .2, .8, .9), function(x){rbinom(1200, 1, x)}, USE.NAMES = FALSE)
colnames(non.random) <- paste("non random", 1:ncol(non.random))
correlated <- rbinom(1200, 1, .4)
correlates <- cbind(correlated, correlated, correlated, correlated)
colnames(correlates) <- paste("correlates", 1:ncol(correlates))

for (i in 1:ncol(correlates)){
    change.idx <- sample(c(TRUE, FALSE), 1200, replace = TRUE, prob = c(0.4, 0.6))
    replacements <- rbinom(sum(change.idx), 1, .4)
    correlates[change.idx, i] <- replacements
}

total <- cbind(random, non.random, correlated, correlates)
tree <- apply(total, 2, as.numeric)


entropy.fake_data <- CalculateEntropy(tree, correction.method = "info", detailed = c("correlated", "correlates 1", "correlates 2", "non random 4"))
plots.fake_data <- BranchPlotter(entropy.fake_data, 4, TRUE)
grid.arrange(arrangeGrob(plots.fake_data[[6]], plots.fake_data[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.fake_data[[2]], plots.fake_data[[3]], plots.fake_data[[4]], plots.fake_data[[5]], ncol = 2), widths = c(1, 1.5))

PermuteDecisionTree(tree, 100)


tree <- apply(df.wide.glioma[, -1], 2, as.numeric)
tree <- as.data.frame(tree, stringsAsFactors = FALSE)

## work iteratively down branches of each tree
entropy.branch1 <- CalculateEntropy(tree[, ], correction.method = "info", detailed = c("IDH1", "10q", "BRAF", "7q"))
plots.branch1 <- BranchPlotter(entropy.branch1, 4)
pdf("Decision Tree Branch 1.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch1[[6]], plots.branch1[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch1[[2]], plots.branch1[[3]], plots.branch1[[4]], plots.branch1[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch1.permuted <- PermuteDecisionTree(tree, 200)


entropy.branch2 <- CalculateEntropy(tree[tree$IDH1 == 1, ], correction.method = "info", detailed = c("TP53", "19p", "19q", "1p"))
plots.branch2 <- BranchPlotter(entropy.branch2, 4)
pdf("Decision Tree Branch 2.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch2[[6]], plots.branch2[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch2[[2]], plots.branch2[[3]], plots.branch2[[4]], plots.branch2[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch2.permuted <- PermuteDecisionTree(tree[tree$IDH1 == 1, colnames(tree) != "IDH1"],200)
sapply(sort(entropy.branch2[[1]], decreasing = TRUE)[1:4], function(x){sum(branch2.permuted > x) / length(branch2.permuted)})

entropy.branch3 <- CalculateEntropy(tree[tree$IDH1 == 0, ], correction.method = "info", detailed = c("BRAF", "10q", "7q", "7p"))
plots.branch3 <- BranchPlotter(entropy.branch3, hit.number = 4, incidence = TRUE)
pdf("Decision Tree Branch 3.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch3[[6]], plots.branch3[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch3[[2]], plots.branch3[[3]], plots.branch3[[4]], plots.branch3[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch3.permuted <- PermuteDecisionTree(tree[tree$IDH1 == 0, colnames(tree) != "IDH1"],200)
sapply(sort(entropy.branch3[[1]], decreasing = TRUE)[1:6], function(x){sum(branch3.permuted > x) / length(branch3.permuted)})


entropy.branch4 <- CalculateEntropy(tree[tree$IDH1 == 1 & tree$TP53 == 1, ], correction.method = "info", detailed = c("19p", "CDK4", "7q", "ATRX"))
plots.branch4 <- BranchPlotter(entropy.branch4, hit.number = 4, incidence = TRUE)
pdf("Decision Tree Branch 4.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch4[[6]], plots.branch4[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch4[[2]], plots.branch4[[3]], plots.branch4[[4]], plots.branch4[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch4.permuted <- PermuteDecisionTree(tree[tree$IDH1 == 1 & tree$TP53 == 1, !(colnames(tree) %in% c("IDH1", "TP53"))],200)
sapply(sort(entropy.branch4[[1]], decreasing = TRUE)[1:6], function(x){sum(branch4.permuted > x) / length(branch4.permuted)})


entropy.branch5 <- CalculateEntropy(tree[tree$IDH1 == 1 & tree$TP53 == 0, ], correction.method = "info", detailed = c("19p", "1p", "19q", "PIK3R1"))
plots.branch5 <- BranchPlotter(entropy.branch5, hit.number = 4, incidence = TRUE)
pdf("Decision Tree Branch 5.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch5[[6]], plots.branch5[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch5[[2]], plots.branch5[[3]], plots.branch5[[4]], plots.branch5[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch5.permuted <- PermuteDecisionTree(tree[tree$IDH1 == 1 & tree$TP53 == 0, !(colnames(tree) %in% c("IDH1", "TP53", "CDK4"))],200)
sapply(sort(entropy.branch5[[1]], decreasing = TRUE)[1:6], function(x){sum(branch5.permuted > x) / length(branch5.permuted)})


entropy.branch6 <- CalculateEntropy(tree[tree$IDH1 == 0 & tree$BRAF == 1, ], correction.method = "info", detailed = c("1p", "ARID1B", "7q", "CDKN2A"))
plots.branch6 <- BranchPlotter(entropy.branch6, hit.number = 4, incidence = TRUE)
pdf("Decision Tree Branch 6.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch6[[6]], plots.branch6[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch6[[2]], plots.branch6[[3]], plots.branch6[[4]], plots.branch6[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch6.permuted <- PermuteDecisionTree(tree[tree$IDH1 == 1 & tree$BRAF == 1, !(colnames(tree) %in% c("IDH1", "BRAF", "CDK4", "PIK3C2B", "19p"))],200)
sapply(sort(entropy.branch6[[1]], decreasing = TRUE)[1:6], function(x){sum(branch6.permuted > x) / length(branch6.permuted)})


entropy.branch7 <- CalculateEntropy(tree[tree$IDH1 == 0 & tree$BRAF == 0, ], correction.method = "info", detailed = c("10q", "FGFR1", "7q", "EGFR"))
plots.branch7 <- BranchPlotter(entropy.branch7, hit.number = 4, incidence = TRUE)
pdf("Decision Tree Branch 7.pdf", height = 7, width = 12)
grid.arrange(arrangeGrob(plots.branch7[[6]], plots.branch7[[1]], nrow = 2, heights = c(1, 3)), 
             arrangeGrob(plots.branch7[[2]], plots.branch7[[3]], plots.branch7[[4]], plots.branch7[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()
branch7.permuted <- PermuteDecisionTree(tree[tree$IDH1 == 1 & tree$BRAF == 0, !(colnames(tree) %in% c("IDH1", "BRAF"))],200)
sapply(sort(entropy.branch7[[1]], decreasing = TRUE)[1:6], function(x){sum(branch7.permuted > x) / length(branch7.permuted)})




## generate decision tree to show order of selected features
pyramid <- data.frame(c(3,2,2,1,1,1,1), c(4,2.5,5.5,2,3,5,6), c("Decision Tree (1)", "IDH1-mt (2)", "IDH1-wt (3)", 
                                                                "TP53-mt (4)", "TP53-wt (5)", "BRAF-mt (6)", "BRAF-wt (7)"))
colnames(pyramid) <- c("height", "width", "names")
colnames(tree)[1] <- "SAMPLE_ACCESSION_NBR"
tree <- merge(tree, master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Specific")])

tree$branch <- NA
tree$branch[tree$IDH1 == 1 & tree$TP53 == 1] <- "Branch4"
tree$branch[tree$IDH1 == 1 & tree$TP53 == 0] <- "Branch5"
tree$branch[tree$IDH1 == 0 & tree$BRAF == 1] <- "Branch6"
tree$branch[tree$IDH1 == 0 & tree$BRAF == 0] <- "Branch7"

temp <- tree[tree$Cancer_Type_Specific %in% names(sort(table(tree$Cancer_Type_Specific), decreasing = TRUE))[1:9], ]
pdf("Branch Subtype Enrichment.pdf")
dist.colors <- distinctColorPalette(length(unique(temp$Cancer_Type_Specific)))
names(dist.colors) <- unique(temp$Cancer_Type_Specific)
ggplot(data = temp[, c("branch", "Cancer_Type_Specific")], aes(x = branch, fill = Cancer_Type_Specific)) + geom_bar() +
    scale_fill_manual(values = dist.colors)
dev.off()

pdf("Tree diagram.pdf", width = 8, height = 5)
ggplot(data = pyramid, aes (x = width, y = height)) + geom_point(color = "grey") + geom_text(aes(label = names)) + scale_x_continuous(limits = c(1.5, 6.5))
dev.off()








