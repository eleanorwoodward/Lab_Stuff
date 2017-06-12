## Decision Tree / Random forest classification for data interprtation
## http://www.saedsayad.com/decision_tree.htm


## generate fake data to make sure that initial split doesn't corrupt data
tree <- df.wide.glioma


random <- sapply(1:20, function(x){rbinom(300,1, .5)}, USE.NAMES = FALSE)
non.random <- sapply(c(.1, .2, .8, .9), function(x){rbinom(300, 1, x)}, USE.NAMES = FALSE)
colnames(non.random) <- rep("non random", ncol(non.random))
correlated <- rbinom(300, 1, .3)
correlates <- cbind(correlated, correlated, correlated, correlated)
colnames(correlates) <- rep("correlates", ncol(correlates))

for (i in 1:ncol(correlates)){
    change.idx <- sample(c(TRUE, FALSE), 300, replace = TRUE, prob = c(0.4, 0.6))
    replacements <- rbinom(sum(change.idx), 1, .3)
    correlates[change.idx, i] <- replacements
}

total <- cbind(random, non.random, correlated, correlates)
tree <- total

## calculates entropy for a given split

entropy.plot1 <- CalculateEntropy(tree[, -1], correction.method = "info", detailed = c("IDH1", "10q", "BRAF", "7q"))

plots.branch1 <- BranchPlotter(entropy.plot1, 4)

pdf("Decision Tree Branch 1.pdf", height = 7, width = 12)
grid.arrange(plots.branch1[[1]], arrangeGrob(plots.branch1[[2]], plots.branch1[[3]], plots.branch1[[4]], plots.branch1[[5]], ncol = 2), widths = c(1, 1.5))
dev.off()

entropy.IDH1_pos <- CalculateEntropy(tree[tree$IDH1 == 1, -1], correction.method = "info", detailed = c("TP53", "19p", "19q", "1p"))


pdf("Decision Tree Branch 2b.pdf", height = 7, width = 12)
grid.arrange(x[[1]], arrangeGrob(x[[2]], x[[3]], x[[4]], x[[5]], ncol = 2), ncol = 2, widths = c(1, 2))
dev.off()



entropy.IDH1_neg <- CalculateEntropy(tree[tree$IDH1 == 0, -1], "info")
sort(entropy.IDH1_neg)
