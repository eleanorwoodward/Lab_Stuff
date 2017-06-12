#SOM analysis
## example taken from https://gist.github.com/dgrapov/f67d0696c4fb02731f55da3e1b9e8c4d

library(kohonen)
data(wines)
set.seed(7)

#create SOM grid
sommap <- som(scale(wines), grid = somgrid(3, 3, "hexagonal"))

## use hierarchical clustering to cluster the codebook vectors
groups<-3
som.hc <- cutree(hclust(dist(sommap$codes[[1]])), groups)

#plot
plot(sommap, type="codes", bgcol=rainbow(groups)[som.hc])

#cluster boundaries
add.cluster.boundaries(sommap, som.hc)


## walkthrough from file:///C:/Users/Noah/Downloads/v21i05.pdf

wines.sc <- scale(wines)
set.seed(7)

wine.som <- som(X = wines.sc, grid = somgrid(3, 2, "hexagonal"))
plot(wine.som, main = "Wine data")
wine.som$codes
wine.som$unit.classif



test.2 <- nmf(testing, 2)

testing.r.3 <- nmf(input.group, 2:7, nrun = 20, seed = 123456)
plot(hello)
plot(testing.r.5)

saveRDS(testing.r.3, "../Upload to xchip/testing.object.RDS")
hello <- readRDS("../Upload to xchip/NMF_OUPUT.RDS")

coefmap(minfit(testing.r.3), annCol = colnames(input.group), annColors = list(c("red", "yellow", "blue", "green4", "purple", "orange")), Colv = "basis")


consensusmap(testing.r.5)
