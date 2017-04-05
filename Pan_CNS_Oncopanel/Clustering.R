## Clustering analysis

## generate distances matrix using dist(). Distances are computed between rows: need to transpose if samples in columns


matrix1 <- matrix(rnorm(352, 2, 5), nrow = 8)
matrix2 <- matrix(rnorm(308, -1, 5), nrow = 7)
matrix3 <- matrix(rnorm(220, 1.5, 5), nrow = 5)

matrix <- rbind(matrix1, matrix2, matrix3)
chromosomes <- paste("chr", sort(rep(1:22, 2)), c("p", "q"), sep = "")
rownames(matrix) <- c(paste("sample", 1:20, sep =  ""))
colnames(matrix) <- chromosomes

d <- dist(matrix)
clust <- hclust(d)
plot(clust)
matrix.df <- data.frame(matrix)
matrix.df$samples <-  rownames(matrix)

samples <- rownames(matrix)
long.df <- melt(matrix.df, id = "samples")

long.df$samples <- factor(long.df$samples, levels = samples[clust$order])
long.df$value[abs(long.df$value) < 3] <- NA
long.df$value[long.df$value > 3] <- 1
long.df$value[long.df$value < -3] <- -1

mut <- ggplot(long.df, aes(x=samples, y=variable, height=0.9, width=0.9)) + 
    geom_tile(aes(fill=value)) +
    scale_fill_gradient(low = "blue", high = "red", na.value = "white", limits = c(-1, 1)) +
    xlab("Subject") +
    ggtitle("Trial Heatmap") +
    theme(
        legend.key = element_rect(fill='NA'),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size=8, face="bold"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.text.y=element_text(colour="Black"),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.background=element_blank()
    )
print(mut) 


## create legend
legend.df <- data.frame(rep(samples, 2), c(rep("tumor", length(samples)), rep("sex", length(samples))), 
                        c(rep(c("glioma", "meingioma"), length(samples)/2), rep(c("male", "female"), length(samples)/2)))
colnames(legend.df) <- c("samples", "phenotype", "values")
l <- ggplot(legend.df, aes(x=samples, y = phenotype)) + coord_fixed(ratio = .2) + geom_tile(aes(fill = values))







fd=data.frame(x = rep(c("x","y","z"),3), 
              y=c("a","b","c","b","c","a","c","a","b"),
              z=c(0,1,0,1,1,1,0,0,1))

# plot
p <- ggplot(fd, aes(x, y)) + geom_tile(aes(fill = z)) 
+ scale_fill_gradient(low = "white",high = "steelblue", limits=c(0,1)) 
+ theme_grey() 
+ labs(x = "", y= "") 
+ coord_fixed(ratio= 1)
+ scale_x_discrete(expand = c(0, 0)) 
+ scale_y_discrete(expand = c(0, 0)) 
+ theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size=12, angle=90, hjust=0, colour="black"))

require(gridExtra)
download.packages("gridExtra")
multiplot(mut, l, cols = 1)
