# first - lets get random gene names from biomaRt
# with 40 genes of interest in 100 subjects
# https://benchtobioinformatics.wordpress.com/2015/05/25/how-to-make-a-co-mutation-plot/
df <- expand.grid(gene=1:15, subj=1:100)
df$Mutation <- sample(c("missense", "nonsense", rep(NA, 10), "nonstop"), replace = T, 1500)
df$gene <- rep(1:15, 100)


df_sub <- subset(df, !is.na(df$Mutation))
df_count <- count(df_sub, "gene")

ord <- df_count[order(df_count$freq), ]$gene

# re-order mutation data.farme by order of frequent mutations
df$gene <- factor(df$gene, levels = df$gene[ord])
df <- df[order(df$gene), ]

table(df$gene)


# now for a Comut plot with ggplot2
mut <- ggplot(df, aes(x=subj, y=gene, height=0.8, width=0.8))
(mut <- mut + geom_tile(aes(fill=Mutation)) +
    scale_fill_brewer(palette = "Set1", na.value="Grey90") +
    xlab("Subject") +
    ggtitle("Example Comut plot") +
    theme(
        legend.key = element_rect(fill='NA'),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size=8, face="bold"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="Black"),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.background=element_blank()
    ))

ggsave(mut,file="2015-5-24-ExampleComutplot.png",width=10,height=8)

values <- rnorm(200, 5, 10)
values.matrix <- matrix(values, 20, 10)
values.d <- dist(values.matrix)
values.clust <- hclust(values.d)
plot(values.clust, main = "copy number clustering-signed", hang = -1, labels = FALSE)
