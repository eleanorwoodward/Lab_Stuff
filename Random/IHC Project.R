## analysis of IHC data from CNS tumors
setwd("C:/Users/Noah/Syncplicity Folders/Immune (Linda Bi)/Path Core Staining/aperio scope data/367-pituitary tumor/")
values <- read.delim("C:/Users/Noah/Syncplicity Folders/Immune (Linda Bi)/Path Core Staining/aperio scope data/367-pituitary tumor/Input_to_R_data.txt", stringsAsFactors = FALSE)
values <- values[ -nrow(values), ]

## calculate the variance among each group of three samples


col1.idx.1 <- rep(c(T, F, F), nrow(values) / 3)
col1.idx.2 <- rep(c(F, T, F), nrow(values) / 3)
col1.idx.3 <- rep(c(F, F, T), nrow(values) / 3)

current.marker <- values$VISTA

current.marker.df <- as.data.frame(cbind(current.marker[col1.idx.1], current.marker[col1.idx.2], current.marker[col1.idx.3]), stringsAsFactors = FALSE)
current.marker.df <- current.marker.df[-c(32, 64, 96)]


current.marker.df$variance <- as.numeric(apply(current.marker.df, 1, var))
current.marker.df$mean <- as.numeric(apply(current.marker.df[, 1:3], 1, mean))

raw.data <- as.numeric(unlist(current.marker.df[, 1:3]))
data2 <- as.data.frame(t(sapply(1:96, function(x){sample(raw.data, 3)})))
data2$variance <- as.numeric(apply(data2, 1, var))
data2$mean <-  as.numeric(apply(data2[, 1:3], 1, mean))

current.marker.df$cohort <- "original"
data2$cohort <- "control"
combined <- rbind(current.marker.df, data2)

#ggplot(data = data, aes(x = variance)) + geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 30)
pdf("LAG3 Comparison.pdf")
ggplot(data = combined, aes(x = cohort, y = variance)) + geom_dotplot(binaxis = "y", stackdir = "center")
dev.off()

pdf("VISTA Comparison Mean.pdf")
ggplot(data = combined, aes(x = cohort, y = mean)) + geom_dotplot(binaxis = "y", stackdir = "center")
dev.off()

