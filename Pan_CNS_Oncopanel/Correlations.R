## correlations analysis


## construct master comut
all.samples <- unique(master.sheet[master.sheet$PANEL_VERSION > 0, ]$SAMPLE_ACCESSION_NBR)
all.genes <- unique(all.mutations.tier1.4$BEST_EFF_GENE)
df <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
colnames(df) <- c("samples", "genes")

master.maf <- melt(all.mutations.tier1.4[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")], id = c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE"))
master.maf <- master.maf[-3]
colnames(master.maf) <- c("samples", "genes", "mutations")
df <- merge(df, master.maf, c("samples", "genes"), all.x = TRUE)
df[!is.na(df$mutations), ]$mutations <- 1
df[is.na(df$mutations), ]$mutations <- 0

df.pheno <- melt(master.sheet[, c("SAMPLE_ACCESSION_NBR", "Cancer_Type_Broad")], id = "SAMPLE_ACCESSION_NBR")
colnames(df.pheno) <- c("samples", "genes", "mutations")
df <- rbind(df, df.pheno)

df.wide <- reshape(df, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")

## set broad types
df.wide$simple.pathology <- "other"
df.wide$simple.pathology[df.wide$mutations.Cancer_Type_Broad == "Glioma"] <- "Glioma"
df.wide$simple.pathology[df.wide$mutations.Cancer_Type_Broad == "Meningioma"] <- "Meningioma"
df.wide$simple.pathology[df.wide$mutations.Cancer_Type_Broad == "Metastasis"] <- "Metastasis"
df.wide$simple.pathology[df.wide$mutations.Cancer_Type_Broad == "Pituitary_adenoma"] <- "Pituitary_adenoma"


## for 1065 samples, what percent of samples need to be mutated before we see significant result for overlap?
fisher.test(matrix(c(1,9,9,81), nrow= 2))


## given starting conditions, produces overlap that would be expected by default

row.marginal <- c(980, 85, 1065)
col.marginal <- c(1013, 52, 1065)

matrix.1 <- floor(matrix(c(row.marginal[1] * (col.marginal[1] / col.marginal[3]), row.marginal[2] * (col.marginal[1] / col.marginal[3]), 
              row.marginal[1] * (col.marginal[2] / col.marginal[3]), row.marginal[2] * (col.marginal[2] / col.marginal[3])), nrow = 2))


fisher.test(matrix.1)
fisher.test(matrix(c(5,5,5,980), nrow = 2))

## for co-occurence, essentially no lower bound on what we can detect, because expectation is that there is essentially no overlap. However,
## to detect anti-corrlelations, expectation needs to be that there is some overlap

fisher.test(matrix(c(0,100,100,800), nrow= 2))
fisher.test(matrix(c(0,50,50,900), nrow= 2))

df.wide.glioma <- df.wide[df.wide$simple.pathology == "Glioma", ]
idx <- apply(df.wide.glioma[, c(-1, -521, -522)], 2, function(x){sum(as.numeric(x))})
df.wide.glioma <- df.wide.glioma[, c(TRUE, idx > 50, FALSE, FALSE)]

colnames(df.wide.glioma) <- c(gsub("mutations.", "", colnames(df.wide.glioma)))

## compute all pairwise comparisons for co-occurence
gene1 <- c()
gene2 <- c()
p.vals <- c()
or <- c()


for (i in 2:(ncol(df.wide.glioma) - 2)){
    ## get gene names
    current.name <- colnames(df.wide.glioma)[i]
    current <- df.wide.glioma[, i]
    targets <- colnames(df.wide.glioma)[-c(1:i)]
    gene1 <- c(gene1, rep(current.name, length(targets)))
    gene2 <- c(gene2, targets)
    p.vals <- c(p.vals, as.numeric(apply(df.wide.glioma[, -c(1:i)], 2, function(x){(fisher.test(table(x, current)))$p.value})))
    or <- c(or, as.numeric(apply(df.wide.glioma[, -c(1:i)], 2, function(x){(fisher.test(table(x, current)))$estimate})))
}

## do last one manually
i <- ncol(df.wide.glioma) - 1
current.name <- colnames(df.wide.glioma)[i]
current <- df.wide.glioma[, i]
targets <- colnames(df.wide.glioma)[-c(1:i)]
gene1 <- c(gene1, rep(current.name, length(targets)))
gene2 <- c(gene2, targets)
p.vals <- c(p.vals, fisher.test(table(df.wide.glioma[, i + 1], current))$p.value)
or <- c(or, fisher.test(table(df.wide.glioma[, i + 1], current))$estimate)

correlations <- data.frame(gene1, gene2, p.vals, or)
correlations$q.vals <- p.adjust(correlations$p.vals, "fdr")

View(correlations[correlations$or < 1, ])


## generate function to perform permutations testing
shuffle <- df[df$samples %in% pathologies[["glioma"]], ]
shuffle <- shuffle[shuffle$mutations == 1, ]
shuffle <- shuffle[, -3]
shuffle <- shuffle[!duplicated(shuffle), ]


PermuteData <- function(gene1, gene2, reps, data){
    output <- rep(0, reps)
    for (i in 1:reps){
        data$genes <- sample(data$genes, nrow(data))
        output[i] <- sum(unique(data[data$genes == gene1, ]$samples) %in% data[data$genes == gene2, ]$samples)
    }
    return(output)
}


## insert permutated p values for top hits
correlations <- correlations[order(correlations$q.vals), ]
correlations$permuted.p.value <- NA

for (i in 100:115){
    expected <- PermuteData(correlations$gene1[i], correlations$gene2[i], 10000, shuffle)
    observed <- sum(unique(shuffle[shuffle$genes == correlations$gene1[i], ]$samples) %in% shuffle[shuffle$genes == correlations$gene2[i], ]$samples)
    correlations$permuted.p.value[i] <- sum(expected > observed) / 10000
}


write.csv(correlations, "C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/Analysis/correlations.csv", row.names = F)



