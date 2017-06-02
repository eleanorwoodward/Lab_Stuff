## correlations analysis


## generate list of mutations
master.snvs <- all.mutations.tier1.4[all.mutations.tier1.4$SAMPLE_ACCESSION_NBR %in% pathologies[[4]], ]
tier.3.muts <- all.mutations.tier1.3[all.mutations.tier1.3$SAMPLE_ACCESSION_NBR %in% pathologies[[4]], ]
master.snvs <- master.snvs[master.snvs$BEST_EFF_GENE %in% tier.3.muts$BEST_EFF_GENE, ]
master.snvs <- master.snvs[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]

## focal CNVs
master.cnv.focal <- all.cnv.high[all.cnv.high$SAMPLE_ACCESSION_NBR %in% pathologies[[4]], c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD")]
master.cnv.focal <- all.cnv.high[, c("SAMPLE_ACCESSION_NBR", "GENE", "CNV_TYPE_CD")]
colnames(master.cnv.focal) <- c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")

master.svs <- all.svs.formatted
master.svs <- master.svs[, c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]

master.alterations <- rbind(master.snvs, master.svs, master.cnv.focal)
master.alterations <- master.alterations[master.alterations$BEST_EFF_GENE %in% names(sort(table(master.alterations$BEST_EFF_GENE), decreasing = TRUE))[1:25], ]

## arm level cnvs
master.cnv.broad <- abs(all.cnvs.broad)
master.cnv.broad <- master.cnv.broad[, colnames(master.cnv.broad) %in% pathologies[[4]]]
master.cnv.broad <- master.cnv.broad[rownames(master.cnv.broad) %in% c("1p", "7p", "7q", "10q", "19q", "19p", "22q"), ]


## clinical correlates
master.phenotypes <- master.sheet[master.sheet$Cancer_Type_Broad == "Glioma", ]
master.phenotypes <- master.phenotypes[!(master.phenotypes$exclude %in% c("duplicate", 1)), ]
master.phenotypes <- master.phenotypes[, c("SAMPLE_ACCESSION_NBR", "Primary", "Age_Interval")]


master.alterations <- master.alterations[master.alterations$SAMPLE_ACCESSION_NBR %in% master.phenotypes$SAMPLE_ACCESSION_NBR, ]
master.cnv.broad <- master.cnv.broad[, colnames(master.cnv.broad) %in% master.phenotypes$SAMPLE_ACCESSION_NBR]

## generate data
glioma.alterations <- PlotComut(
    mut.maf1 = master.alterations, 
    broad.cnv.maf = master.cnv.broad,
    samples = master.phenotypes$SAMPLE_ACCESSION_NBR,
    input.samples = "SAMPLE_ACCESSION_NBR", 
    input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification",
    dimensions = c(18,10), 
    col.vector = c("damaging mutation" = "orange", "focal loss" = "blue", "focal gain" = "red", "rearrangement" = "green", "missense" = "skyblue",
                   "arm-level gain" = "red", "arm-level loss" = "blue","Glioblastoma" = "black", "AnaplasticAstro" = "green", "Pediatric" = "purple",
                   "Angiocentric" = "orange", "Astro" = "green", "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                   "0" = "gray90", "1" = "gray70", "2" = "gray50", "3" = "gray25", "4" = "gray1"),
    #fixed.order = clust$labels[clust$order],
    #manual.order = c("Cancer_Type_Specific", "TP53", "EGFR"),
    #file.path = "../Analysis/GBM Comut All Alterations by Cancer Type TP53 EGFR.pdf",
    return.matrix = TRUE,
    #phenotypes = clinical[, -c(4:5)]
    )


helper <- function(x){
    sapply(x, function(y){if (y=="wt"){return(0)}else {return(1)}}, USE.NAMES = FALSE)
}

glioma.mutations.wide <- reshape(glioma.alterations, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")
df.wide.glioma <- glioma.mutations.wide
df.wide.glioma <- df.wide.glioma[!is.na(df.wide.glioma$mutations.7p), ]
df.wide.glioma <- as.data.frame(cbind(as.character(df.wide.glioma[, 1]), apply(df.wide.glioma[, -1], 2, helper)), stringsAsFactors = FALSE)

colnames(df.wide.glioma) <- c(gsub("mutations.", "", colnames(df.wide.glioma)))


## converts back to binary
df.wide.glioma$Primary <- 0
df.wide.glioma$Primary[glioma.mutations.wide[!is.na(glioma.mutations.wide$mutations.7p), ]$mutations.Primary == 1] <- 1
df.wide.glioma$Age_interval <- 0
df.wide.glioma$Age_interval[glioma.mutations.wide[!is.na(glioma.mutations.wide$mutations.7p), ]$mutations.Age_interval == "Old"] <- 1


## automated loop for comparing correlations from two different conditions
gbms <- master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Specific == "Glioblastoma"]
input.1 <- df.wide.glioma
input.1 <- df.wide.glioma[(df.wide.glioma$V1 %in% gbms), ]

input.2 <- input.1[(input.1$V1 %in% gbms), ]
input.2 <- input.1[input.1$TP53 == 0, -2]
input.2 <- input.1[input.1$IDH1 == 1, -4]
input.2 <- input.1[input.1$EGFR == 0, -3]

inputs <- list(input.1, input.2)
path <- "correlations_heatmap_gbm_vs_EGFR_wt.pdf"
for (i in 1:length(inputs)){
    
    ## initialize values
    gene1 <- c()
    gene2 <- c()
    p.vals <- c()
    or <- c()
    
    df.correlations.input <- inputs[[i]]
    ## remove columns with no entries
    df.matrix <- apply(df.correlations.input[, -1], 2, as.numeric)
    df.correlations.input <- df.correlations.input[, c(TRUE, as.numeric(colSums(df.matrix)) != 0)]
    
    for (j in 2:(ncol(df.correlations.input) - 2)){
        ## get gene names
        current.name <- colnames(df.correlations.input)[j]
        current <- df.correlations.input[, j]
        targets <- colnames(df.correlations.input)[-c(1:j)]
        gene1 <- c(gene1, rep(current.name, length(targets)))
        gene2 <- c(gene2, targets)
        
        ## calculate p vals across all columns for each phenotype iteratively
        p.vals <- c(p.vals, as.numeric(apply(df.correlations.input[, -c(1:j)], 2, function(x){(fisher.test(table(x, current)))$p.value})))
        or <- c(or, as.numeric(apply(df.correlations.input[, -c(1:j)], 2, function(x){(fisher.test(table(x, current)))$estimate})))
    }
    
    ## do last one manually
    j <- ncol(df.correlations.input) - 1
    current.name <- colnames(df.correlations.input)[j]
    current <- df.correlations.input[, j]
    targets <- colnames(df.correlations.input)[-c(1:j)]
    gene1 <- c(gene1, rep(current.name, length(targets)))
    gene2 <- c(gene2, targets)
    p.vals <- c(p.vals, fisher.test(table(df.correlations.input[, j + 1], current))$p.value)
    
    ## generate data frame
    or <- c(or, fisher.test(table(df.correlations.input[, j + 1], current))$estimate)
    correlations <- data.frame(gene1, gene2, p.vals, or)
    correlations$q.vals <- p.adjust(correlations$p.vals, "fdr")
    
    
    ## heatmap for correlations
    temp <- correlations
    temp$q.vals.plotting <- -log10(temp$q.vals)
    temp$q.vals.plotting[temp$q.vals.plotting > 30] <- 30
    temp$q.vals.plotting[temp$q.vals.plotting < 2] <- NA
    temp$q.vals.plotting[temp$or < 1] <- -temp$q.vals.plotting[temp$or < 1]
    temp$gene1 <- factor(temp$gene1, levels = colnames(df.correlations.input)[-1])
    temp$gene2 <- factor(temp$gene2, levels = colnames(df.correlations.input)[-1])
    
    
    if (i == 1){
        ## for baseline comparison, save the following for default values
        baseline <- temp
        baseline <- baseline[!is.na(baseline$q.vals.plotting), ]
        gene1.values <- as.character(unique(baseline$gene1))
        gene2.values <- as.character(unique(baseline$gene2))
    }else{
        ## for modified comparison, remove genes that aren't in baseline unless q.value is significant (new correlation)
        modified <- temp
        modified <- modified[(modified$gene1 %in% gene1.values) | !is.na(modified$q.vals.plotting), ]
        modified <- modified[(modified$gene2 %in% gene2.values) | !is.na(modified$q.vals.plotting), ]
        
        ## of genes remaining, figure out which gene-gene pairs were significant in baseline
        dups <- duplicated(rbind(baseline[, 1:2], modified[, 1:2]))
        dups <- dups[-(1:nrow(baseline))]
        modified <- modified[dups | !is.na(modified$q.vals.plotting), ]
        modified$q.vals.plotting[is.na(modified$q.vals.plotting)] <- 0
        
        ## determine which significant genes in baseline are no longer significant
        missing.comp <- modified[modified$q.vals.plotting == 0, 1:2]
        missing.idx <- duplicated(rbind(missing.comp, baseline[, 1:2]))
        missing.idx <- missing.idx[-c(1:nrow(missing.comp))]
        
        ## add gene pairs which are missing because gene was removed in modified in order to control for presense/absense of mutations
        missing.2 <- !duplicated(rbind(modified[, 1:2], baseline[, 1:2]))
        missing.2 <- missing.2[-c(1:nrow(modified))]
        missing.idx <- missing.2 | missing.idx
        
        ## annotate with character value corresponding to plotting shape
        baseline$missing <- NA
        baseline$missing[missing.idx] <- 16
        
        ## check and see if any new correlations emerge in modified
        added <- modified[modified$q.vals.plotting != 0, 1:2]
        added.idx <- !duplicated(rbind(baseline[, 1:2], added))
        added.idx <- added.idx[-c(1:nrow(baseline))]
        
        ## find index in original dataframe of any discovered correlations
        if (sum(added.idx) > 0){
            added.idx <- rownames(modified) %in% rownames(added[added.idx, ])
            change.idx <- (nrow(baseline) + 1):(nrow(baseline) + sum(added.idx))
            baseline[change.idx, 1:6] <- modified[added.idx, ]
            baseline[change.idx, 7] <- 17
        }
        
        pdf(path, height = 6, width = 7)
        heatmap <- ggplot(data = baseline, aes(x = gene1, y = gene2)) + geom_tile(aes(fill = q.vals.plotting)) + 
            scale_shape_identity() + geom_point(aes(shape = missing, color = factor(missing))) + scale_color_manual(values = c("purple", "orange")) +
            scale_fill_gradient2(limits = c(-30, 30), low = "blue", mid = "white", high = "red") + rameen_theme + theme(panel.background = element_rect(fill = "grey96"))
        print(heatmap)
        dev.off()
        print(heatmap)
    }
}



## manual version for individual comparison

## compute all pairwise comparisons for co-occurence
gene1 <- c()
gene2 <- c()
p.vals <- c()
p.vals.gbm <- c()
or <- c()

gbms <- master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Specific == "Glioblastoma"]
df.correlations.input <- df.wide.glioma[!(df.wide.glioma$V1 %in% gbms), ]

## restrict to a subset of samples
df.correlations.input <- df.correlations.input[df.correlations.input$TP53 == 0, -2]

## remove columns with no entries
df.matrix <- apply(df.correlations.input[, -1], 2, as.numeric)
df.correlations.input <- df.correlations.input[, c(TRUE, as.numeric(colSums(df.matrix)) != 0)]

for (i in 2:(ncol(df.correlations.input) - 2)){
    ## get gene names
    current.name <- colnames(df.correlations.input)[i]
    current <- df.correlations.input[, i]
    current.gbm <- df.correlations.input[df.correlations.input$V1 %in% gbms, i]
    current.non.gbm <- df.correlations.input[!(df.correlations.input$V1 %in% gbms), i]
    targets <- colnames(df.correlations.input)[-c(1:i)]
    gene1 <- c(gene1, rep(current.name, length(targets)))
    gene2 <- c(gene2, targets)
    
    ## calculate p vals across all columns for each phenotype iteratively
    p.vals <- c(p.vals, as.numeric(apply(df.correlations.input[, -c(1:i)], 2, function(x){(fisher.test(table(x, current)))$p.value})))
    or <- c(or, as.numeric(apply(df.correlations.input[, -c(1:i)], 2, function(x){(fisher.test(table(x, current)))$estimate})))
    
    ## alternatively, do two indepdent t-tests based on subset of data
    #p.vals.gbm <- c(p.vals.gbm, exp((log(as.numeric(apply(df.correlations.input[df.correlations.input$V1 %in% gbms, -c(1:i)], 2, function(x){(fisher.test(table(x, current.gbm)))$p.value}))) + 
    #           log(as.numeric(apply(df.correlations.input[!(df.correlations.input$V1 %in% gbms), -c(1:i)], 2, function(x){(fisher.test(table(x, current.non.gbm)))$p.value})))) / 2))
}

## do last one manually
i <- ncol(df.correlations.input) - 1
current.name <- colnames(df.correlations.input)[i]
current <- df.correlations.input[, i]
targets <- colnames(df.correlations.input)[-c(1:i)]
gene1 <- c(gene1, rep(current.name, length(targets)))
gene2 <- c(gene2, targets)
p.vals <- c(p.vals, fisher.test(table(df.correlations.input[, i + 1], current))$p.value)

current.gbm <- df.correlations.input[df.correlations.input$V1 %in% gbms, i]
current.non.gbm <- df.correlations.input[!(df.correlations.input$V1 %in% gbms), i]
p.vals.gbm <- c(p.vals.gbm, exp((log(fisher.test(table(df.correlations.input[df.correlations.input$V1 %in% gbms, i + 1], current.gbm))$p.value) + 
                                     log(fisher.test(table(df.correlations.input[!(df.correlations.input$V1 %in% gbms), i + 1], current.non.gbm))$p.value)) / 2))


or <- c(or, fisher.test(table(df.correlations.input[, i + 1], current))$estimate)
correlations <- data.frame(gene1, gene2, p.vals, or)
correlations$q.vals <- p.adjust(correlations$p.vals, "fdr")

correlations$q.vals.gbm <- p.adjust(correlations$p.vals.gbm, "fdr")


## heatmap for correlations
temp <- correlations
temp$q.vals.plotting <- -log10(temp$q.vals)
temp$q.vals.plotting[temp$q.vals.plotting > 30] <- 30
temp$q.vals.plotting[temp$q.vals.plotting < 2] <- NA
temp$q.vals.plotting[temp$or < 1] <- -temp$q.vals.plotting[temp$or < 1]
temp$gene1 <- factor(temp$gene1, levels = colnames(df.correlations.input)[-1])
temp$gene2 <- factor(temp$gene2, levels = colnames(df.correlations.input)[-1])





## volcano plot of correlations

temp <- correlations
temp$or[temp$or > 32] <- 32
temp$or <- log2(temp$or)
temp$q.vals.plotting <- -log10(temp$q.vals)
temp$q.vals.plotting[temp$q.vals.plotting > 25] <- 25
temp$names <- paste(temp$gene1, temp$gene2, sep = "__")

temp$names[temp$q.vals.plotting < 3 ] <- NA
temp$names[abs(temp$or) < 1] <- NA

temp$names[!(temp$q.vals.plotting > 5 | abs(temp$or) > 2.2)] <- NA

pdf("correlations_volcano_gliomas_all.pdf", width = 20, height = 20)
ggplot(data = temp, aes(x = or, y = q.vals.plotting, label = names)) + geom_point() + geom_text(aes(label=names), size = 3) + scale_y_continuous(limits = c(0, 50)) +
    scale_x_continuous(limits = c(-5, 5)) + theme(panel.background = element_blank())
dev.off()
write.csv(temp, "../Analysis/correlations_raw_values.csv", row.names = FALSE)

## generate mini comuts based on interesting results


mini <- PlotComut(
    mut.maf1 = master.snvs[master.snvs$BEST_EFF_GENE == "IDH1", ],
    focal.cnv.maf = master.cnv.focal[master.cnv.focal$GENE == "EGFR", ],
    broad.cnv.maf = master.cnv.broad[rownames(master.cnv.broad) %in% c("7p"), ],
    samples = master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Broad == "Glioma"],
    input.samples = "SAMPLE_ACCESSION_NBR", 
    input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification", gene.cutoff = 60, 
    dimensions = c(18,10),
    title = "", y.axis.font.size = 6, legend.font.size = 6,
    col.vector = c("frameshift_indel" = "red", "missense" ="skyblue" , "nonsense" = "orange", "splice_site" = "yellow", "stop_codon" = "purple", "in_frame_indel" = "gold",
                   "other" = "cyan", "arm-level gain" = "red", "arm-level loss" = "blue", "2DEL" = "blue", "HA" = "red",
                   "Glioblastoma" = "black", "AnaplasticAstro" = "green", "Angiocentric" = "orange", "Astro" = "green", 
                   "DiffuseAstro" = "green", "Oligo" = "yellow", "OligoAstro" = "yellow", "PilocyticAstro" = "green", 
                   "0" = "gray90", "1" = "gray65", "2" = "gray40", "3" = "gray20", "4" = "gray1"),
    return.plot = TRUE
    #phenotypes = master.phenotypes
    )


## plot multiple plots together

grid.arrange(bp, dp, vp, sc, ncol = 2, nrow = 2)

## nested grid arrange: each grob is treated as a single plot

grid.arrange(bp, arrangeGrob(bp, dp, bp, bp, nrow =2, ncol = 2), ncol =2 )


## to specify exact width and height for each: set width and height of rows, columns
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
            ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))


## plot multiple plots
pdf("Figure 2 combined draft.pdf", width = 15, height = 10)
grid.arrange(vp, arrangeGrob(hm, arrangeGrob(mini, mini, ncol = 1), ncol = 1), ncol = 2)
dev.off()




### permutations testing

## convert to long format for permutations testing
shuffle.wide <- df.correlations.input
shuffle <- melt(shuffle.wide, id.vars = "V1")
colnames(shuffle) <- c("samples", "genes", "mutations")
shuffle$genes <- as.character(shuffle$genes)
shuffle$samples <- as.character(shuffle$samples)

#shuffle$mutations[shuffle$mutations %in% c("Middle", "wt", "Young")] <- 0
#shuffle$mutations[!(shuffle$mutations %in% c("", 0))] <- 1
shuffle <- shuffle[shuffle$mutations == 1, ]
shuffle <- shuffle[shuffle$samples %in% df.wide.glioma$V1, ]


PermuteData <- function(gene1, gene2, reps, data, split = NA){
    ## takes 2 gene names, number of reps of that gene, plus a data frame in long format
    ## if split is provided, splits data into two separate groups which are internally permutated, and the combined distribution is then analyzed
    output <- rep(0, reps)
    if (missing(split)){
        for (i in 1:reps){
            data$genes <- sample(data$genes, nrow(data))
            output[i] <- sum(unique(data[data$genes == gene1, ]$samples) %in% data[data$genes == gene2, ]$samples)
        }
    }else{
        data.1 <- data[data$samples %in% split, ]
        data.2 <- data[!(data$samples %in% split), ]
        for (i in 1:reps){
            data.1$genes <- sample(data.1$genes, nrow(data.1))
            data.2$genes <- sample(data.2$genes, nrow(data.2))
            output[i] <- sum(c(unique(data.1[data.1$genes == gene1, ]$samples), unique(data.2[data.2$genes == gene1, ]$samples)) %in% 
                                 c(data.1[data.1$genes == gene2, ]$samples, data.2[data.2$genes == gene2, ]$samples))
        }
    }
    return(output)
}



## insert permutated p values for top hits
correlations <- correlations[order(correlations$q.vals), ]
correlations$permuted.p.value <- NA

for (i in 1:400){
    expected <- PermuteData(correlations$gene1[i], correlations$gene2[i], 1000, shuffle)
    observed <- sum(unique(shuffle[shuffle$genes == correlations$gene1[i], ]$samples) %in% shuffle[shuffle$genes == correlations$gene2[i], ]$samples)
    correlations$permuted.p.value[i] <- sum(expected >= observed) / 1000
}

correlations$permuted.q.value <- p.adjust(correlations$permuted.p.value, "fdr")
write.csv(correlations, "C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/Analysis/correlations.csv", row.names = F)


data <- shuffle
data$genes <- sample(data$genes, nrow(data))
sum(data$genes == "PIK3C2B-cnv")

## check specific values


check.expected <- PermuteData("7p", "Primary", 10000, shuffle, master.sheet$SAMPLE_ACCESSION_NBR[master.sheet$Cancer_Type_Specific == "Glioblastoma"])
check.observed <- sum(unique(shuffle$samples[shuffle$genes == "7p"]) %in% unique(shuffle$samples[shuffle$genes == "Primary"]))
table(check.expected > check.observed)

shuffle.gbm <- shuffle[shuffle$samples %in% gbms, ]
check.expected.gbm <- PermuteData("7p", "Primary", 10000, shuffle.gbm)
check.observed.gbm <- sum(unique(shuffle.gbm$samples[shuffle.gbm$genes == "7p"]) %in% unique(shuffle.gbm$samples[shuffle.gbm$genes == "Primary"]))
table(check.expected.gbm > check.observed.gbm)

shuffle.lgg <- shuffle[!(shuffle$samples %in% gbms), ]
check.expected.lgg <- PermuteData("7p", "Primary", 10000, shuffle.lgg)
check.observed.lgg <- sum(unique(shuffle.lgg$samples[shuffle.lgg$genes == "7p"]) %in% unique(shuffle.lgg$samples[shuffle.lgg$genes == "Primary"]))
table(check.expected.lgg > check.observed.lgg)

fisher.test(df.wide.glioma[!(df.wide.glioma$V1 %in% gbms), ]$`7p`, df.wide.glioma[!(df.wide.glioma$V1 %in% gbms), ]$'Primary')
fisher.test(df.wide.glioma[df.wide.glioma$V1 %in% gbms, ]$`7p`, df.wide.glioma[df.wide.glioma$V1 %in% gbms, ]$'Primary')
fisher.test(df.wide.glioma$`7p`, df.wide.glioma$'Primary')






