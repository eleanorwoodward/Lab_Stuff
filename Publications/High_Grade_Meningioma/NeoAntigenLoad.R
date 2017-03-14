## Analysis of neo-antigen load

neo.snp.folder <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Neoantigen_Load/snv"

neo.snps <- NULL
for (i in 1:length(list.files(neo.snp.folder, pattern = "mutation"))){
    temp <- read.delim(paste(neo.snp.folder, list.files(neo.snp.folder, pattern = "mutation")[i], sep = "/"),
                       stringsAsFactors = F, comment.char = "#")
    name <- strsplit(list.files(neo.snp.folder, pattern = "mutation")[i], "_")[[1]][1]
    temp <- cbind(temp, rep(name, nrow(temp)))
    neo.snps <- rbind(neo.snps, temp)    
}


neo.indel.folder <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Neoantigen_Load/indel/"

neo.indels <- NULL
for (i in 1:length(list.files(neo.indel.folder, pattern = "mutation"))){
    temp <- read.delim(paste(neo.indel.folder, list.files(neo.indel.folder, pattern = "mutation")[i], sep = "/"),
                       stringsAsFactors = F, comment.char = "#")
    name <- strsplit(list.files(neo.indel.folder, pattern = "mutation")[i], "_")[[1]][1]
    temp <- cbind(temp, rep(name, nrow(temp)))
    neo.indels <- rbind(neo.indels, temp)    
}

write.table(neo.snps, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Neoantigen_Load/total.snps.txt", 
            row.names = F, quote = F, sep = "\t")
write.table(neo.indels, "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Neoantigen_Load/total.indels.txt", 
            row.names = F, quote = F, sep = "\t")

neo.snps <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Neoantigen_Load/total.snps.txt", stringsAsFactors = F)
neo.indels <- read.delim("C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/Neoantigen_Load/total.indels.txt", stringsAsFactors = F)

small.snps <- cbind(neo.snps[, c("chr", "pos", "gene.symbol", "mut_median_IC50", "rep.name..nrow.temp..")], rep("snp", nrow(neo.snps)))
colnames(small.snps) <- c("chr", "pos", "gene.symbol", "median_IC50", "sample", "variant")
small.indels <- cbind(neo.indels [c("chr", "start", "gene_symbol", "median_IC50", "rep.name..nrow.temp..")], rep("indel", nrow(neo.indels)))
colnames(small.indels) <- c("chr", "pos", "gene.symbol", "median_IC50", "sample", "variant")

neo.snindels <- rbind(small.indels, small.snps)


## load ccf data across entire meningioma cohort
neo.absolute.folder <- "C:/Users/Noah/Syncplicity Folders/Meningioma (Linda Bi)/ABSOLUTE/all.all.mafs/"
neo.ccfs <- NULL
for (i in 1:length(list.files(neo.absolute.folder))){
    temp <- readRDS(paste(neo.absolute.folder, list.files(neo.absolute.folder)[i], sep = "/"))
    neo.ccfs <- rbind(neo.ccfs, temp[, c("pair_id", "Hugo_Symbol", "Chromosome", "Start_position", "Variant_Classification", "i_tumor_f", "ccf_hat")])
}

neo.ccfs <- neo.ccfs[neo.ccfs$Variant_Classification %in% c(snp.variants, indel.variants), ]
neo.ccfs <- neo.ccfs[!is.na(neo.ccfs$i_tumor_f), ]


neoantigen.combined <- neo.ccfs
neoantigen.combined$ic50 <- NA
#neoantigen.combined <- neoantigen.combined[!(neoantigen.combined$pair_id %in% master.table[master.table$Sequencing == "WGS", ]$Pair.Name), ]
neoantigen.combined[neoantigen.combined$pair_id == "MEN0122_RBM-P1", ]$pair_id <- "MEN0122-RBM-P1"
neoantigen.combined[neoantigen.combined$pair_id == "MEN0122_RBM-P2", ]$pair_id <- "MEN0122-RBM-P2"
neoantigen.combined[neoantigen.combined$pair_id == "MEN0122_RBM-P3", ]$pair_id <- "MEN0122-RBM-P3"
neoantigen.combined[neoantigen.combined$pair_id == "MEN0122_RBM-P4", ]$pair_id <- "MEN0122-RBM-P4"
neoantigen.combined[neoantigen.combined$pair_id == "MEN0122_RBM-P5", ]$pair_id <- "MEN0122-RBM-P5"

for (i in 1:nrow(neoantigen.combined)){
    sample <- neoantigen.combined$pair_id[i]
    chr <- neoantigen.combined$Chromosome[i]
    pos <- neoantigen.combined$Start_position[i]
    
    hit <- neo.snindels[neo.snindels$sample == sample & neo.snindels$chr == chr & neo.snindels$pos == pos, ]
    if (nrow(hit) == 0){
        #nothing
    }else if (nrow(hit) == 1){
        neoantigen.combined$ic50[i] <- hit$median_IC50
    }else{
        neoantigen.combined$ic50[i] <- min(hit$median_IC50)
    }
}

## missing values are all frameshifts or start before coding sequence of the protein

unique(neoantigen.combined[is.na(neoantigen.combined$ic50), ]$pair_id)
unique(neoantigen.combined[!(is.na(neoantigen.combined$ic50)), ]$pair_id)

table(is.na(neoantigen.combined$ic50))

clean.neoantigen <- neoantigen.combined[!is.na(neoantigen.combined$ic50), ]

clean.neoantigen[clean.neoantigen$pair_id == "MEN0122-RBM-P1", ]$pair_id <- "MEN0122_RBM-P1"
clean.neoantigen[clean.neoantigen$pair_id == "MEN0122-RBM-P2", ]$pair_id <- "MEN0122_RBM-P2"
clean.neoantigen[clean.neoantigen$pair_id == "MEN0122-RBM-P3", ]$pair_id <- "MEN0122_RBM-P3"
clean.neoantigen[clean.neoantigen$pair_id == "MEN0122-RBM-P4", ]$pair_id <- "MEN0122_RBM-P4"
clean.neoantigen[clean.neoantigen$pair_id == "MEN0122-RBM-P5", ]$pair_id <- "MEN0122_RBM-P5"


plot(log2(clean.neoantigen$ic50), clean.neoantigen$ccf_hat)

## high CCF mutations more likely to be immunogenic

immunogenic <- clean.neoantigen[clean.neoantigen$Variant_Classification %in% snp.variants, ]$ic50 < 500
clonal <- clean.neoantigen[clean.neoantigen$Variant_Classification %in% snp.variants, ]$ccf_hat == 1
fisher.test(table(immunogenic, clonal))
371/(1441+371)
25/(199+25)

clonal.neoantigen <- clean.neoantigen[clean.neoantigen$ccf_hat == 1, ]
mean(sort(table(clonal.neoantigen$pair_id)))

## Force call neo-antigen mutations
neo.forecall.list <- list(c("MEN0042.TumorA", "MEN0042.TumorB", "MEN0042.TumorC"), 
                   c("MEN0045.TumorA", "MEN0045.TumorB", "MEN0045.TumorC", "MEN0045.TumorD", "MEN0045.TumorE"), c("MEN0048.TumorA", "MEN0048.TumorD", "MEN0048.TumorB", "MEN0048.TumorC"),
                   c("MEN0093.TumorC", "MEN0093.TumorD", "MEN0093.TumorA", "MEN0093.TumorE"), c("MEN0097.Tumor", "MEN0097.TumorA", "MEN0097.TumorB", "MEN0097.TumorC"),
                   c("MEN0101.TumorB", "MEN0101.Tumor"), c("MEN0118.TumorA", "MEN0118.TumorB"), c("MEN0119.TumorA", "MEN0119.TumorB"), c("MEN0120.Tumor", "MEN0120.TumorB"), 
                   c("MEN0121.TumorA", "MEN0121.TumorB", "MEN0121.TumorC", "MEN0121.TumorD"), c("MEN0122_RBM.TumorA", "MEN0122_RBM.TumorB", "MEN0122_RBM.TumorC",
                                                                                                "MEN0122_RBM.TumorD", "MEN0122_RBM.TumorE"))
neo.mafs <- list()
for (i in 1:length(neo.forecall.list)){
    pairs <- c()
    for (j in 1:length(neo.forecall.list[[i]])){
        ## sets up relevant info for each sample
        sample.name <- neo.forecall.list[[i]][j]
        pairs <- c(pairs, master.table[master.table$Tumor.Name == sample.name, ]$Pair.Name)
    }
    temp.maf <- FilterMaf(clean.neoantigen, pairs,"pair_id")
    forcefull <- ForceCallMutations(temp.maf, "pair_id","Hugo_Symbol", "Start_position")
    neo.mafs[[i]] <- forcefull
}
    
    
## plot mutations per patient
ggplot(data=data.frame(sort(table(neoantigen.combined$pair_id), decreasing = T)), aes(x = Var1, y = Freq)) +
geom_bar(stat = "identity") + rameen_theme


neoantigen.combined$pair_id <- factor(neoantigen.combined$pair_id, levels = names(sort(table(neoantigen.combined$pair_id), decreasing = T)))
ggplot(data=neoantigen.combined, aes(x = pair_id)) +
    geom_bar() + rameen_theme

## Classify into low and high subgroups for heterogeneity





## Creates stacked bar plot with private, shared, ubiquitous neo-antigens

## Figs
## create stacked bar plot for each patient with totally private SCNAS, totally shared SCNAS, and anything in between 
sample.vector <- c()
clonal.vector <- c()
clonal.muts <- list()

for (i in 1:length(neo.mafs)){
    maf <- neo.mafs[[i]]
    maf <- apply(maf[, -(1:2)],2, as.numeric)
    clonal.vals <- rowSums(maf)
    samples <- ncol(maf)
    clonal.vals <- clonal.vals[clonal.vals != 0]
    current.clonal <- rep("private", length(clonal.vals))
    current.clonal[clonal.vals == samples] <- "clonal"
    current.clonal[clonal.vals < samples & clonal.vals > 1] <- "subclonal"
    clonal.vector <- c(clonal.vector, current.clonal)
    sample.vector <- c(sample.vector, rep(colnames(maf)[1], length(clonal.vals)))
}

clonal.df <- data.frame(sample.vector, clonal.vector)
clonal.df$sample.vector <- factor(clonal.df$sample.vector, levels=names(sort(table(clonal.df$sample.vector), decreasing = T)))
clonal.df$clonal.vector <- factor(clonal.df$clonal.vector, levels = c("clonal", "subclonal", "private"))

ggplot(data=clonal.df, aes(x=sample.vector, fill=clonal.vector)) + geom_bar() + scale_fill_hue(l = "45") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + rameen_theme


## plot percent immunogenic alterations in primary vs recurrence

vals <- master.table[master.table$Pair.Name %in% c(as.character(unique(neo.snindels$sample)), "MEN0122_RBM-P1", "MEN0122_RBM-P2","MEN0122_RBM-P3",
                                                     "MEN0122_RBM-P4", "MEN0122_RBM-P5"), ]$percent.immunogenic
#keepers.idx <- c(1, 3:4, 6:9, 11:13, 18:20, 24:26, 31, 49, 51,53,55:57, 59:62)
keepers.idx <- c(1, 3:4, 11:13, 18:20, 31, 49, 51,53,55:57, 59:62)

vals.x <- vals[keepers.idx]
vals.y <- vals[keepers.idx + 1]

plot(vals.x, vals.y, ylim = c(0, 1), xlim = c(0, 1))
# Add fit lines
abline(lm(vals.y~vals.x), col="red") # regression line (y~x) 

mut.model <- lm(vals.y ~ vals.x)
summary(mut.model)


mean(PatientAverage(neo.mafs[[5]], 2))

men97 <- clean.neoantigen[grep("MEN0097", clean.neoantigen$pair_id), ]
men97.clonal <- men97[men97$ccf_hat > .8, ]
men97.subclonal <- men97[men97$ccf_hat < .8, ]

men97.forecalled <- neo.mafs[[5]]
men97.forecalled <- cbind(men97.forecalled, rowSums(apply(men97.forecalled[, -(1:2)], 2, as.numeric)))
colnames(men97.forecalled)[7] <- "total"

fisher.test(table(men97.forecalled$Hugo_Symbol %in% men97.clonal$Hugo_Symbol, men97.forecalled$total > 1))

## Revision
## check for recurrence of specific neo-antigens
## remove non-unique patients
bye.bye <- c("MEN0030-P1", "MEN0042-P2", "MEN0042G-P1","MEN0045-P1", "MEN0045-P2", "MEN0045-P3", "MEN0045-P4", "MEN0048-P2", "MEN0048-P3", "MEN0048G-P1", 
             "MEN0093-P1", "MEN0093-P4", "MEN0093G-P2", "MEN0097-P1", "MEN0097-P2", "MEN0097-P3", "MEN0101-P1", "MEN0118-P1", "MEN0119-P1", "MEN0120-P1",
             "MEN0121-P1", "MEN0121-P2", "MEN0121-P3", "MEN0121-P4", "MEN0122-RBM-P2", "MEN0122-RBM-P3", "MEN0122-RBM-P4", "MEN0122-RBM-P5")
neo.snps.unique <- neo.snps[!(neo.snps$rep.name..nrow.temp.. %in% bye.bye), ]
hits <- neo.snps.unique[duplicated(neo.snps.unique$mut_peptide), ]$mut_peptide
View(neo.snps[neo.snps$mut_peptide %in% hits, ])
table(duplicated(neo.snps.unique$mut_peptide))
52/(52+ 2715)


