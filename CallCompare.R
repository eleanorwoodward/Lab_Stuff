## Noah Greenwald
## Compare mutations found in same sample via different calling methods, filtering, or other criteria
## Adapted from R Cookbook

source("C:/Users/Noah/OneDrive/Work/R/Scripts/MafFunctions.R")
source("C:/Users/Noah/OneDrive/Work/R/Scripts/ExacFilter.R")
source("C:/Users/Noah/OneDrive/Work/R/Scripts/ESPFilter.R")
source("C:/Users/Noah/OneDrive/Work/R/Scripts/PonFilter.R")

## Gets all unpaired samples
folder <- ("C:/Users/Noah/OneDrive/Work/Meningioma/Analysis/PoN")
up.snp <- CombineMaf(folder)

paired.list <- c("M2-tumor", "M17-tumor", "M26-tumor", "M44-tumor", "M45-tumor", "M133-tumor", "M148-tumor", "M159-tumor", "M203-tumor", "M224-tumor", 
                 "M226-tumor", "M255-tumor", "M267-tumor", "M274-tumor", "M276-tumor")

## Gets paired samples
p.snp <- FilterMaf(val.snp.no.filter, paired.list, "Tumor_Sample_Barcode")

if (sum(duplicated(up.snp)) | sum(duplicated(p.snp)) > 0){
  print("There are duplicate rows in the data frames")
}

## Make tumor sample names the same for comparison
up.snp$Tumor_Sample_Barcode <- sapply(up.snp$Tumor_Sample_Barcode, PairSetFormat, 5, USE.NAMES = FALSE)

## Values to run optimization over
exac.afs <- c(.000005, .00005, .0005, .005, .01, .05)
pon.logs <- c(-3.0, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1.0)

## Sets empty list to hold outputs
pairs <- c()
cutoff <- vector()
coding.variants <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Nonstop_Mutation", 
                     "De_Novo_Start_OutOfFrame")

for (i in 1:length(pon.logs)){
    
    #up.snp <- run.exac(up.snp, .05)
    #up.snp <- run.esp(up.snp)
    up.snp <- run.pon(up.snp, pon.logs[i])
    
    ## Create indices to use to find filtered mutations for given combinations
    exac.germs <- up.snp$germline == TRUE
    esp.germs <- up.snp$esp_germline == TRUE
    pon.germs <- (up.snp$pon_germline == TRUE)
    all.germs <- (up.snp$esp_germline == TRUE | up.snp$pon_germline == TRUE)
    
    ##Use correct index to generate plots here
    up.snp.muts <- up.snp[!pon.germs, ]
    up.snp.germ <- up.snp[pon.germs, ]

    ## Elminate unnecessary/potentially non-homologous columns so that duplicates can be found
    keepers <- c("Hugo_Symbol", "Chromosome", "Start_position", "Variant_Classification", "i_tumor_f", "Tumor_Sample_Barcode")

    ## Original paired sample
    small.p.snp <- p.snp[, keepers]

    ## Germline calls from unpaired
    small.up.snp.germ <- up.snp.germ[, keepers]

    ## Non-germline calls from unpaired
    small.up.snp.muts <- up.snp.muts[, keepers]
    
    small.up.snp <- up.snp[, keepers]
    
    ## Combine databses together to find duplicates
    somatic.flagged <- duplicated(rbind(small.up.snp.germ, small.p.snp))
    somatic.called <- duplicated(rbind(small.p.snp, small.up.snp))

    GermlineRemoved <- function(called.germline, total.muts, actual.somatic, false.flagged.somatic){
      pct <- (called.germline - false.flagged.somatic) / (total.muts - actual.somatic)
      round(pct, 3) * 100
    }
    
    SomaticRemoved <- function (false.flagged.somatic, actual.somatic){
      pct <- false.flagged.somatic / actual.somatic
      round(pct, 3) * 100
    }
    
    ## Get number of mutations for analysis
    n.up.total <- nrow(up.snp)
    n.up.snp.germ <- nrow(up.snp.germ)
    n.p.snp <- nrow(p.snp)
    
    paste("These samples had a total of " , n.up.total, " mutations detected by MuTect. Of these, ExAc flagged ", 
         n.up.snp.germ, "as germline, missing ",
         100 - GermlineRemoved(n.up.snp.germ, n.up.total, n.p.snp, sum(somatic.flagged)),
          "% of the germline events. There were a total of ", n.p.snp, 
          " somatic mutations. Of these, exac flagged ", sum(somatic.flagged),
          " as germline events, incorrectly removing ", SomaticRemoved(sum(somatic.flagged), n.p.snp), 
          "% of the true mutations. This was with an allelic fraction cutoff of ", pon.logs[i])
    
    ## Store values for barplot
    pairs <- c(pairs, 100 - GermlineRemoved(n.up.snp.germ, n.up.total, sum(somatic.called), sum(somatic.flagged)), 
              SomaticRemoved(sum(somatic.flagged), sum(somatic.called)))
    cutoff <- c(cutoff, pon.logs[i])
}

error.matrix <- matrix(pairs, 2, length(pairs)/2)

barplot(error.matrix, beside = TRUE, names.arg = cutoff, xlab = "Log Likelihood Cutoff Value", ylab = "Percent" ,
        legend.text = c("Germline Mutations Missed", "Somatic Mutations Removed"), 
        main = "False negatives vs false positives with PoN likelihood filtering", args.legend = list(x = "topleft"))


## For multiple different filters
combined <- rbind(error.matrix1, error.matrix2)

barplot(combined, beside = TRUE, names.arg = cutoff, xlab = "Log Likelihood Cutoff Value", ylab = "Percent" ,
        legend.text = c("% Germline Mutations Missed", "% Somatic Mutations Removed", "with exac", "with exac"), 
        main = "False negatives vs false positives with Log Likelihood Filtering", args.legend = list(x = "topleft"))


