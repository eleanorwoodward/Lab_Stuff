## Noah Greenwald
## Compare mutations found in same sample via different calling methods, filtering, or other criteria
## Adapted from R Cookbook

run.compare <- function(vals, filter){
  if(filter %in% c("PoN", "Exac", "ESP") == FALSE){
    stop("Unsupported filter specified")
  }
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


## Sets empty list to hold outputs
pairs <- c()
cutoff <- vector()
coding.variants <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Nonstop_Mutation", 
                     "De_Novo_Start_OutOfFrame")

for (i in 1:length(vals)){
    if (filter == "PoN"){
      up.snp <- run.pon(up.snp, vals[i])
      pon.germs <- (up.snp$pon_germline == TRUE)
      up.snp.muts <- up.snp[!pon.germs, ]
      up.snp.germ <- up.snp[pon.germs, ]
    }else if (filter == "Exac"){
      up.snp <- run.exac(up.snp, vals[i])
      exac.germs <- up.snp$germline == TRUE
      up.snp.muts <- up.snp[!exac.germs, ]
      up.snp.germ <- up.snp[exac.germs, ]
    }else{
      up.snp <- run.esp(up.snp)
      esp.germs <- up.snp$esp_germline == TRUE
      up.snp.muts <- up.snp[!esp.germs, ]
      up.snp.germ <- up.snp[esp.germs, ]
    }
    
   
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
    real.muts <-rbind(small.p.snp, small.up.snp)
    somatic.called <- duplicated(real.muts)
    present <- real.muts[somatic.called, ]
    comb <- rbind(present, small.p.snp)
    idx <- duplicated(comb)
    missing <- comb[!idx, ]
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
    
    print(paste("These samples had a total of " , n.up.total, " mutations detected by MuTect. Of these, ExAc flagged ", 
         n.up.snp.germ, "as germline, missing ",
         100 - GermlineRemoved(n.up.snp.germ, n.up.total, n.p.snp, sum(somatic.flagged)),
          "% of the germline events. There were a total of ", n.p.snp, 
          " somatic mutations, though only", sum(somatic.called), " were seen. Of these, exac flagged ", sum(somatic.flagged),
          " as germline events, incorrectly removing ", SomaticRemoved(sum(somatic.flagged), sum(somatic.called)), 
          "% of the true mutations. This was with an allelic fraction cutoff of ", pon.logs[i]))
    
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
# combined <- rbind(error.matrix1, error.matrix2)
# 
# barplot(combined, beside = TRUE, names.arg = cutoff, xlab = "Log Likelihood Cutoff Value", ylab = "Percent" ,
#         legend.text = c("% Germline Mutations Missed", "% Somatic Mutations Removed", "with exac", "with exac"), 
#         main = "False negatives vs false positives with Log Likelihood Filtering", args.legend = list(x = "topleft"))
# 

}
