## Merge data 

Master <- read.delim("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/Master_Sheet_For_R.txt", stringsAsFactors = F)

PDL <- read.delim("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/PDL_for_R.txt", stringsAsFactors = F)

Master[, 6:13] <- "no data"
for (i in 1:nrow(Master)){
    idx <- PDL$mrn %in% Master$mrn[i]
    if (sum(idx) == 1){
        row <- PDL[idx,4:11]
        Master[i, 6:13] <- row
        
    }else if (sum(idx) > 1){
        Master[i, ] <- "Duplicates"
    }
}

write.csv(Master, "C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/R_ouput.csv", row.names = F)

total.data <- read.delim("C:/Users/Noah/Syncplicity Folders/Radiomics-NG LB only/", stringsAsFactors = F)

