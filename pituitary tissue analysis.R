## formatting for pituitary tumor spreadsheet. Makes sure that all samples with paired blood, normals, etc, are properly matched. 

master.sheet <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/Master_Sheet_R.txt",
                           stringsAsFactors = F, header = T)

tissue.bank <- read.delim("C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/Tissue_Bank_inventory_R.txt",
                          stringsAsFactors = F, header = T)


output.sheet <- master.sheet[, c(1,2,3,4)]
output.sheet[, 5:8] <- 0
colnames(output.sheet)[5:8] <- c("Tissue.Available", "Bank.Number", "Matched.Blood", "Bank.Number.2")

for (i in 1:nrow(output.sheet)){
    samples <- FilterMaf(tissue.bank, output.sheet$BWH_MRN[i], "MRN")
    if (nrow(samples) > 0){
        idx <- samples$SampleType %in% "Tissue"
        
        if (sum(idx) > 0){
            if (sum(idx) > 1){
                output.sheet[i, 5:6] <- c("multiple", "multiple")
            }else{
                output.sheet[i, 5:6] <- c(samples$Date[idx], samples$SampleSetID[idx])
            }
            
            blood <- samples$SampleType %in% c("DNA", "Buffy Coat")
            if (sum(blood) > 0){
                blood.type <- c()
                blood.label <- c()
                for (j in 1:sum(blood)){
                    blood.type <- paste(blood.type, samples$SampleType[blood][j], sep = ",")
                    blood.label <- paste(blood.label, samples$SampleSetID[blood][j], sep = ",") 
                }
                    output.sheet[i, 7:8] <- c(blood.type, blood.label)
            }
        }
    }
}



write.csv(output.sheet, file = "C:/Users/Noah/OneDrive/Work/Pituitary Tumor Paper/tissue_banking.csv", row.names = F)
