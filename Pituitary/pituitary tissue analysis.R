## formatting for pituitary tumor spreadsheet. Makes sure that all samples with paired blood, normals, etc, are properly matched. 

master.sheet <- read.delim("C:/Users/Noah/Syncplicity Folders/Pit (Linda Bi)/WGS project/Tissue Bank/Combined Pituitary Cohort 2017.03.06.txt",
                           stringsAsFactors = F, header = T)

tissue.bank <- read.delim("C:/Users/Noah/Syncplicity Folders/Pit (Linda Bi)/WGS project/Tissue Bank/All cases.txt",
                          stringsAsFactors = F, header = T)

master.sheet[, 13:14] <- 0
colnames(master.sheet)[13:14] <- c("Tissue.Available", "Matched.Blood")

## goes through each entry in master sheet, checking if tissue or blood is banked
for (i in 1:nrow(master.sheet)){
    ## gets all entries associated with MRN
    samples <- FilterMaf(tissue.bank, master.sheet$BWH_MRN[i], "MRN")
    if (nrow(samples) > 0){
        idx <- samples$SampleType %in% "Tissue"
        
        if (sum(idx) > 0){
            if (sum(idx) > 1){
                master.sheet[i, 13] <- "multiple"
            }else{
                master.sheet[i, 13] <- "Yes"
            }
        }
            
            
        blood <- samples$SampleType %in% c("DNA", "Buffy Coat", "Plasma")
        if (sum(blood) > 0){
            if (sum(blood) > 1){
                master.sheet[i, 14] <- "multiple"
            }else{
                master.sheet[i, 14] <- "Yes"
            }            
        }
    
    }
}

## calculate which patients have both blood and tissue available
master.sheet$paired <- master.sheet$Tissue.Available != "0" & master.sheet$Matched.Blood != "0"

## make sure we're not missing any cases
table(unique(tissue.bank$MRN) %in% master.sheet$BWH_MRN)
sum(master.sheet$Tissue.Available != "0" | master.sheet$Matched.Blood != "0")

## investigate potential cases
new.cases <- tissue.bank[!(tissue.bank$MRN %in% master.sheet$BWH_MRN), ]
tissue.mrn <- unique(new.cases[new.cases$SampleType == "Tissue", ]$MRN)
blood.mrn <- unique(new.cases[new.cases$SampleType %in% c("DNA", "Buffy Coat", "Plasma"), ]$MRN)
combined.mrn <- c(tissue.mrn, blood.mrn)
write.csv(combined.mrn[duplicated(combined.mrn)], "C:/Users/Noah/Syncplicity Folders/Pit (Linda Bi)/WGS project/Tissue Bank/new_cases.csv")


write.csv(master.sheet, file = "C:/Users/Noah/Syncplicity Folders/Pit (Linda Bi)/WGS project/Tissue Bank/master_sheet_annotated.csv", row.names = F)
