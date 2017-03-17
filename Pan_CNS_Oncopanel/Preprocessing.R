## Data cleaning and pre-processing

## Generate master table for all samples from various data sources
setwd("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/OncDRS data/")

master.sheet <- read.csv("REQ_ID08_65337_ONCOPANEL_SPECIMEN.csv", stringsAsFactors = F)
oncomap <- read.csv("REQ_ID08_65337_ONCOMAP_MUT_RESULTS_PROFILE.csv", stringsAsFactors = F)
oncomap$na_column <- NA
cancer.diag <- read.csv("REQ_ID08_65337_CANCER_DIAGNOSIS_CAREG.csv", stringsAsFactors = F)
identifiers <- read.csv("REQ_ID08_65337_PT_INFO_STATUS_REGISTRATION.csv", stringsAsFactors = F)
master.sheet.short <- master.sheet[, c("PATIENT_ID", "SAMPLE_ACCESSION_NBR", "PRIMARY_CANCER_DIAGNOSIS", "ORIGINAL_PATH_DIAGNOSIS", 
                                       "BIOPSY_SITE", "BIOPSY_SITE_TYPE", "TUMOR_PURITY", "PANEL_VERSION", "REPORT_DT", "REPORT_COMMENT", "SNV_COUNT")]

## combine oncopanel with oncomap information
oncomap.unique <- oncomap[!duplicated(oncomap$PATIENT_ID), ]
oncomap.short <- oncomap.unique[, c("PATIENT_ID", "MOLE_SPEC_ACCESSION_NBR", "SPECIMEN_HISTOLOGY_NM", "PATH_DIAGNOSIS_NM",
                             "SPECIMEN_ANATOMIC_SITE_NM", "SPECIMEN_TYPE_NM", "na_column", "na_column", "CALENDAR_DT", "na_column", "na_column"), ]

colnames(oncomap.short) <- colnames(master.sheet.short)
oncomap.short$TUMOR_PURITY <- NA
oncomap.short$PANEL_VERSION <- 0
oncomap.short$REPORT_COMMENT <- NA
oncomap.short$SNV_COUNT <- NA

master.sheet.short <- rbind(master.sheet.short, oncomap.short)


## convert from current date format to numeric 
temp <- sapply(master.sheet.short$REPORT_DT, as.Date, "%d-%b-%y", USE.NAMES = F)

## can't use SAPPLY, simplifies format
DATE <- rep(as.Date(temp[1], origin = "1970-01-01"), length(temp))
for (i in 1:length(temp)){
    DATE[i] <- as.Date(temp[i], origin = "1970-01-01" )
}
master.sheet.short$DATE <- DATE

## marks which samples will have CNV data from oncopanel
master.sheet.short$CNV_ONC <- DATE > "2014-04-10"

## marks which patients have multiple tumors
master.sheet.short$multiple_samples <- master.sheet.short$PATIENT_ID %in% master.sheet.short[(duplicated(master.sheet.short$PATIENT_ID)),]$PATIENT_ID 

## checks which patients have multiple cancer diagnosis
cancer.diag$multiple_samples <- cancer.diag$PATIENT_ID %in% cancer.diag[(duplicated(cancer.diag$PATIENT_ID)),]$PATIENT_ID 


## import information from cancer diagnosis sheet
master.sheet.short$Cancer_Diagnosis_Detailed <- NA
for (i in 1:nrow(master.sheet.short)){
    ## copies information over from those cases present in detailed cancer diagnosis spreadsheet
    if (master.sheet.short$PATIENT_ID[i] %in% cancer.diag$PATIENT_ID){
        ## indicates duplicates exist those samples that have duplicates
        if (master.sheet.short[i, "multiple_samples"] == T | cancer.diag[cancer.diag$PATIENT_ID == master.sheet.short$PATIENT_ID[i], ]$multiple_samples[1] == T){
            master.sheet.short$Cancer_Diagnosis_Detailed[i] <- "Multiple samples"
        }else{
            master.sheet.short$Cancer_Diagnosis_Detailed[i] <- cancer.diag[cancer.diag$PATIENT_ID == master.sheet.short$PATIENT_ID[i], ]$HISTOLOGY_DESCR
        }
    }
}

master.sheet.short$MRN <- NA
for (i in 1:nrow(master.sheet.short)){
    id <- master.sheet.short$PATIENT_ID[i]
    mrn <- identifiers[identifiers$PATIENT_ID == id, ]$BWH_MRN
    master.sheet.short$MRN[i] <- mrn
}


## TODO: add code for processing aCGH data



## write cleaned up table for manual review
master.sheet.ordered <- master.sheet.short[, c("PATIENT_ID", "MRN", "SAMPLE_ACCESSION_NBR", "PRIMARY_CANCER_DIAGNOSIS", "ORIGINAL_PATH_DIAGNOSIS", 
                                               "Cancer_Diagnosis_Detailed", colnames(master.sheet.short)[6:14])]

## remove /t symbols to enable effecient writing to excel
master.sheet.ordered$REPORT_COMMENT<- sapply(1:nrow(master.sheet.ordered), function(x){gsub("\t", "", master.sheet.ordered$REPORT_COMMENT[x])})


write.table(master.sheet.ordered, "../Analysis/master.sheet.tsv", row.names = F, sep = "\t")

