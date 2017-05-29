## Data cleaning and pre-processing

## latest iteration of master sheet
master.sheet <- read.delim("../OncDRS data/Master Sheet 20170524.txt", stringsAsFactors = F)
master.sheet <- merge(master.sheet, identifiers[, c("CHILDRENS_MRN", "PATIENT_ID"),], all.x = TRUE)
write.table(master.sheet, "../OncDRS data/master master 20170524.tsv", row.names = FALSE, sep = "\t")
## update date of death
ids <- master.sheet$PATIENT_ID[master.sheet$Date_of_death == "can't find"]

death.info <- identifiers[identifiers$PATIENT_ID %in% ids, c("PATIENT_ID", "DERIVED_DEATH_DT")]
death.info$DERIVED_DEATH_DT_fixed <- as.Date(death.info$DERIVED_DEATH_DT, "%d-%b-%y")

colnames(death.info)[3] <- "DoD"
master.sheet$DoD <- master.sheet$Date_of_death

## censored values: double check once we have all dates that this is correct
master.sheet$DoD <- as.Date(master.sheet$DoD, "%m/%d/%Y")
master.sheet <- merge(master.sheet, death.info[, -2], "PATIENT_ID", all.x = TRUE)
master.sheet$DoD.x[!is.na(master.sheet$DoD.y)] <- master.sheet$DoD.y[!is.na(master.sheet$DoD.y)]
colnames(master.sheet)[32] <- "DoD"
master.sheet <- master.sheet[, -33]

write.table(master.sheet, "../OncDRS data/Master Sheet 20170512.tsv", sep = "\t", row.names = F)
## Generate master table for all samples from various data sources
setwd("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/OncDRS data/")

master.sheet <- read.csv("REQ_ID08_65337_ONCOPANEL_SPECIMEN.csv", stringsAsFactors = F)
oncomap <- read.csv("REQ_ID08_65337_ONCOMAP_MUT_RESULTS_PROFILE.csv", stringsAsFactors = F)
oncomap$na_column <- NA
cancer.diag <- read.csv("REQ_ID08_65337_CANCER_DIAGNOSIS_CAREG.csv", stringsAsFactors = F)
identifiers <- read.csv("../OncDRS data/REQ_ID08_65337_PT_INFO_STATUS_REGISTRATION.csv", stringsAsFactors = F)
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


## second round cleaning, after manual annotation of pathology subtypes

master.sheet <- read.delim("../OncDRS data/Master_Sheet_R_Upload 20170502.txt", stringsAsFactors = F)
master.sheet[master.sheet$MRN %in% master.sheet$MRN[duplicated(master.sheet$MRN)], ]$multiple_samples <- TRUE

## add block number
opanel.sheet <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_SPECIMEN.csv", stringsAsFactors = F)
opanel.sheet <- opanel.sheet[, 4:5]
master.sheet <- merge(master.sheet, opanel.sheet, "SAMPLE_ACCESSION_NBR", all.x = T)
colnames(master.sheet)[which(colnames(master.sheet) == "BLOCK_ACCESSION_NBR.x")] <- "BLOCK_ACCESSION_NBR"
master.sheet <- master.sheet[, -(which(colnames(master.sheet) == "BLOCK_ACCESSION_NBR.y"))]

## remove trailing block number
master.sheet$BLOCK_ACCESSION_NBR <- sapply(master.sheet$BLOCK_ACCESSION_NBR, function(x){gsub("-[A-Z][0-9]$", "", x)}, USE.NAMES = FALSE)
master.sheet$BLOCK_ACCESSION_NBR <- sapply(master.sheet$BLOCK_ACCESSION_NBR, function(x){gsub(" [A-Z][0-9]$", "", x)}, USE.NAMES = FALSE)


## read in annotated version of other master sheets to add meningioma/pituitary/glioma pathology information
pit.data <- read.delim("../OncDRS data/Pituitary_R_Upload_20170501.txt", stringsAsFactors = F)
pit.data <- pit.data[, c("Surg.Path..", "BWH_MRN", "Recurrent.", "Histochemical.Expression", "Disrupted")]
colnames(pit.data) <- c("BLOCK_ACCESSION_NBR", "MRN", "Primary", "Cancer_Type_Specific", "Cancer_Type_Details")
pit.data$Primary <- !pit.data$Primary
master.sheet <- merge(master.sheet, pit.data, "MRN", all.x = TRUE)

## checks to see which ones disagree, remove data from those. 
master.sheet[which(master.sheet$BLOCK_ACCESSION_NBR.y != master.sheet$BLOCK_ACCESSION_NBR.x), c("exclude")] <-  "check"

## because merge doesn't work when rows have different values, have to combine values for common columns, even if NA's as placeholders
master.sheet$Primary.x[!is.na(master.sheet$Primary.y) & !master.sheet$multiple_samples] <- master.sheet$Primary.y[!is.na(master.sheet$Primary.y) & !master.sheet$multiple_samples]
master.sheet$Cancer_Type_Specific.x[!is.na(master.sheet$Cancer_Type_Specific.y) & !master.sheet$multiple_samples] <- 
    master.sheet$Cancer_Type_Specific.y[!is.na(master.sheet$Cancer_Type_Specific.y) & !master.sheet$multiple_samples]
master.sheet$Cancer_Type_Details.x[!is.na(master.sheet$Cancer_Type_Details.y)] <- master.sheet$Cancer_Type_Details.y[!is.na(master.sheet$Cancer_Type_Details.y)]
colnames(master.sheet)[colnames(master.sheet) %in% c("BLOCK_ACCESSION_NBR.x", "Cancer_Type_Specific.x", "Cancer_Type_Details.x", "Primary.x")] <- 
    c("BLOCK_ACCESSION_NBR", "Cancer_Type_Specific", "Cancer_Type_Details", "Primary")
master.sheet <- master.sheet[, -c(which(colnames(master.sheet) %in% c("Cancer_Type_Specific.y", "Primary.y", "Cancer_Type_Details.y", "BLOCK_ACCESSION_NBR.y")))]
write.table(master.sheet, "master.sheet.tsv", sep = "\t", row.names = F)
## save point
master.sheet.pit <- master.sheet
master.sheet <- master.sheet.pit

## same thing for meningioma data
men.data <- read.delim("../OncDRS data/archive/Meningioma_R_Upload.txt", stringsAsFactors = F)
men.data <- men.data[, c("surg_path", "mrn", "primary", "tumor_subtype", "grade")]
colnames(men.data) <- c("BLOCK_ACCESSION_NBR", "MRN", "Primary", "Cancer_Type_Specific", "Grade")
men.data$Primary <- men.data$Primary == TRUE
men.data$BLOCK_ACCESSION_NBR <- sapply(men.data$BLOCK_ACCESSION_NBR, function(x){gsub("(BS)([0-9])", "\\1-\\2", x)}, USE.NAMES = FALSE)
men.data$BLOCK_ACCESSION_NBR <- sapply(men.data$BLOCK_ACCESSION_NBR, function(x){gsub("(BS-[0-9][0-9])([A-z0-9])", "\\1-\\2", x)}, USE.NAMES = FALSE)
master.sheet <- merge(master.sheet, men.data, "MRN", all.x = TRUE)

## checks to see which ones disagree, remove data from those
master.sheet[which(master.sheet$BLOCK_ACCESSION_NBR.y != master.sheet$BLOCK_ACCESSION_NBR.x), c("exclude")] <- "check"

master.sheet$Primary.x[!is.na(master.sheet$Primary.y) & !master.sheet$multiple_samples] <- master.sheet$Primary.y[!is.na(master.sheet$Primary.y) & !master.sheet$multiple_samples]
master.sheet$Cancer_Type_Specific.x[!is.na(master.sheet$Cancer_Type_Specific.y) & !master.sheet$multiple_samples] <- 
    master.sheet$Cancer_Type_Specific.y[!is.na(master.sheet$Cancer_Type_Specific.y) & !master.sheet$multiple_samples]
master.sheet$Grade.x[!is.na(master.sheet$Grade.y) & !master.sheet$multiple_samples] <- master.sheet$Grade.y[!is.na(master.sheet$Grade.y) & !master.sheet$multiple_samples]
colnames(master.sheet)[colnames(master.sheet) %in% c("BLOCK_ACCESSION_NBR.x", "Cancer_Type_Specific.x", "Primary.x", "Grade.x")] <- 
    c("BLOCK_ACCESSION_NBR", "Cancer_Type_Specific", "Grade", "Primary")
master.sheet <- master.sheet[, -c(which(colnames(master.sheet) %in% c("Cancer_Type_Specific.y", "Primary.y", "BLOCK_ACCESSION_NBR.y", "Grade.y")))]


## same thing for glioma data
lgg.data <- read.delim("../OncDRS data/LGG_R_Upload.txt", stringsAsFactors = F)
lgg.data <- lgg.data[, c("Surg.Path", "BWH.MRN", "Diagnosis..........Path.Report.", "Primary.vs.Recurrent")]
lgg.data$Grade <- as.integer(sapply(lgg.data$Diagnosis..........Path.Report., function(x){gsub("[A-z]", "", x)}, USE.NAMES = FALSE))
lgg.data$Cancer_Type_Specific <- sapply(lgg.data$Diagnosis..........Path.Report., function(x){gsub("[0-9]", "", x)}, USE.NAMES = FALSE)
lgg.data$Cancer_Type_Specific[lgg.data$Cancer_Type_Specific == "O"] <- "Oligo"
lgg.data$Cancer_Type_Specific[lgg.data$Cancer_Type_Specific == "A"] <- "Astro"
lgg.data$Cancer_Type_Specific[lgg.data$Cancer_Type_Specific == "DA"] <- "DiffuseAstro"
lgg.data$Cancer_Type_Specific[lgg.data$Cancer_Type_Specific == "AA"] <- "AnaplasticAstro"
lgg.data$Primary <- lgg.data$Primary.vs.Recurrent == "Primary"
lgg.data <- lgg.data[, -c(3,4)]
colnames(lgg.data) <- c("BLOCK_ACCESSION_NBR", "MRN", "Grade", "Cancer_Type_Specific", "Primary")
master.sheet <- merge(master.sheet, lgg.data, "MRN", all.x = TRUE)

## checks to see which ones disagree, remove data from those
master.sheet[which(master.sheet$BLOCK_ACCESSION_NBR.y != master.sheet$BLOCK_ACCESSION_NBR.x), c("exclude")] <- "check"


master.sheet$Primary.x[!is.na(master.sheet$Primary.y) & !master.sheet$multiple_samples & master.sheet$Primary.x == ""] <- 
    master.sheet$Primary.y[!is.na(master.sheet$Primary.y) & !master.sheet$multiple_samples & master.sheet$Primary.x == ""]
master.sheet$Cancer_Type_Specific.x[!is.na(master.sheet$Cancer_Type_Specific.y) & !master.sheet$multiple_samples & master.sheet$Cancer_Type_Specific.x == ""] <- 
    master.sheet$Cancer_Type_Specific.y[!is.na(master.sheet$Cancer_Type_Specific.y) & !master.sheet$multiple_samples & master.sheet$Cancer_Type_Specific.x == ""]
master.sheet$Grade.x[!is.na(master.sheet$Grade.y) & !master.sheet$multiple_samples & master.sheet$Grade.x == ""] <- 
    master.sheet$Grade.y[!is.na(master.sheet$Grade.y) & !master.sheet$multiple_samples & master.sheet$Grade.x == ""]
colnames(master.sheet)[colnames(master.sheet) %in% c("BLOCK_ACCESSION_NBR.x", "Cancer_Type_Specific.x", "Primary.x", "Grade.x")] <- 
    c("BLOCK_ACCESSION_NBR", "Cancer_Type_Specific", "Grade", "Primary")
master.sheet <- master.sheet[, -c(which(colnames(master.sheet) %in% c("Cancer_Type_Specific.y", "Primary.y", "BLOCK_ACCESSION_NBR.y", "Grade.y")))]

## add in dob, gender, and dod if available
master.sheet <- merge(master.sheet, identifiers[, c("PATIENT_ID", "DERIVED_DEATH_IND","GENDER_NM", "BIRTH_DT")], all.x =)


## clean up prior to writing to table for manual review
master.sheet$Grade[is.na(master.sheet$Grade)] <- ""
master.sheet$Cancer_Type_Specific[is.na(master.sheet$Cancer_Type_Specific)] <- ""
master.sheet$Primary[is.na(master.sheet$Primary)] <- ""
master.sheet$Location_CNS[is.na(master.sheet$Location_CNS)] <- ""
master.sheet$Location_detailed[is.na(master.sheet$Location_detailed)] <- ""
write.table(master.sheet, "../OncDRS data/master.sheet.tsv", row.names = F, sep = "\t")


## read in raw mutations, copy number, and rearrangement data
## start with mutations
all.mutations <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_MUTATION_RESULTS.csv", stringsAsFactors = F)

## rename to consistent gene name
all.mutations$BEST_EFF_GENE[all.mutations$BEST_EFF_GENE == "MLL"] <- "KMT2A"
all.mutations$BEST_EFF_GENE[all.mutations$BEST_EFF_GENE == "MLL2"] <- "KMT2D"
all.mutations$BEST_EFF_GENE[all.mutations$BEST_EFF_GENE == "TERT" & all.mutations$CANONICAL_VARIANT_CLASS == "TERT_Promoter"] <- "TERT_Promoter"

## standardizes mutation classification
all.mutations$variant_classification <- "other"
all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% c("Missense", "Missense_Mutation", "protein_altering", "coding_sequence")] <- "missense"
all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% 
                                         c("In_Frame_Ins",  "Inframe_Del",  "Inframe_Ins", "In_Frame_Del")] <- "in_frame_indel"

all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% c("Frame_Shift_Del","Frame_Shift_Ins","Frameshift")] <- "frameshift_indel"
all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% c("Splice_Acceptor", "Splice_Donor", "Splice_Region", "Splice_Site")] <- "splice_site"
all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% c("Stop_Lost", "Nonstop_Mutation", "incomplete_terminal_codon")]<- "stop_codon"
all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% c("Nonsense", "Nonsense_Mutation")] <- "nonsense"
#all.mutations$variant_classification[all.mutations$BEST_EFF_VARIANT_CLASS %in% c("Translation_Start_Site", "Initiator_Codon")] <- "TSS"


all.mutations.tier1.3 <- all.mutations[all.mutations$TIER_ID < 4, ]
all.mutations.tier1.4 <- all.mutations[all.mutations$TIER_ID < 5, ]
all.cnvs <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_CNV_RESULTS.csv", stringsAsFactors = F)


## read in coverage information
gene.list <- read.delim("../OncDRS data/archive/Gene_lists.txt", stringsAsFactors = F)
not.covered <- list(gene.list$Gene[-1][!as.numeric(gene.list$OncoPanel.v1[-1])], gene.list$Gene[-1][!as.numeric(gene.list$OncoPanel.v2[-1])],
                    gene.list$Gene[-1][!as.numeric(gene.list$OncoPanel.v3[-1])])
not.covered.map <- gene.list$Gene[-1][!as.numeric(gene.list$OncoMap[-1])]



## preprocessing for rearrangements
all.svs <- read.csv("../OncDRS data/REQ_ID08_65337_ONCOPANEL_SV_RESULTS.csv", stringsAsFactors = F)

all.svs$empty <- FALSE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[3], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[18], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[44], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[71], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[112], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[126], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[121], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[133], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[206], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[238], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[277], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[360], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[375], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[401], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[402], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[430], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[447], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[456], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[595], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[606], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[638], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[697], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[767], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[820], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[916], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[956], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1020], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1078], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1097], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1119], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1153], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1202], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1205], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1248], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1249], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1250], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1257], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1278], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1296], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1298], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1306], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1314], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1316], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1317], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1321], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1332], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1335], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1336], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1344], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1377], ]$empty <- TRUE
all.svs[all.svs$STRUCTURAL_VARIANT_TEXT == all.svs$STRUCTURAL_VARIANT_TEXT[1388], ]$empty <- TRUE

## write to excel sheet for manual gene name integration

write.csv(all.svs[, c("PATIENT_ID", "SAMPLE_ACCESSION_NBR", "BLOCK_ACCESSION_NBR", "STRUCTURAL_VARIANT_TEXT", "empty")], 
          "../OncDRS data/rearrangements_for_manual_review.csv")

all.svs <- read.delim("C:/Users/Noah/Syncplicity Folders/Pan-CNS Oncopanel/OncDRS data/rearrangements_for_manual_review_R_upload_20170526.txt", stringsAsFactors = F)
all.svs <- all.svs[all.svs$variant != "", ]
all.sv.indels <- all.svs[all.svs$variant == "indel", ]
all.svs <- all.svs[all.svs$variant == "rearrangement", ]

all.sv.indels <- all.sv.indels[, -7]

## create merged list of indels and SVs for combining with mutations
all.svs.formatted <- all.sv.indels
colnames(all.svs.formatted)[6:7] <- c("BEST_EFF_GENE", "variant_classificaiton")
all.svs.breakpoints <- all.svs
all.svs.breakpoints<- all.svs.breakpoints[, -6]
colnames(all.svs.breakpoints)[6] <- "Gene1"
all.svs.breakpoints <- rbind(all.svs.breakpoints, all.svs[, -7])
colnames(all.svs.breakpoints)[6:7] <- c("BEST_EFF_GENE", "variant_classification")
all.svs.breakpoints <- all.svs.breakpoints[!(all.svs.breakpoints$BEST_EFF_GENE %in% c("intergenic", "not_given", "not_giveN")), ]
colnames(all.svs.formatted) <- colnames(all.svs.breakpoints)
all.svs.formatted <- rbind(all.svs.formatted, all.svs.breakpoints)




