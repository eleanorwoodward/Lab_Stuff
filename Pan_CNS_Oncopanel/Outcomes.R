## outcomes analysis
## adapted from https://www.openintro.org/download.php?file=survival_analysis_in_R&amp%3Breferrer

install.packages("OIsurv")
library(OIsurv)

can.reg <- read.csv("../OncDRS data/REQ_ID08_65337_CANCER_DIAGNOSIS_CAREG.csv", stringsAsFactors = F)
pat.info <- read.csv("../OncDRS data/REQ_ID08_65337_PT_INFO_STATUS_REGISTRATION.csv", stringsAsFactors = F)


temp <- sapply(can.reg$DIAGNOSIS_DT, as.Date, "%d-%b-%y", USE.NAMES = F)

## can't use SAPPLY, simplifies format
DATE <- rep(as.Date(temp[1], origin = "1970-01-01"), length(temp))
for (i in 1:length(temp)){
    DATE[i] <- as.Date(temp[i], origin = "1970-01-01" )
}
can.reg$DIAGNOSIS_DT <- DATE


temp <- sapply(pat.info$DERIVED_DEATH_DT, as.Date, "%d-%b-%y", USE.NAMES = F)

## can't use SAPPLY, simplifies format
DATE <- rep(as.Date(temp[1], origin = "1970-01-01"), length(temp))
for (i in 1:length(temp)){
    DATE[i] <- as.Date(temp[i], origin = "1970-01-01" )
}
pat.info$DERIVED_DEATH_DT <- DATE

outcomes <- merge(can.reg[, c("PATIENT_ID", "DIAGNOSIS_DT")], pat.info[, c("PATIENT_ID", "DERIVED_DEATH_DT")])
outcomes$observed <- !is.na(outcomes$DERIVED_DEATH_DT)
outcomes$DERIVED_DEATH_DT[is.na(outcomes$DERIVED_DEATH_DT)] <- as.Date("2017-03-21")
outcomes$days <- as.numeric(outcomes$DERIVED_DEATH_DT - outcomes$DIAGNOSIS_DT)


outcomes.glioma <- subset(outcomes, PATIENT_ID %in% master.sheet[master.sheet$SAMPLE_ACCESSION_NBR %in% pathologies[["glioma"]], ]$PATIENT_ID)
outcomes.glioma$IDH1 <- outcomes.glioma$PATIENT_ID %in% all.mutations.tier1.3[all.mutations.tier1.3$BEST_EFF_GENE == "IDH1", ]$PATIENT_ID

my.surv.object <- Surv(as.numeric(outcomes.glioma$days), outcomes.glioma$observed)

my.fit <- survfit(my.surv.object ~ 1)


my.fit1 <- survfit(Surv(outcomes.glioma$days, outcomes.glioma$observed) ~ outcomes.glioma$IDH1)
plot(my.fit1, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("IDH1 wt", "IDH1 mutant"), lty = 1:3)
