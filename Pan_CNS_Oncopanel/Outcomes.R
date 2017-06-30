## outcomes analysis
## adapted from https://www.openintro.org/download.php?file=survival_analysis_in_R&amp%3Breferrer
## https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Cox-Regression.pdf
## https://courses.nus.edu.sg/course/stacar/internet/st3242/handouts/notes1.pdf


install.packages("OIsurv")
library(OIsurv)

# ## can't use SAPPLY, simplifies format
# DATE <- rep(as.Date(temp[1], origin = "1970-01-01"), length(temp))
# for (i in 1:length(temp)){
#     DATE[i] <- as.Date(temp[i], origin = "1970-01-01" )
# }
# can.reg$DIAGNOSIS_DT <- DATE
# 
# 
# temp <- sapply(pat.info$DERIVED_DEATH_DT, as.Date, "%d-%b-%y", USE.NAMES = F)
# 
# ## can't use SAPPLY, simplifies format
# DATE <- rep(as.Date(temp[1], origin = "1970-01-01"), length(temp))
# for (i in 1:length(temp)){
#     DATE[i] <- as.Date(temp[i], origin = "1970-01-01" )
# }
# pat.info$DERIVED_DEATH_DT <- DATE

outcome <- master.sheet[master.sheet$Cancer_Type_Broad == "Glioma" & master.sheet$exclude %in% c("", "bch") & master.sheet$Deceased %in% c("0", "1"), ]

## censored values: double check once we have all dates that this is correct
outcome$Outcomes_Date[outcome$Outcomes_Date == ""] <- "05/30/2017"
outcome$Outcomes_Date <- as.Date(outcome$Outcomes_Date, "%m/%d/%Y")

## format date, get age
outcome$DoS <- outcome$Date_of_surgery
outcome$DoS <- as.Date(outcome$DoS, "%m/%d/%Y")
outcome$time <- as.numeric(outcome$Outcomes_Date - outcome$DoS)
outcome <- outcome[!is.na(outcome$time), ]
outcome <- outcome[outcome$time > 1, ]
outcome$Date_of_birth <- outcome$DOB
outcome$DOB <- gsub("([0-9]*-[A-z]*)-([0-9])", "\\1-19\\2", outcome$DOB)
outcome$DOB <- as.Date(outcome$DOB, "%d-%b-%Y")
outcome$Age <- round(as.numeric(outcome$DoS - outcome$DOB) / 365)
outcome$Age_interval <- "Young"
outcome$Age_interval[outcome$Age >= 40 & outcome$Age < 65] <- "Middle"
outcome$Age_interval[outcome$Age >= 65] <- "Old"
outcome$Age_elderly <- outcome$Age_interval
outcome$Age_elderly[outcome$Age_elderly == "Middle"] <- "Young"

outcomes.glioma <- outcome[outcome$Cancer_Type_Broad == "Glioma", ]
master.sheet <- merge(master.sheet, outcomes.glioma[, c("SAMPLE_ACCESSION_NBR", "Age_interval")], all.x = TRUE)


## integrate mutation data from comut plot generated file
colnames(df.wide.glioma)[1] <- "SAMPLE_ACCESSION_NBR"
outcomes.glioma <- merge(outcomes.glioma, df.wide.glioma[, !(colnames(df.wide.glioma) == "Age_interval")], "SAMPLE_ACCESSION_NBR")
outcomes.glioma$GBM <- outcomes.glioma$Cancer_Type_Specific == "Glioblastoma"
outcomes.glioma$MGMT_Binary <- outcomes.glioma$MGMT
outcomes.glioma$MGMT_Binary[outcomes.glioma$MGMT_Binary == "not given" | outcomes.glioma$MGMT_Binary == ""] <- NA
outcomes.glioma$'1p19q' <- outcomes.glioma$'1p' == "1" & outcomes.glioma$'19q' == "1"

## generate data.frame for subsequent analysis
outcomes.GBM <- outcomes.glioma[outcomes.glioma$Cancer_Type_Specific == "Glioblastoma", ]
outcomes.LGG <- outcomes.glioma[outcomes.glioma$Cancer_Type_Specific != "Glioblastoma", ]

## survival objects for subsequent analysis
surv.gliomas <- Surv(outcomes.glioma$time, as.numeric(outcomes.glioma$Deceased))
surv.LGG <- Surv(outcomes.LGG$time, as.numeric(outcomes.LGG$Deceased))
surv.GBM <- Surv(outcomes.GBM$time, as.numeric(outcomes.GBM$Deceased))

## first analyze difference between LGG and GBM
fit.lgg_vs_gbm <- survfit(surv.gliomas ~ outcomes.glioma$GBM)

## assess statistical significance in difference btween curves
survdiff(surv.gliomas ~ outcomes.glioma$GBM, subset = T, na.action = options()$na.action)

plot(fit.lgg_vs_gbm, xlab = "Days", ylab = "Overall Surival", lty = 1:3, mark.time = TRUE)
legend(100, .9, c("LGG", "GBM"), lty = 1:3)

summary(fit.lgg_vs_gbm)


## baseline prediction for glioma recurrence based on histology
cox_fit.lgg_vs_gbm <- coxph(surv.gliomas ~ factor(outcomes.glioma$GBM) + factor(outcomes.glioma$Age_elderly))
summary(cox_fit.lgg_vs_gbm)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.lgg_vs_gbm)
plot(survfit(cox_fit.lgg_vs_gbm))


## compare to adding IDH1 status: use log likelihood ratio test for nested models only
cox_fit.lgg_vs_gbm_idh1 <- coxph(surv.gliomas ~ factor(outcomes.glioma$GBM) + factor(outcomes.glioma$IDH1) + factor(outcomes.glioma$Age_elderly))
summary(cox_fit.lgg_vs_gbm_idh1)
anova(cox_fit.lgg_vs_gbm, cox_fit.lgg_vs_gbm_idh1)

cox_fit.IDH1 <- coxph(surv.gliomas ~ factor(outcomes.glioma$IDH1) + factor(outcomes.glioma$CDKN2A) + factor(outcomes.glioma$`10q`) + factor(outcomes.glioma$Age_elderly))
summary(cox_fit.IDH1)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.IDH1)
plot(survfit(cox_fit.IDH1))



## subtype specific analysis



## IDH1 status
fit.GBM_IDH1 <- survfit(surv.GBM ~ outcomes.GBM$IDH1)

## assess statistical significance in difference btween curves
survdiff(surv.GBM ~ outcomes.GBM$IDH1, subset = T, na.action = options()$na.action)
pdf("Kaplan GBM IDH1.pdf")
plot(fit.GBM_IDH1, xlab = "Days", ylab = "Overall Surival", lty = 1:3, mark.time = TRUE)
legend(100, .9, c("IDH1 mut", "IDH1 wt"), lty = 1:3)
dev.off()

## ATRX within IDH1 mutant tumors
surv.GBM_IDH1_mt <- Surv(outcomes.GBM$time[outcomes.GBM$IDH1 != "z"], 
                         as.numeric(outcomes.GBM$Deceased[outcomes.GBM$IDH1 != "z"]))
fit.GBM_IHD1_mt_ATRX <- survfit(surv.GBM_IDH1_mt ~ outcomes.GBM$ATRX[outcomes.GBM$IDH1 != "z"])
plot(fit.GBM_IHD1_mt_ATRX, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("ATRX mut", "ATRX wt"), lty = 1:3)
survdiff(surv.GBM_IDH1_mt ~ outcomes.GBM$ATRX[outcomes.GBM$IDH1 != "z"], subset = T, na.action = options()$na.action)



## MGMT status
fit.GBM_MGMT <- survfit(surv.GBM[!is.na(outcomes.GBM$MGMT_Binary)] ~ outcomes.GBM$MGMT[!is.na(outcomes.GBM$MGMT_Binary)])

## assess statistical significance in difference btween curves
survdiff(surv.GBM[!is.na(outcomes.GBM$MGMT_Binary)] ~ outcomes.GBM$MGMT[!is.na(outcomes.GBM$MGMT_Binary)], subset = T, na.action = options()$na.action)
pdf("Kaplan GBM MGMT.pdf")
plot(fit.GBM_MGMT, xlab = "Days", ylab = "Overall Surival", lty = 1:3, mark.time = TRUE)
legend(100, .9, c("MGMT methlyated", "MGMT unmethlyated"), lty = 1:2)
dev.off()


## IDH1 status LGG
fit.LGG_IDH1 <- survfit(surv.LGG ~ outcomes.LGG$IDH1)

## assess statistical significance in difference btween curves
survdiff(surv.LGG ~ outcomes.LGG$IDH1, subset = T, na.action = options()$na.action)
pdf("Kaplan LGG IDH1.pdf")
plot(fit.LGG_IDH1, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("IDH1 mut", "IDH1 wt"), lty = 1:3)
dev.off()

## MGMT status LGG
fit.LGG_MGMT <- survfit(surv.LGG[!is.na(outcomes.LGG$MGMT_Binary)] ~ outcomes.LGG$MGMT_Binary[!is.na(outcomes.LGG$MGMT_Binary)])

## assess statistical significance in difference btween curves
survdiff(surv.LGG[!is.na(outcomes.LGG$MGMT_Binary)] ~ outcomes.LGG$MGMT_Binary[!is.na(outcomes.LGG$MGMT_Binary)], subset = T, na.action = options()$na.action)

plot(fit.LGG_MGMT, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("MGMT methlyated", "MGMT unmethlyated"), lty = 1:3)


## cox proportional hazards from http://www.tbrieder.org/epidata/course_e_ex04_task.pdf
## https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Cox-Regression.pdf


## baseline prediction for GBM recurrence based on age, MGMT and IDH1 status
cox_fit.GBM_baseline <- coxph(surv.GBM ~ factor(outcomes.GBM$IDH1) + factor(outcomes.GBM$Age_elderly) + factor(outcomes.GBM$MGMT_Binary))
summary(cox_fit.GBM_baseline)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.GBM_baseline)
plot(survfit(cox_fit.GBM_baseline))


## example code to get univariate predictors from multiple covariates: taken from https://www.r-bloggers.com/cox-proportional-hazards-model/
covariates <- c("MGMT_Binary", "Primary","TP53", "PTEN", "NF1", "EGFR", "ATRX", "RB1")
colnames(outcomes.GBM) <- gsub("-", "_", colnames(outcomes.GBM))
covariates <- colnames(outcomes.GBM)[c(17, 39:64, 74)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, as.numeric(Deceased))~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = outcomes.GBM)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           #HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           #HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           #HR <- paste0(HR, " (", 
                        #                HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                       })



res <- t(as.data.frame(univ_results, check.names = FALSE))
results <- as.data.frame(res)
results$q.value <- p.adjust(results$p.value, "fdr")

results[order(results$q.value), ]


## Kaplan plots for univariate predictive markers

fit.GBM_Age <- survfit(surv.GBM ~ outcomes.GBM$Age_elderly)

## assess statistical significance in difference btween curves
survdiff(surv.GBM ~ outcomes.GBM$Age_elderly, subset = T, na.action = options()$na.action)
pdf("Kaplan GBM Age.pdf")
plot(fit.GBM_Age, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("> 65", " < 65"), lty = 1:3)
dev.off()

fit.GBM_CDKN2A <- survfit(surv.GBM ~ outcomes.GBM$CDKN2A)

## assess statistical significance in difference btween curves
survdiff(surv.GBM ~ outcomes.GBM$CDKN2A, subset = T, na.action = options()$na.action)
pdf("Kaplan GBM CDKN2A.pdf")
plot(fit.GBM_CDKN2A, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("wt", "mut"), lty = 1:3)
dev.off()

fit.GBM_PDGFRA <- survfit(surv.GBM ~ outcomes.GBM$PDGFRA)

## assess statistical significance in difference btween curves
survdiff(surv.GBM ~ outcomes.GBM$PDGFRA, subset = T, na.action = options()$na.action)
pdf("Kaplan GBM CDKN2A.pdf")
plot(fit.GBM_PDGFRA, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("wt", "mut"), lty = 1:3)
dev.off()


cox_fit.GBM_new <- coxph(surv.GBM ~ factor(Age_elderly) + factor(MGMT_Binary) + factor(IDH1) +  
                             factor(CDKN2A), data = outcomes.GBM)


summary(cox_fit.GBM_new)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.GBM_new)
plot(survfit(cox_fit.GBM_new))


cox_fit.GBM_reduced <- coxph(surv.GBM ~ factor(Age_interval) + factor(MGMT_Binary) + factor(NF1) + factor(PRKDC), 
                         data = outcomes.GBM)
#factor(ATRX) +

summary(cox_fit.GBM_reduced)




# do same thing for LGG
colnames(outcomes.LGG) <- gsub("-", "_", colnames(outcomes.LGG))
covariates <- colnames(outcomes.LGG)[c(17, 39:64, 74)]
for (i in 1:length(covariates)){
    if (sum(outcomes.LGG[, colnames(outcomes.LGG) == covariates[i]] == "0") == nrow(outcomes.LGG)){
        covariates[i] <- "remove"
    }
}
covariates <- covariates[covariates != "remove"]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, as.numeric(Deceased))~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = outcomes.LGG)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           #HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           #HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           #HR <- paste0(HR, " (", 
                           #                HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                       })



res <- t(as.data.frame(univ_results, check.names = FALSE))
results <- as.data.frame(res)

results$q.value <- p.adjust(results$p.value, "fdr")
results[order(results$q.value), ]

outcomes.LGG$10q[outcomes.LGG$10q == "arm-level gain"] <- "z"
outcomes.LGG$7p[outcomes.LGG$7p == "arm-level loss"] <- "z"
outcomes.LGG$17p[outcomes.LGG$7p == "arm-level loss"] <- "z"
outcomes.LGG$7p[outcomes.LGG$7p == "arm-level loss"] <- "z"
table(outcomes.LGG$PIK3C2B_cnv)

cox_fit.LGG_new <- coxph(surv.LGG ~ factor(EGFR) + factor(PTEN) + factor(CDKN2A) +  factor(KDR) +  factor(Age_elderly) + factor(RB1) + factor(APC) +
                             factor(IDH1), data = outcomes.LGG) 
                            #+ factor(EGFR_cnv) + factor(APC) + factor(RB1) + factor(CDKN2B_cnv) +
                             #factor(SETD2) + factor(BRCA2) + factor(10q) + factor(Grade) + factor(Age_interval) + factor(NF1) +
                             #factor(TERT) + factor(7p) + factor(7q) + factor(CREBBP) + factor(ARID1B) +
                             #factor(13q) + factor(MDM4_cnv) + factor(ATM), data = outcomes.LGG)


summary(cox_fit.LGG_new)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.LGG_new)
plot(survfit(cox_fit.LGG_new))


cox_fit.LGG_reduced <- coxph(surv.LGG ~ factor(MSH6) + factor(ARID2) + factor(IDH1), data = outcomes.LGG)
#factor(ATRX) +

summary(cox_fit.LGG_reduced)


## 1P19Q status
fit.LGG_1p19q<- survfit(surv.LGG ~ outcomes.LGG$mutations1p19q)

## assess statistical significance in difference btween curves
survdiff(surv.LGG ~ outcomes.LGG$mutations1p19q, subset = T, na.action = options()$na.action)

plot(fit.LGG_1p19q, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("1p19q wt", "1p19q mut"), lty = 1:3)

## 1P19Q status
fit.LGG_age<- survfit(surv.LGG ~ outcomes.LGG$Age_interval)

## assess statistical significance in difference btween curves
survdiff(surv.LGG ~ outcomes.LGG$mutations1p19q, subset = T, na.action = options()$na.action)

plot(fit.LGG_age, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("middle", "old", "young"), lty = 1:3)


## 1P19Q status
fit.glioma_age<- survfit(surv.gliomas ~ outcomes.glioma$Age_elderly)

## assess statistical significance in difference btween curves
survdiff(surv.gliomas ~ outcomes.glioma$Age_elderly, subset = T, na.action = options()$na.action)

plot(fit.glioma_age, xlab = "Days", ylab = "Overall Surival", lty = 1:2)
legend(100, .9, c("Elderly", "Young"), lty = 1:3)




## Decision Tree outcomes analysis
outcomes.GBM$branch3_vs_4[outcomes.GBM$IDH1 == 0 & outcomes.GBM$BRAF == 1] <- "branch3"
outcomes.GBM$branch3_vs_4[outcomes.GBM$IDH1 == 0 & outcomes.GBM$BRAF == 0] <- "branch4"

fit.GBM_decision_tree <- survfit(surv.GBM[!is.na(outcomes.GBM$branch3_vs_4), ] ~ outcomes.GBM[!is.na(outcomes.GBM$branch3_vs_4), "branch3_vs_4"])
survdiff(surv.GBM[!is.na(outcomes.GBM$branch3_vs_4), ] ~ outcomes.GBM[!is.na(outcomes.GBM$branch3_vs_4), "branch3_vs_4"])
plot(fit.GBM_decision_tree, xlab = "Days", ylab = "Overall Surival", lty = 1:2)
legend(100, .9, c("BRAF +", "BRAF -"), lty = 1:3)



## NMF cluster analysis
outcomes.GBM$NMF_cluster[outcomes.GBM$SAMPLE_ACCESSION_NBR %in% rownames(classifier[classifier$subtype == "Glioblastoma" & classifier$basis == 3, ])] <- "typical"
outcomes.GBM$NMF_cluster[outcomes.GBM$SAMPLE_ACCESSION_NBR %in% rownames(classifier[classifier$subtype == "Glioblastoma" & classifier$basis != 3, ])] <- "atypical"

fit.GBM_NMF <- survfit(surv.GBM ~ outcomes.GBM$NMF_cluster)
survdiff(surv.GBM ~ outcomes.GBM$NMF_cluster)
plot(fit.GBM_NMF, xlab = "Days", ylab = "Overall Surival", lty = 1:2)
legend(100, .9, c("Non-GBM GBMs", "GBM GBMs"), lty = 1:3)

