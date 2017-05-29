## outcomes analysis
## adapted from https://www.openintro.org/download.php?file=survival_analysis_in_R&amp%3Breferrer

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

outcome <- master.sheet[master.sheet$Cancer_Type_Broad == "Glioma" & master.sheet$exclude == "" & master.sheet$Deceased %in% c("0", "1"), ]

## censored values: double check once we have all dates that this is correct
outcome$DoD[is.na(outcome$DoD)] <- "05/01/2017"
outcome$DoD <- as.Date(outcome$DoD, "%m/%d/%Y")

outcome$DoS <- outcome$Date_of_surgery
outcome$DoS <- as.Date(outcome$DoS, "%m/%d/%Y")
outcome$time <- as.numeric(outcome$DoD - outcome$DoS)
outcome <- outcome[!is.na(outcome$time), ]
outcome <- outcome[outcome$time > 1, ]
outcome$Date_of_birth <- outcome$DOB
outcome$DOB <- gsub("([0-9]*-[A-z]*)-([0-9])", "\\1-19\\2", outcome$DOB)
outcome$DOB <- as.Date(outcome$DOB, "%d-%b-%Y")
outcome$Age <- round(as.numeric(outcome$DoS - outcome$DOB) / 365)
outcome$Age_interval <- "Young"
outcome$Age_interval[outcome$Age >= 40 & outcome$Age < 65] <- "Middle"
outcome$Age_interval[outcome$Age >= 65] <- "Old"
outcomes.glioma <- outcome[outcome$Cancer_Type_Broad == "Glioma", ]
master.sheet <- merge(master.sheet, outcomes.glioma[, c("SAMPLE_ACCESSION_NBR", "Age_interval")], all.x = TRUE)


## integrate mutation data from comut plot generated file
colnames(glioma.mutations)[1] <- "SAMPLE_ACCESSION_NBR"
outcomes.glioma <- merge(outcomes.glioma, glioma.mutations, "SAMPLE_ACCESSION_NBR")
outcomes.glioma$GBM <- outcomes.glioma$Cancer_Type_Specific == "Glioblastoma"
outcomes.glioma$MGMT_Binary <- outcomes.glioma$MGMT
outcomes.glioma$MGMT_Binary[outcomes.glioma$MGMT_Binary == "not given"] <- NA
outcomes.glioma$mutations1p19q <- outcomes.glioma$mutations.1p == "arm-level loss" & outcomes.glioma$mutations.19q == "arm-level loss"

## generate data.frame for subsequent analysis
outcomes.GBM <- outcomes.glioma[outcomes.glioma$Cancer_Type_Specific == "Glioblastoma", ]
outcomes.LGG <- outcomes.glioma[outcomes.glioma$Cancer_Type_Specific != "Glioblastoma", ]

## survival objects for subsequent analysis
surv.gliomas <- Surv(outcomes.glioma$time, as.numeric(outcomes.glioma$Deceased))
surv.LGG <- Surv(outcomes.LGG$time, as.numeric(outcomes.LGG$Deceased))
surv.GBM <- Surv(outcomes.GBM$time, as.numeric(outcomes.GBM$Deceased))

## first analyze difference between LGG and GBM
fit.glioma_vs_gbm <- survfit(surv.gliomas ~ outcomes.glioma$GBM)

## assess statistical significance in difference btween curves
survdiff(surv.gliomas ~ outcomes.glioma$GBM, subset = T, na.action = options()$na.action)

plot(fit.glioma_vs_gbm, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("LGG", "GBM"), lty = 1:3)


## IDH1 status
fit.GBM_IDH1 <- survfit(surv.GBM ~ outcomes.GBM$mutations.IDH1)

## assess statistical significance in difference btween curves
survdiff(surv.GBM ~ outcomes.GBM$mutations.IDH1, subset = T, na.action = options()$na.action)

plot(fit.GBM_IDH1, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("IDH1 mut", "IDH1 wt"), lty = 1:3)


## ATRX within IDH1 mutant tumors
surv.GBM_IDH1_mt <- Surv(outcomes.GBM$time[outcomes.GBM$mutations.IDH1 != "z"], 
                         as.numeric(outcomes.GBM$Deceased[outcomes.GBM$mutations.IDH1 != "z"]))
fit.GBM_IHD1_mt_ATRX <- survfit(surv.GBM_IDH1_mt ~ outcomes.GBM$mutations.ATRX[outcomes.GBM$mutations.IDH1 != "z"])
plot(fit.GBM_IHD1_mt_ATRX, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("ATRX mut", "ATRX wt"), lty = 1:3)
survdiff(surv.GBM_IDH1_mt ~ outcomes.GBM$mutations.ATRX[outcomes.GBM$mutations.IDH1 != "z"], subset = T, na.action = options()$na.action)



## MGMT status
fit.GBM_MGMT <- survfit(surv.GBM ~ outcomes.GBM$MGMT)

## assess statistical significance in difference btween curves
survdiff(surv.GBM ~ outcomes.GBM$MGMT, subset = T, na.action = options()$na.action)

plot(fit.GBM_MGMT, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("MGMT methlyated", "unknown", "MGMT unmethlyated"), lty = 1:3)



## IDH1 status LGG
fit.LGG_IDH1 <- survfit(surv.LGG ~ outcomes.LGG$mutations.IDH1)

## assess statistical significance in difference btween curves
survdiff(surv.LGG ~ outcomes.LGG$mutations.IDH1, subset = T, na.action = options()$na.action)

plot(fit.LGG_IDH1, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("IDH1 mut", "IDH1 wt"), lty = 1:3)


## MGMT status LGG
fit.LGG_MGMT <- survfit(surv.LGG ~ outcomes.LGG$MGMT_Binary)

## assess statistical significance in difference btween curves
survdiff(surv.LGG ~ outcomes.LGG$MGMT_Binary, subset = T, na.action = options()$na.action)

plot(fit.LGG_MGMT, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("MGMT methlyated", "unknown", "MGMT unmethlyated"), lty = 1:3)


## cox proportional hazards from http://www.tbrieder.org/epidata/course_e_ex04_task.pdf
## https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Cox-Regression.pdf


## baseline prediction for GBM recurrence based on age, MGMT and IDH1 status
cox_fit.GBM_baseline <- coxph(surv.GBM ~ factor(outcomes.GBM$mutations.IDH1) + factor(outcomes.GBM$Age_interval) + factor(outcomes.GBM$MGMT_Binary))
summary(cox_fit.GBM_baseline)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.GBM_baseline)
plot(survfit(cox_fit.GBM_baseline))


## example code to get univariate predictors from multiple covariates: taken from https://www.r-bloggers.com/cox-proportional-hazards-model/
covariates <- c("MGMT_Binary", "Primary","mutations.TP53", "mutations.PTEN", "mutations.NF1", "mutations.EGFR", "mutations.ATRX", "mutations.RB1")
colnames(outcomes.GBM) <- gsub("-", "_", colnames(outcomes.GBM))
covariates <- colnames(outcomes.GBM)[c(17, 30, 37:97, 101)]
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

results[order(results$p.value), ]

outcomes.GBM$mutations.10q[outcomes.GBM$mutations.10q == "arm-level gain"] <- "z"
outcomes.GBM$mutations.19p[outcomes.GBM$mutations.19p == "arm-level loss"] <- "z"
cox_fit.GBM_new <- coxph(surv.GBM ~ factor(Age_interval) + factor(MGMT_Binary) + factor(mutations.IDH1) + factor(mutations.10q) +  
                             factor(mutations.CDKN2A_cnv) + factor(mutations.CDKN2B_cnv) + factor(mutations.PDGFRA_cnv) + factor(mutations.NF1) +
                             factor(mutations.PTEN) + factor(mutations.19p) + factor(mutations.BRCA2) + factor(mutations.PRKDC) + factor(mutations.19q), 
                            data = outcomes.GBM)
                                #factor(mutations.ATRX) +

summary(cox_fit.GBM_new)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.GBM_new)
plot(survfit(cox_fit.GBM_new))


cox_fit.GBM_reduced <- coxph(surv.GBM ~ factor(Age_interval) + factor(MGMT_Binary) + factor(mutations.NF1) + factor(mutations.PRKDC), 
                         data = outcomes.GBM)
#factor(mutations.ATRX) +

summary(cox_fit.GBM_reduced)




# do same thing for LGG
colnames(outcomes.LGG) <- gsub("-", "_", colnames(outcomes.LGG))
covariates <- colnames(outcomes.LGG)[c(17, 30, 37:117, 121:122)]
for (i in 1:length(covariates)){
    if (sum(outcomes.LGG[, colnames(outcomes.LGG) == covariates[i]] == "z") == nrow(outcomes.LGG)){
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

results[order(results$q.value), ]
results$q.value <- p.adjust(results$p.value, "fdr")

outcomes.LGG$mutations.10q[outcomes.LGG$mutations.10q == "arm-level gain"] <- "z"
outcomes.LGG$mutations.7p[outcomes.LGG$mutations.7p == "arm-level loss"] <- "z"
outcomes.LGG$mutations.17p[outcomes.LGG$mutations.7p == "arm-level loss"] <- "z"
outcomes.LGG$mutations.7p[outcomes.LGG$mutations.7p == "arm-level loss"] <- "z"
table(outcomes.LGG$mutations.PIK3C2B_cnv)

cox_fit.LGG_new <- coxph(surv.LGG ~ factor(mutations.EGFR) + factor(mutations.MSH6) + factor(mutations.PTEN) + factor(mutations.ARID2) + factor(mutations.CDKN2A_cnv) +
                             factor(mutations.IDH1), data = outcomes.LGG) 
                            #+ factor(mutations.EGFR_cnv) + factor(mutations.APC) + factor(mutations.RB1) + factor(mutations.CDKN2B_cnv) +
                             #factor(mutations.SETD2) + factor(mutations.BRCA2) + factor(mutations.10q) + factor(mutations.Grade) + factor(Age_interval) + factor(mutations.NF1) +
                             #factor(mutations.TERT) + factor(mutations.7p) + factor(mutations.7q) + factor(mutations.CREBBP) + factor(mutations.ARID1B) +
                             #factor(mutations.13q) + factor(mutations.MDM4_cnv) + factor(mutations.ATM), data = outcomes.LGG)


summary(cox_fit.LGG_new)

## test if prorportional hazards assumption is violated: p <.05 yes for individual variable, GLOBAL for overall effect
cox.zph(cox_fit.LGG_new)
plot(survfit(cox_fit.LGG_new))


cox_fit.LGG_reduced <- coxph(surv.LGG ~ factor(mutations.MSH6) + factor(mutations.ARID2) + factor(mutations.IDH1), data = outcomes.LGG)
#factor(mutations.ATRX) +

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

