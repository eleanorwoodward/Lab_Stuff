## analysis for pediatric meningiomas paper

install.packages("OIsurv")
library(OIsurv)

outcomes.data <- read.delim("C:/Users/Noah/OneDrive/Work/Pediatric Meningioma/Cleaned_Data.txt", stringsAsFactors = F)
outcomes.data <- outcomes.data[-13, ]

## survival objects for subsequent analysis
surv.ped <- Surv(outcomes.data$Time, as.numeric(outcomes.data$Endpoint.reached))

## analyze each condition
fit.cranial <- survfit(surv.ped ~ outcomes.data$Cranial)
fit.dural <- survfit(surv.ped ~ outcomes.data$Dural.invasion)
fit.hydro <- survfit(surv.ped ~ outcomes.data$Hydrocephalus.)
fit.nf2 <- survfit(surv.ped ~ outcomes.data$NF.2.)
fit.radiation <- survfit(surv.ped ~ outcomes.data$Previous.exposure.to.radiation)
fit.gtr <- survfit(surv.ped ~ outcomes.data$Total.resection)
fit.grade <- survfit(surv.ped ~ outcomes.data$Grade.1)

## assess statistical significance in difference btween curves
survdiff(surv.ped ~ outcomes.data$Cranial, subset = T, na.action = options()$na.action)
survdiff(surv.ped ~ outcomes.data$Dural.invasion, subset = T, na.action = options()$na.action)
survdiff(surv.ped ~ outcomes.data$Hydrocephalus., subset = T, na.action = options()$na.action)
survdiff(surv.ped ~ outcomes.data$NF.2., subset = T, na.action = options()$na.action)
survdiff(surv.ped ~ outcomes.data$Previous.exposure.to.radiation, subset = T, na.action = options()$na.action)
survdiff(surv.ped ~ outcomes.data$Total.resection, subset = T, na.action = options()$na.action)
survdiff(surv.ped ~ outcomes.data$Grade.1, subset = T, na.action = options()$na.action)

p.vals <- c(0.0711, 0.658, 0.0422, 0.964, 0.0929, 0.964, 0.062)
q.vals <- p.adjust(p.vals, "fdr")



plot(fit.glioma_vs_ped, xlab = "Days", ylab = "Overall Surival", lty = 1:3)
legend(100, .9, c("ped", "ped"), lty = 1:3)
