rm(list=ls())

####Load data####
#install needed packages
library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)
library(ggsci)

dat <- read.csv(file = "dat.surv.int.csv", header = T)

#time to event = diff_days
surv_object <- Surv(time = dat$diff_days, event = dat$interaction)
surv_object

#survival difference between match / mismatch
fit1 <- survfit(surv_object ~ Camouflage.mismatch, data = dat)
fit1

summary(fit1)

windows()
ggsurvplot(fit1, data = dat, palette = "npg", conf.int = T, linetype = 1,
           pval = T, pval.method = T, surv.median.line = "hv", 
           risk.table = "percentage", ggtheme = theme_classic2(),
           xlab = "Time (days)", legend.title = "Camouflage",
           legend.labs = c("Match", "Mismatch"), break.time.by = 2,
           ylab = "Decoy survival probability",
           font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
           fontsize= 7,
           tables.theme = theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20)))

#survival difference between camouflage in different habitats
dat$cxh <- paste(dat$Camouflage.mismatch, dat$habitat)

fit2 <- survfit(surv_object ~ cxh, data = dat)
fit2
summary(fit2)

ggsurvplot(fit2, data = dat, palette = "npg", conf.int = T, linetype = 1,
           pval = T, pval.method = T, surv.median.line = "hv", 
           risk.table = "percentage", ggtheme = theme_classic2(),
           xlab = "Time (days)", legend.title = "Habitat", break.time.by = 2,
           legend.labs = c("Closed-Match", "Open-Match", "Closed-Mismatch",
                           "Open-Mismatch"),
           ylab = "Decoy survival probability",
           font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
           conf.int.alpha = 0.2, conf.int.style = "ribbon",
           fontsize= 7,
           tables.theme = theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20)))

# Fit a Cox proportional hazards model mismatch & habitat separate
fit.coxph <- coxph(surv_object ~ Camouflage.mismatch + habitat, data = dat)
summary(fit.coxph)
fit.coxph %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph

# Fit a Cox proportional hazards model
dat$cxh <- factor(dat$cxh, levels = c("Match O","Mismatch O",
                                            "Match C", "Mismatch C"))
fit.coxph.1 <- coxph(surv_object ~ cxh, data = dat)
summary(fit.coxph.1)
fit.coxph.1 %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph.1

#surival analysis without repeated measures
dat.2 <- dat  %>% 
  group_by(deploymentID) %>% 
  filter(row_number()==1)

surv_object.1 <- Surv(time = dat.2$diff_days, event = dat.2$interaction)
surv_object.1

#survival difference between match / mismatch
fit3 <- survfit(surv_object.1 ~ Camouflage.mismatch, data = dat.2)
fit3

summary(fit3)

ggsurvplot(fit3, data = dat.2, palette = "npg", conf.int = T, linetype = 1,
           pval = T, pval.method = T, surv.median.line = "hv", 
           risk.table = "percentage", ggtheme = theme_classic2(),
           xlab = "Time (days)", legend.title = "Camouflage",
           legend.labs = c("Match", "Mismatch"), break.time.by = 2,
           ylab = "Decoy survival probability",
           font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
           fontsize= 7,
           tables.theme = theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20)))

# Fit a Cox proportional hazards model mismatch & habitat separate
fit.coxph.2 <- coxph(surv_object.1 ~ Camouflage.mismatch + habitat,
                     data = dat.2)
summary(fit.coxph.2)
fit.coxph.2 %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph.2

#survival difference between camouflage in different habitats without rep int
dat.2$cxh <- paste(dat.2$Camouflage.mismatch, dat.2$habitat)

fit4 <- survfit(surv_object.1 ~ cxh, data = dat.2)
fit4
summary(fit4)

ggsurvplot(fit4, data = dat.2, palette = "npg", conf.int = T, linetype = 1,
           pval = T, pval.method = T, surv.median.line = "hv", 
           risk.table = "percentage", ggtheme = theme_classic2(),
           xlab = "Time (days)", legend.title = "Habitat", break.time.by = 2,
           legend.labs = c("Closed-Match", "Open-Match", "Closed-Mismatch",
                           "Open-Mismatch"),
           ylab = "Decoy survival probability",
           font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
           conf.int.alpha = 0.2, conf.int.style = "ribbon",
           fontsize= 7,
           tables.theme = theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20)))

# Fit a Cox proportional hazards model
dat.2$cxh <- factor(dat.2$cxh, levels = c("Match O","Mismatch O",
                                            "Match C", "Mismatch C"))
fit.coxph.3 <- coxph(surv_object.1 ~ cxh, data = dat.2)
summary(fit.coxph.3)
fit.coxph.3 %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph.3

#survival analysis with detection, not only interaction
dat.3 <- read.csv(file = "dat.surv.dec.csv", header = T)

#time = diff_days 
surv_object.2 <- Surv(time = dat.3$diff_days, event = dat.3$event)
surv_object.2

#survival difference between match / mismatch
fit5 <- survfit(surv_object.2 ~ Camouflage.mismatch, data = dat.3)
fit5
summary(fit5)

ggsurvplot(fit5, data = dat.3, palette = "npg", conf.int = T, linetype = 1,
           pval = T, pval.method = T, surv.median.line = "hv", 
           risk.table = "percentage", ggtheme = theme_classic2(),
           xlab = "Time (days)", legend.title = "Camouflage",
           legend.labs = c("Match", "Mismatch"), break.time.by = 2,
           ylab = "Decoy survival probability",
           font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
           fontsize= 7,
           tables.theme = theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20)))

fit.coxph.4 <- coxph(surv_object.2 ~ Camouflage.mismatch + habitat,
                     data = dat.3)
summary(fit.coxph.4)
fit.coxph.4 %>%
  gtsummary::tbl_regression(exp = T)

#survival difference for mismatch / match within certain habitat
dat.3$cxh <- paste(dat.3$Camouflage.mismatch, dat.3$habitat)

dat.3$cxh <- factor(dat.3$cxh, levels = c("Match O","Mismatch O",
                                                "Match C", "Mismatch C"))

fit6 <- survfit(surv_object.2 ~ cxh, data = dat.3)
fit6
summary(fit6)

ggsurvplot(fit6, data = dat.3, palette = "npg", conf.int = T, linetype = 1,
           pval = T, pval.method = T, surv.median.line = "hv", 
           risk.table = "percentage", ggtheme = theme_classic2(),
           xlab = "Time (days)", legend.title = "Habitat", break.time.by = 2,
           legend.labs = c("Closed-Match", "Open-Match", "Closed-Mismatch",
                           "Open-Mismatch"),
           ylab = "Decoy survival probability",
           font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
           conf.int.alpha = 0.2, conf.int.style = "ribbon",
           fontsize= 7,
           tables.theme = theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20)))

fit.coxph.5 <- coxph(surv_object.2 ~ cxh, data = dat.3)
summary(fit.coxph.3)
fit.coxph.5 %>%
  gtsummary::tbl_regression(exp = T)


####END####
