rm(list=ls())

####Load data####
#install needed packages
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(readxl)
library(survival)
library(survminer)
library(jtools)
library(multcomp)
library(car)
library(gtsummary)
library(ggsci)
library(rcompanion)
library(cowplot)
library(interactions)

#set working directory
setwd("C:/Users/piete/OneDrive/Documenten/Study/Master/Master research (2) camouflage mismatch/Data") 

#Csv cam trap data (change nr. 45050 to 2022-05-02)
dat <- read.csv(file="observations_v4 [17-05] (3).csv", header = T)

####Process data for analysis####
dat <- dat %>%
  mutate(predator = ifelse(commonName %in% c("Red Fox", "Domestic Cat",
                           "Domestic Dog", "European Badger",
                           "Common Raven", "Hooded Crow", "Short-eared Owl"),
                           1, 0))

dat.pred <- dat %>%
  filter(predator == 1) %>%
  group_by(predator) %>%
  summarise(n = n())

#filter for deployments with no detection
dat.x <- dat %>%
  dplyr::group_by(deploymentID) %>%
  dplyr::summarise(detect = max(predator))
sum(dat.x$detect)

dat <- left_join(dat, dat.x, by = "deploymentID")

dat.x <- dat %>%
  dplyr::filter(detect == 0) 

#calc first and last row of a deployment with no detection
first <- c(1)
for(i in 2:nrow(dat.x))
{
  add <- ifelse(dat.x$deploymentID[i]==dat.x$deploymentID[i-1],0,1)
  first <- c(first, add)
}
dat.x$first <- first
dat.x$last <- if_else(dat.x$comments == "Pickup", 1 , 0)

dat.x <- dat.x %>%
  filter(first == 1 | last == 1)# data list with deployments with no detection

#filter for predators & join with undetected deployments
dat.1 <- dat %>%
  filter(predator == 1)
dat.2 <- bind_rows(dat.1, dat.x)

#deploy_date
depl.dat <- read_excel(path = "Camera_dat.xlsx", sheet = "deploy_date")
depl.dat <- depl.dat %>%
  rename(deploymentID = deployment_id)

#join dates
dat.2 <- left_join(dat.2, depl.dat, by = "deploymentID")

####Create usefull columns####
dat.2 <- dat.2 %>%
  separate(deploymentID, c("round", "deployment"), remove = F) %>%
  separate(deployment, c("week", "treatment"), sep = 3) %>%
  separate(treatment, c("habitat", "treatment"), sep = 1) %>%
  separate(treatment, c("cluster", "treatment"), sep = 2) %>%
  separate(deploymentID, c("cluster_id"), sep = 9, remove = F)

dat.2 <- dat.2 %>%
  separate(timestamp, c("datetime", "GMT"), sep = 19) %>%
  separate(datetime, into = c("date", "time"), sep = "T", remove = FALSE)

dat.2 <- dat.2 %>% 
  mutate(interaction = ifelse(Interaction_with_lure == "True", 1, 0))

dat.2$delta_date <- as.double(difftime(lubridate::ymd(dat.2$date),
                        lubridate::ymd(dat.2$deploy_date),
                        unit = "days"))

dat.2$first[is.na(dat.2$first)] <- 0
dat.3 <- dat.2 %>%
  filter(first == 0)#only keep end dates of undetected deployments

dat.4 <- dat.3 %>%
  group_by(deploymentID, cluster_id, round, week, habitat, cluster, treatment,
           date, sequenceID, commonName, Snow.cover, Camouflage.mismatch,
           Interaction_with_lure, deploy_date, interaction, delta_date) %>%
  summarise(n = max(count))

dat.4$deploymentID <- as.factor(dat.4$deploymentID)

#calculate time diff between sequences per deployment
dat.5 <- dat.4 %>%
  arrange(deploymentID, delta_date) %>%
  group_by(deploymentID) %>%
  mutate(diff = delta_date - lag(delta_date),
         diff_days = as.numeric(diff, units = 'days')) %>%
  mutate(diff_days = coalesce(diff_days, delta_date))

####Detection difference between treatments####
dat.4 <- dat.3 %>%
  dplyr::group_by(deploymentID, cluster_id, round, week, habitat, cluster,
                  treatment, commonName,sequenceID) %>%
  dplyr::summarise(n = max(count))

dat.5$n[dat.5$commonName == ""] <- 0
dat.5$n[dat.5$commonName == "Human"] <- 0

#add cam_effort to data
depl.dat <- read.csv(file="depl_dat.csv", header = T)
depl.dat1 <- depl.dat %>%
  dplyr::select(deploymentID, cam_effort)

dat.5 <- left_join(dat.5, depl.dat1, by = "deploymentID")

dat.6 <- dat.5 %>%
  dplyr::group_by(deploymentID, cluster_id, round, week, habitat, cluster,
                  treatment, cam_effort) %>%
  dplyr::summarise(n = sum(n))
dat.6$treatment <- factor(dat.6$treatment, c("C", "W", "B"))

#Calc overall camouflage per deployment
calc_mode <- function(x){
  
  distinct_values <- unique(x)
  
  distinct_tabulate <- tabulate(match(x, distinct_values))
  
  distinct_values[which.max(distinct_tabulate)]
}

dat.y <- dat %>%
  filter(Camouflage.mismatch != "")

dat.y.1 <- dat.y %>%
  group_by(deploymentID) %>%
  mutate(Camouflage.mismatch.2 = calc_mode(Camouflage.mismatch)) %>%
  summarise(Camouflage.mismatch.2 = unique(Camouflage.mismatch.2))

miss.dat <- read_excel(path = "Camera_dat.xlsx", sheet = "empty_camouflage")
miss.dat <- miss.dat %>%
  rename(deploymentID = deployment_id)

dat.y.1 <- bind_rows(dat.y.1, miss.dat)

dat.tot <- depl.dat %>%
  select(-cam_effort) %>%
  left_join(dat.y.1, by = "deploymentID")

dat.tot$Camouflage.mismatch.2[is.na(dat.tot$Camouflage.mismatch.2)] <- "Control"

dat.6 <- left_join(dat.6, dat.tot, by = "deploymentID")

dat.6$Camouflage.mismatch.2 <- factor(dat.6$Camouflage.mismatch.2,
                          c("Control", "Match", "Mismatch"))

#visualization of raw data (not used in report)
value_max <- dat.6 %>%
  dplyr::group_by(Camouflage.mismatch.2) %>% 
  dplyr::summarize(max_value = max(n))
value_max$labels <- c("AB", "A", "B") #labels from Tukey test in line 162
sum(dat.6$n)

windows()
ggplot(data = dat.6, aes(x = Camouflage.mismatch.2, y = n)) +
  geom_boxplot(data = dat.6, aes(x = Camouflage.mismatch.2,
                                 y = n, fill = Camouflage.mismatch.2),
               alpha = 0.5, size = 1.2) + 
  geom_point(data = dat.6, aes(x= Camouflage.mismatch.2, y = n,
                               fill = Camouflage.mismatch.2), 
             position = position_jitterdodge(jitter.width=0.3,jitter.height=0),
             size = 5, shape = 21) +
  geom_text(data=value_max,aes(x=Camouflage.mismatch.2, y=15,
                               label=labels), size = 8) +
  scale_y_sqrt(name = "Predators (n)") + 
  scale_x_discrete(labels = c("Control", "Match","Mismatch"),
                   name = "Camouflage") +
  scale_fill_npg() + 
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=30))

#statistics all treatments separate(?)
dat.6 <- dat.6 %>% 
  mutate(snow.cover = ifelse(
    treatment == "B" & Camouflage.mismatch.2 == "Match", 0,
    ifelse(treatment == "B" & Camouflage.mismatch.2 == "Mismatch", 1, 
    ifelse(treatment == "W" & Camouflage.mismatch.2 == "Match", 1, 
    ifelse(treatment == "W" & Camouflage.mismatch.2 == "Mismatch", 0, NA)))))

dat.x1 <- dat.6[,c(1,2,5,7,8,15)]

m1 <- glmer(n ~ Camouflage.mismatch.2 + habitat + offset(log(cam_effort)) + (1|cluster_id), data = dat.6, family = poisson)
summary(m1)

#tukey HSD test for differences between groups
m1a <- glht(m1, linfct= mcp(Camouflage.mismatch.2 = "Tukey"))
summary(m1a) #significant difference between groups B & W but not between rest

m2a <- glht(m1, linfct= mcp(habitat = "Tukey"))
summary(m2a)

#filter C out of data
dat.6.cam <- dat.6 %>%
  filter(treatment != "C")

#test for interaction between treatment & snow.cover
dat.6.cam$snow.cover <- factor(dat.6.cam$snow.cover, c(1,0))
dat.6.cam$treatment <- factor(dat.6.cam$treatment, c("B", "W"))
m2 <- glmer(n ~ treatment + treatment:snow.cover +
              offset(log(cam_effort)) + (1|cluster_id), data = dat.6.cam, family = poisson)

summary(m2)
anova(m2)

#plot predicted values based on GLMM
p1 <- effect_plot(m1, pred = Camouflage.mismatch.2, interval = T,
                  x.label = "Camouflage", y.label = "Predators (n)",
                  data = dat.6, colors = c( '#0072b2', '#FFC107', '#d55e00'),
                  line.thickness = 3, cat.pred.point.size = 5) + 
  theme_nice() +
  theme(text = element_text(size=30))

p2 <- cat_plot(m2, data = dat.6.cam, pred = treatment, modx = snow.cover, interval = T, colors = c('#0072b2', '#FFC107'),
         x.label = "Decoy colour", y.label = "Predators (n)", pred.labels = c("Brown", "White"), modx.labels = c("Snow\npresent", "Snow\nabsent"),
         legend.main = "Snow\ncover", line.thickness = 3, pred.point.size = 5) + theme_nice() + theme(text = element_text(size=30))

windows()
plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, label_size = 20)

#calc mean difference between match & mismatch with control
dat.6.1 <- dat.6 %>%
  select(cluster_id, Camouflage.mismatch.2, n)
dat.6.1 <- dat.6.1[,c(7:9)]

dat.6.1 <- dat.6.1 %>%
  group_by(cluster_id) %>%
  mutate(delta = n - n[Camouflage.mismatch.2 == "Control"])

dat.6.1.avg <- dat.6.1 %>%
  group_by(Camouflage.mismatch.2) %>%
  summarise(avg = mean(delta))

dat.6.2.avg <- dat.6 %>%
  group_by(habitat) %>%
  summarise(avg = mean(n))

####GLMM for interaction predictors####
#use dat.5
dat.5 <- dat.5 %>%
  mutate(mammal = ifelse(commonName %in% c("Red Fox", "Domestic Cat",
                                      "Domestic Dog", "European Badger"),1, 0))

dat.7.int <- dat.5 %>%
  filter(treatment != "C" & commonName != "" & commonName != "Human")

m4 <- glmer(interaction ~ Camouflage.mismatch + mammal + (1|cluster_id), 
            data = dat.7.int, family = binomial, na.action = "na.fail")
summary(m4)


#separate mammalian and avian predators
dat.7.mam <- dat.5 %>%
  filter(mammal == 1) %>%
  filter(treatment != "C")
dat.7.mam$Camouflage.mismatch <- as.character(dat.7.mam$Camouflage.mismatch)
dat.7.mam$Camouflage.mismatch <- as.factor(dat.7.mam$Camouflage.mismatch)

dat.7.av <- dat.5 %>%
  filter(mammal == 0) %>%
  filter(treatment != "C") %>%
  filter(commonName != "" & commonName != "Human")
dat.7.av$Camouflage.mismatch <- as.character(dat.7.av$Camouflage.mismatch)
dat.7.av$Camouflage.mismatch <- as.factor(dat.7.av$Camouflage.mismatch)

m5 <- glmer(interaction ~ Camouflage.mismatch + (1|cluster_id), 
            data = dat.7.mam, family = binomial)
summary(m5)

m50 <- glmer(interaction ~ (1|cluster_id), 
            data = dat.7.mam, family = binomial)
summary(m50)
AIC(m5, m50)

m6 <- glmer(interaction ~ Camouflage.mismatch + (1|cluster_id), 
            data = dat.7.av, family = binomial)
summary(m6)
m60 <- glmer(interaction ~ (1|cluster_id), 
             data = dat.7.av, family = binomial)


windows()      
#Jtools effect plot, large overlap in CI
#mammalian predators
effect_plot(m5, pred = Camouflage.mismatch, interval = T,
            x.label = "Camouflage", y.label = "Decoy interaction probability"
            , data = dat.7.mam, cat.geom = "line", line.thickness = 3,
            cat.pred.point.size = 7, colors = c('#DDAA33', '#BB5566')) + 
  theme_nice() +
  theme(text = element_text(size=30))
  
#avian predators
effect_plot(m6, pred = Camouflage.mismatch, interval = T,
            x.label = "Camouflage", y.label = "Decoy interaction probability"
            , data = dat.7.av, cat.geom = "line", line.thickness = 3,
            cat.pred.point.size = 7, colors = c('#DDAA33', '#BB5566')) + 
  theme_nice() +
  theme(text = element_text(size=30))

#summary of interaction data (raw data)
dat.7.sum <- dat.7.int %>%
  group_by(mammal, Camouflage.mismatch) %>%
  summarise(sum = sum(interaction))

####Survival analysis####
#filter out control
dat.8 <- dat.5 %>%
  filter(treatment != "C")

#create dataset with interactions and tot deploy time for cams with no int
dat.9 <- dat.8 %>%
  filter(interaction == 1)

n_distinct(dat.9$deploymentID)

dat.test <- dat.8 %>%
  filter(interaction != 1)
n_distinct(dat.test$deploymentID)
 
dat.noint <- dat %>%
  filter(Interaction_with_lure == "False")

first <- c(1)
for(i in 2:nrow(dat.noint))
{
  add <- ifelse(dat.noint$deploymentID[i]==dat.noint$deploymentID[i-1],0,1)
  first <- c(first, add)
}
dat.noint$first <- first
dat.noint$last <- if_else(dat.noint$comments == "Pickup", 1 , 0)

dat.noint <- dat.noint %>%
  filter(first == 1 | last == 1)

#join dates
dat.noint <- left_join(dat.noint, depl.dat, by = "deploymentID")

dat.noint <- dat.noint %>%
  separate(timestamp, c("datetime", "GMT"), sep = 19) %>%
  separate(datetime, into = c("date", "time"), sep = "T", remove = FALSE)

dat.noint <- dat.noint %>%
  filter(last == 1)

dat.noint$delta_date <- as.double(difftime(lubridate::ymd(dat.noint$date),
                                           lubridate::ymd(dat.noint$deploy_date),
                                           unit = "days"))
dat.time <- dat.noint %>%
  select(deploymentID, delta_date)

dat.9 <- dat.8 %>%
  group_by(deploymentID) %>%
  filter(any(interaction==1)) 

test.2 <- dat.8 %>%
  filter(!deploymentID %in% dat.9$deploymentID)
n_distinct(test.2$deploymentID)

test.2 <- test.2 %>%
  group_by(deploymentID, habitat, treatment) %>%
  summarise(interaction = max(interaction))
test.2 <- left_join(test.2, dat.time, by = "deploymentID")
test.2 <- left_join(test.2, dat.tot[,c(2:7)], by = "deploymentID")
test.2 <- test.2 %>%
  rename(Camouflage.mismatch = Camouflage.mismatch.2) %>%
  rename(diff_days = delta_date) %>%
  select(deploymentID, habitat, treatment, interaction, diff_days,
         Camouflage.mismatch)

dat.10 <- dat.8 %>%
  filter(interaction == 1)

test.3 <- dat.10 %>%
  select(deploymentID, habitat, treatment, interaction, diff_days,
         Camouflage.mismatch)

#data set for suvival analysis interaction with repeated measures
dat.11 <- bind_rows(test.2, test.3) 
dat.11[58,5] <- 4
dat.11$Camouflage.mismatch <- as.factor(dat.11$Camouflage.mismatch)
dat.11$habitat <- as.factor(dat.11$habitat)
dat.11$treatment <- as.factor(dat.11$treatment)

#time to event = diff_days
surv_object <- Surv(time = dat.11$diff_days, event = dat.11$interaction)
surv_object

#survival difference between match / mismatch
fit1 <- survfit(surv_object ~ Camouflage.mismatch, data = dat.11)
fit1

summary(fit1)

windows()
ggsurvplot(fit1, data = dat.11, palette = "npg", conf.int = T, linetype = 1,
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
dat.11$cxh <- paste(dat.11$Camouflage.mismatch, dat.11$habitat)

fit2 <- survfit(surv_object ~ cxh, data = dat.11)
fit2
summary(fit2)

ggsurvplot(fit2, data = dat.11, palette = "npg", conf.int = T, linetype = 1,
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
fit.coxph <- coxph(surv_object ~ Camouflage.mismatch + habitat, data = dat.11)
summary(fit.coxph)
fit.coxph %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph
ggforest(fit.coxph, data = dat.11)

# Fit a Cox proportional hazards model
dat.11$cxh <- factor(dat.11$cxh, levels = c("Match O","Mismatch O",
                                            "Match C", "Mismatch C"))
fit.coxph.1 <- coxph(surv_object ~ cxh, data = dat.11)
summary(fit.coxph.1)
fit.coxph.1 %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph.1
ggforest(fit.coxph.1, data = dat.11)

#surival analysis without repeated measures
dat.12 <- dat.11  %>% 
  group_by(deploymentID) %>% 
  filter(row_number()==1)

surv_object.1 <- Surv(time = dat.12$diff_days, event = dat.12$interaction)
surv_object.1

#survival difference between match / mismatch
fit3 <- survfit(surv_object.1 ~ Camouflage.mismatch, data = dat.12)
fit3

summary(fit3)

ggsurvplot(fit3, data = dat.12, palette = "npg", conf.int = T, linetype = 1,
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
                     data = dat.12)
summary(fit.coxph.2)
fit.coxph.2 %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph.2
ggforest(fit.coxph.2, data = dat.12)

#survival difference between camouflage in different habitats without rep int
dat.12$cxh <- paste(dat.12$Camouflage.mismatch, dat.12$habitat)

fit4 <- survfit(surv_object.1 ~ cxh, data = dat.12)
fit4
summary(fit4)

ggsurvplot(fit4, data = dat.12, palette = "npg", conf.int = T, linetype = 1,
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
dat.12$cxh <- factor(dat.12$cxh, levels = c("Match O","Mismatch O",
                                            "Match C", "Mismatch C"))
fit.coxph.3 <- coxph(surv_object.1 ~ cxh, data = dat.12)
summary(fit.coxph.3)
fit.coxph.3 %>%
  gtsummary::tbl_regression(exp = T)

fit.coxph.3
ggforest(fit.coxph.3, data = dat.12)

#survival analysis with detection, not only interaction
dat.dect <- dat.8 %>%
  mutate(event = ifelse(commonName == "", 0, 1)) %>%
  select(-Camouflage.mismatch)

dat.dect <- dat.dect %>% 
  select(-deploy_date) %>%
  left_join(dat.tot, by = "deploymentID")

dat.dect <- dat.dect %>%
  rename(Camouflage.mismatch = Camouflage.mismatch.2)

dat.dect$Camouflage.mismatch <- as.factor(dat.dect$Camouflage.mismatch)
dat.dect$habitat <- as.factor(dat.dect$habitat)
dat.dect$treatment <- as.factor(dat.dect$treatment)

#time = diff_days 
surv_object.2 <- Surv(time = dat.dect$diff_days, event = dat.dect$event)
surv_object.2

#survival difference between match / mismatch
fit5 <- survfit(surv_object.2 ~ Camouflage.mismatch, data = dat.dect)
fit5
summary(fit5)

ggsurvplot(fit5, data = dat.dect, palette = "npg", conf.int = T, linetype = 1,
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
                   data = dat.dect)
summary(fit.coxph.4)
fit.coxph.4 %>%
  gtsummary::tbl_regression(exp = T)

#survival difference for mismatch / match within certain habitat
dat.dect$cxh <- paste(dat.dect$Camouflage.mismatch, dat.dect$habitat)

dat.dect$cxh <- factor(dat.dect$cxh, levels = c("Match O","Mismatch O",
                                                "Match C", "Mismatch C"))

fit6 <- survfit(surv_object.2 ~ cxh, data = dat.dect)
fit6
summary(fit6)

ggsurvplot(fit6, data = dat.dect, palette = "npg", conf.int = T, linetype = 1,
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

fit.coxph.5 <- coxph(surv_object.2 ~ cxh, data = dat.dect)
summary(fit.coxph.3)
fit.coxph.5 %>%
  gtsummary::tbl_regression(exp = T)

####General cam data####
#species vs habitat / treatment table
dat.13 <- dat.3 %>%
  dplyr::group_by(deploymentID, cluster_id, round, week, habitat, cluster,
                  treatment, commonName,sequenceID,
                  Camouflage.mismatch) %>%
  dplyr::summarise(n = max(count)) %>%
  dplyr::filter(commonName != "")
dat.13 <- dat.13 %>%
  filter(commonName != "Human")
dat.13$Camouflage.mismatch[dat.13$Camouflage.mismatch == ""] <- "control"

dat.13 <- dat.13 %>%
  dplyr::group_by(habitat, treatment, Camouflage.mismatch, commonName) %>%
  dplyr::summarise(n = sum(n))

dat.13 <- dat.13 %>%
  pivot_wider(names_from = commonName, values_from = n)

dat.13[is.na(dat.13)] <- 0

dat.13 <- dat.13[c(1:3,10,5,6,7,4,8,9)]

dat.13 <- dat.13 %>%
  arrange(factor(treatment, levels = c("C", "B", "W"))) %>%
  arrange(habitat)

dat.13$sum = rowSums(dat.13[,c(4:10)])

library(writexl)
write_xlsx(dat.13,"summary_pub.xlsx")

#general cam data:
dat.14 <- dat %>%
  group_by(observationType) %>%
  summarise(n = n())
dat.14[2,2] <- dat.14[2,2] + dat.14[4,2]
dat.14 <- dat.14[-4,]

dat.14[2,2] <- dat.14[2,2] + 52195

sum(dat.14$n)

dat.14$percentage <- ((dat.14$n / sum(dat.14$n) * 100))

dat.15 <- dat %>%
  filter(observationType == "animal") %>%
  group_by(predator) %>%
  summarise(n = n())
dat.15$percentage <- ((dat.15$n / sum(dat.15$n) * 100))

####snow cover####
dat.snow <- dat %>%
  filter(observationType == "animal")

dat.snow <- dat.snow %>%
  separate(timestamp, c("datetime", "GMT"), sep = 19) %>%
  separate(datetime, into = c("date", "time"), sep = "T", remove = FALSE)
dat.snow$date <- as.Date(dat.snow$date )

dat.snow <- dat.snow %>%
  separate(deploymentID, c("round", "deployment"), remove = F) %>%
  separate(deployment, c("week", "treatment"), sep = 3) %>%
  separate(treatment, c("habitat", "treatment"), sep = 1) %>%
  separate(treatment, c("cluster", "treatment"), sep = 2) %>%
  separate(deploymentID, c("cluster_id"), sep = 9, remove = F)

dat.snow <- dat.snow %>%
  group_by(deploymentID, date, habitat) %>%
  summarise(Snow.cover = max(Snow.cover))

unique(dat.snow$deploymentID)

windows()
ggplot() +
  geom_jitter(data = dat.snow, aes(x = date, y = Snow.cover, color = habitat),
              width = 0.6, height = 0.1, size = 2) +
  geom_smooth(data = dat.snow, aes(x = date, y = Snow.cover, color = habitat),
              method = "loess", size = 3) + 
  scale_y_continuous(name = "Snow cover (%)", limits = c(0,100)) + 
  scale_x_date(name = "Date", date_labels = "%d-%m") + 
  scale_color_npg(name = "Habitat", labels = c("Closed", "Open")) + 
  theme_bw() +
  theme(text = element_text(size=30))

####END####

#transparant background
windows()
p <- ggsurvplot(fit1, data = dat.11, palette = "npg", conf.int = T, linetype = 1,
                conf.int.aplha = 0.7,pval = T, pval.method = T, surv.median.line = "hv", 
                risk.table = "percentage",
                xlab = "Time (days)", legend.title = "Camouflage",
                legend.labs = c("Match", "Mismatch"), break.time.by = 2,
                ylab = "Decoy survival probability",
                font.x = 20, font.y = 20, font.tickslab = 20, font.legend = 20,
                fontsize= 7,
                ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                plot.background = element_rect(fill='transparent',
                                                               color=NA),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                legend.background = element_rect(fill='transparent'),
                                legend.box.background = element_rect(fill='transparent')))


grid.draw.ggsurvplot <- function(x) survminer:::print.ggsurvplot(x, newpage = FALSE)

ggsave("test_survplot.png", p, bg = "transparent", width = 12, height = 7)

p1 <- effect_plot(m1, pred = Camouflage.mismatch.2, interval = T,
                  x.label = "Camouflage", y.label = "Predators (n)",
                  data = dat.6, colors = c( '#0072b2', '#FFC107', '#d55e00'),
                  line.thickness = 3, cat.pred.point.size = 5) + 
  theme_nice() +
  theme(text = element_text(size=30), panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent',
                                       color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))
)

ggsave("test_effect.png", p1, bg = "transparent", width = 12, height = 7)



