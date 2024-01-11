rm(list=ls())

####Load data####
#install needed packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(jtools)
library(multcomp)
library(ggsci)
library(rcompanion)
library(cowplot)

#load data
dat <- read.csv(file="dat.pred.csv", header = T)

#test for predator difference between habitats
dat.control <- dat %>%
  filter(Camouflage.mismatch.2 == "Control")
m1 <- glm(n ~ habitat + offset(log(cam_effort)), data = dat.control, family = poisson)
summary(m1)

#Run model with desired variables
m2 <- glmer(n ~ Camouflage.mismatch.2 + habitat + offset(log(cam_effort)) + (1|cluster_id), data = dat, family = poisson)
summary(m2)

#tukey HSD test for differences between groups
m1a <- glht(m2, linfct= mcp(Camouflage.mismatch.2 = "Tukey"))
summary(m1a) #significant difference between groups B & W but not between rest

m2a <- glht(m2, linfct= mcp(habitat = "Tukey"))
summary(m2a)

#filter C out of data
dat.cam <- dat %>%
  filter(treatment != "C")

#test for interaction between treatment & snow.cover
dat.cam$snow.cover <- factor(dat.cam$snow.cover, c(1,0))
dat.cam$treatment <- factor(dat.cam$treatment, c("B", "W"))
m3 <- glmer(n ~ treatment + treatment:snow.cover +
              offset(log(cam_effort)) + (1|cluster_id), data = dat.cam, family = poisson)

summary(m3)

#plot predicted values based on GLMM
p1 <- effect_plot(m2, pred = Camouflage.mismatch.2, interval = T,
                  x.label = "Camouflage", y.label = "Predators (n)",
                  data = dat, colors = c( '#0072b2', '#FFC107', '#d55e00'),
                  line.thickness = 3, cat.pred.point.size = 5) + 
  theme_nice() +
  theme(text = element_text(size=30))

p2 <- cat_plot(m3, data = dat.cam, pred = treatment, modx = snow.cover, interval = T, colors = c('#0072b2', '#FFC107'),
               x.label = "Decoy colour", y.label = "Predators (n)", pred.labels = c("Brown", "White"), modx.labels = c("Snow\npresent", "Snow\nabsent"),
               legend.main = "Snow\ncover", line.thickness = 3, pred.point.size = 5) + theme_nice() + theme(text = element_text(size=30))

windows()
plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, label_size = 20)

#test for difference between predator type
dat.pred.type <- read.csv(file = "dat.pred.type.csv", header = T)

m4 <- glmer(n~Camouflage.mismatch.2 + habitat + pred.type + offset(log(cam_effort)) + (1|cluster_id), data = dat.pred.type, family = poisson)
summary(m4)

#####GLMM's for interaction####
dat.int <- read.csv(file="dat.interaction.csv", header = T)

dat.int.2 <- dat.int %>%
  filter(treatment != "C" & commonName != "" & commonName != "Human")

m4 <- glmer(interaction ~ Camouflage.mismatch + (1|cluster_id), 
            data = dat.int.2, family = binomial)
summary(m4)

#separate mammalian and avian predators
dat.int.mam <- dat.int %>%
  filter(mammal == 1) %>%
  filter(treatment != "C")
dat.int.mam$Camouflage.mismatch <- as.character(dat.int.mam$Camouflage.mismatch)
dat.int.mam$Camouflage.mismatch <- as.factor(dat.int.mam$Camouflage.mismatch)

dat.int.av <- dat.int %>%
  filter(mammal == 0) %>%
  filter(treatment != "C") %>%
  filter(commonName != "" & commonName != "Human")
dat.int.av$Camouflage.mismatch <- as.character(dat.int.av$Camouflage.mismatch)
dat.int.av$Camouflage.mismatch <- as.factor(dat.int.av$Camouflage.mismatch)

m5 <- glmer(interaction ~ Camouflage.mismatch + (1|cluster_id), 
            data = dat.int.mam, family = binomial)
summary(m5)

m6 <- glmer(interaction ~ Camouflage.mismatch + (1|cluster_id), 
            data = dat.int.av, family = binomial)
summary (m6)
####end####


