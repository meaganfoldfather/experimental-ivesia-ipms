#Intent: CH2 - Look at vital rate responses to experimental treatment
#20200507

#### Library ####
library(lme4)
library(lmerTest)
library(ggplot2)
library(car)
library(popdemo)
library(dplyr)
library(tidyr)
library(viridis)
library(MASS)
library(MuMIn)
library(rstan)
library(brms)
library(purrr)
library(broom)
library(devtools)
library(furrr)
library(tidyverse)
library(scales)
library(cowplot)


#### Bring in brms model for the vital rates ####
#load in models
# read in models
survival_model <- readRDS("surv_mod_additive.rds")
growth_model <- readRDS("grwth_mod_additive.rds")
establishment_model <- readRDS("recruit_mod_additive.rds")
hurdleRep <- readRDS("hurdle_mod_additive.rds")

#### make survival predictions ####
surv_newdata <- data.frame(expand.grid(size.s = 0, degree.days = seq(to = -1.5, from = 2, length.out = 71) , vwc = seq(to = -1.5, from = 1, length.out = 51),  heat = c(0,1), water = c(0,1), t1 = NA))
surv_pred <- fitted(survival_model, newdata = surv_newdata, re_formula = NA)
surv <- cbind(surv_newdata, surv_pred)
head(surv)
dim(surv)

#dd
vr1 <- surv %>% 
  filter(vwc == 0) %>% 
ggplot(aes(degree.days, Estimate, col = interaction(heat,water), fill = interaction(heat,water)))+
  geom_smooth(se=F, alpha = .5)+
  geom_ribbon(aes(ymin = Estimate - Est.Error, ymax = Estimate + Est.Error), alpha = .1,lwd = .01)+
  xlab("Degree-days")+
  ylab("Survival Probability")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())
vr1

#vwc
vr2<- surv %>% 
  filter(degree.days == 0) %>% 
ggplot(aes(vwc, Estimate, col = interaction(heat,water), fill = interaction(heat,water)))+
  geom_smooth(se=F)+
    geom_ribbon(aes(ymin = Estimate - Est.Error, ymax = Estimate + Est.Error), alpha = .1, lwd = .1)+
  xlab("Soil Moisture")+
  ylab("Survival Probability")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
    scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())
vr2
plot_grid(vr1, vr2, nrow = 2)

#### make growth predictions ####
grwth_newdata <- data.frame(expand.grid(size.s = 0, degree.days = seq(to = -1.5, from = 2, length.out = 71) , vwc = seq(to = -1.5, from = 1, length.out = 51) ,heat = c(0,1), water = c(0,1),t1 = NA))
grwth_pred <- predict(growth_model, newdata = grwth_newdata, re_formula = NA)
grwth <- cbind(grwth_newdata, grwth_pred)
head(grwth)
#dd
vr3 <- grwth %>% 
  filter(vwc == 0) %>% 
ggplot(aes(degree.days, Estimate, col = interaction(heat,water), fill = interaction(heat,water)))+
  geom_smooth(se=F)+
    geom_ribbon(aes(ymin = Estimate - Est.Error, ymax = Estimate + Est.Error), alpha = .1, lwd = .1)+
  xlab("Degree-days")+
  ylab("Standerdized Growth")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())+
  geom_abline(slope = 0, intercept = 0, color = "darkgrey") 
vr3

#vwc
vr4<- grwth %>% 
  filter(degree.days == 0) %>% 
ggplot(aes(vwc, Estimate, color = interaction (heat, water), fill = interaction (heat, water)))+
  geom_smooth(se=T)+
  xlab("Soil Moisture")+
  ylab("Standerdized Growth")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())+
  ylim(-.25, .25)+
  geom_abline(slope = 0, intercept = 0, color = "darkgrey") 
plot_grid(vr3, vr4, nrow = 2)

#### make establishment predictions --> not working ####
est_newdata <- data.frame(expand.grid(seeds.avaliable = 100, degree.days = seq(to = -2, from = 2, length.out = 51) , vwc = seq(to = -2, from = 2, length.out = 51) , heat = c(0,1), water = c(0,1), ))
est_pred <- predict(establishment_model, newdata = est_newdata, re_formula = NA)
est<- cbind(est_newdata, est_pred)
head(est)
#dd
vr1 <- est %>% 
  filter(vwc == 0) %>% 
ggplot(aes(degree.days, Estimate, color = interaction (heat, water)))+
  geom_smooth(se=F)+
  xlab("Degree-days")+
  ylab("Establishment Probability")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=16), legend.title = element_blank())
  #ylim(.6, 1)
vr1
#vwc
vr2<- est %>% 
  filter(degree.days == 0) %>% 
ggplot(aes(vwc, Estimate, color = interaction (heat, water)))+
  geom_smooth(se=F)+
  xlab("Soil Moisture")+
  ylab("Establishment Probability")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=16), legend.title = element_blank())
  #ylim(.6, 1)
plot_grid(vr1, vr2, nrow = 2)

#### make reproduction predictions ####
rep_newdata <- data.frame(expand.grid(sizeNext.s = 0, degree.days = seq(to = -1.5, from = 2, length.out = 71) , vwc = seq(to = -1.5, from = 1, length.out = 51) , heat = c(0,1), water = c(0,1), t1 = 2015:2017))
rep_pred <- predict(hurdleRep, newdata = rep_newdata, re_formula = NA)
rep <- cbind(rep_newdata, rep_pred)
head(rep)
#dd
vr5 <- rep %>% 
  filter(vwc == 0) %>% 
ggplot(aes(degree.days, Estimate, color = interaction (heat, water), fill = interaction (heat, water)))+
  geom_smooth(se=T)+
  xlab("Degree-days")+
  ylab("Reproduction")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())+
  ylim(0, 20)
vr5

#vwc
vr6<- rep %>% 
  filter(degree.days == 0) %>% 
ggplot(aes(vwc, Estimate, color = interaction (heat, water), fill = interaction (heat, water)))+
  geom_smooth(se=T)+
  xlab("Soil Moisture")+
   ylab("Reproduction")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())+
  ylim(0, 20)
vr6


#### Make full plot ####
pdf(file = "vr.pdf")
prow <- plot_grid(
vr1 + theme(legend.position="none"),
vr2 + theme(legend.position="none"),
vr3 + theme(legend.position="none"),
vr4 + theme(legend.position="none"),
vr5 + theme(legend.position="none"),
vr6 + theme(legend.position="none"), nrow = 3)
legend_b <- get_legend(vr1 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
dev.off()
