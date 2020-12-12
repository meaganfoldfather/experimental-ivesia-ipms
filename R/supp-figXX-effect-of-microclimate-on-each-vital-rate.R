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


dir.create("data/data_output", recursive = TRUE, showWarnings = FALSE)

# remote_source is where we can find the vital rate models to download
remote_source <- "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms"
# remote_target is where we put the lambda estimates from the IPM
remote_target <- "s3://earthlab-mkoontz/experimental-ivesia-ipms"

# Desired filenames for the new vital rate models
# If they are the same as previously-uploaded .rds files at the `remote_target`,
# AND the `overwrite` variable is TRUE (it is FALSE by default), then those
# vital rate models will be overwritten

surv_mod_fname <- "surv_mod_additive.rds"
growth_mod_fname <- "grwth_mod_additive.rds"
establishment_mod_fname <- "recruit_mod_additive.rds"
hurdle_mod_fname <- "hurdle_mod_additive.rds"

lambda_df_fname <- "mc_vr_effects_additive.csv"

overwrite <- FALSE

##### Get vital rate models from remote_source ####

if(!file.exists(file.path("data/data_output", surv_mod_fname))) {
  download.file(url = glue::glue("{remote_source}/{surv_mod_fname}"), 
                destfile = file.path("data/data_output", surv_mod_fname))
}

if(!file.exists(file.path("data/data_output", growth_mod_fname))) {
  download.file(url = glue::glue("{remote_source}/{growth_mod_fname}"),
                destfile = file.path("data/data_output", growth_mod_fname))
}

if(!file.exists(file.path("data/data_output", establishment_mod_fname))) {
  
  download.file(url = glue::glue("{remote_source}/{establishment_mod_fname}"),
                destfile = file.path("data/data_output", establishment_mod_fname))
}

if(!file.exists(file.path("data/data_output", hurdle_mod_fname))) {
  download.file(url = glue::glue("{remote_source}/{hurdle_mod_fname}"),
                destfile = file.path("data/data_output", hurdle_mod_fname))
}

# If the files STILL don't exist on disk, then they weren't available from the
# remote source, so we have to run the script to build the vital rate models
# using the filenames declared above
if(!file.exists(file.path("data/data_output", surv_mod_fname)) |
   !file.exists(file.path("data/data_output", growth_mod_fname)) |
   !file.exists(file.path("data/data_output", hurdle_mod_fname)) |
   !file.exists(file.path("data/data_output", hurdle_mod_fname))) {
  
  source("R/CH2_VRModels_Revised.R")
}

survival_model <- readr::read_rds(file.path("data/data_output", surv_mod_fname))
growth_model <- readRDS(file.path("data/data_output", growth_mod_fname))
establishment_model <- readRDS(file.path("data/data_output", establishment_mod_fname))
hurdleRep <- readRDS(file.path("data/data_output", hurdle_mod_fname))

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
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .1,lwd = .01)+
  xlab("Degree-days")+
  ylab("Survival Probability")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())
vr1

# vr1 <- surv %>% 
#   filter(vwc == 0) %>% 
# ggplot(aes(degree.days, Estimate, col = interaction(heat,water), fill = interaction(heat,water)))+
#   geom_smooth(se=F, alpha = .5)+
#   geom_ribbon(aes(ymin = Estimate - Est.Error, ymax = Estimate + Est.Error), alpha = .1,lwd = .01)+
#   xlab("Degree-days")+
#   ylab("Survival Probability")+
#   theme_bw()+
#   scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
#   scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
#   theme(text = element_text(size=14), legend.title = element_blank())
# vr1

#vwc
vr2<- surv %>% 
  filter(degree.days == 0) %>% 
  ggplot(aes(vwc, Estimate, col = interaction(heat,water), fill = interaction(heat,water)))+
  geom_smooth(se=F)+
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .1, lwd = .1)+
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
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .1, lwd = .1)+
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
    geom_smooth(se=F)+
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .1, lwd = .1)+
  xlab("Soil Moisture")+
  ylab("Standerdized Growth")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())+
  geom_abline(slope = 0, intercept = 0, color = "darkgrey") 
plot_grid(vr3, vr4, nrow = 2)

#### make reproduction predictions ####
rep_newdata <- data.frame(expand.grid(sizeNext.s = 0, degree.days = seq(to = -1.5, from = 2, length.out = 71) , vwc = seq(to = -1.5, from = 1, length.out = 51) , heat = c(0,1), water = c(0,1), t1 = NA))
rep_pred <- predict(hurdleRep, newdata = rep_newdata, re_formula = NA)
rep <- cbind(rep_newdata, rep_pred)
head(rep)
#dd
vr5 <- rep %>% 
  filter(vwc == 0) %>% 
ggplot(aes(degree.days, Estimate, color = interaction (heat, water), fill = interaction (heat, water)))+
     geom_smooth(se=F)+
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .1, lwd = .1)+
  xlab("Degree-days")+
  ylab("Reproduction")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())
vr5

#vwc
vr6<- rep %>% 
  filter(degree.days == 0) %>% 
ggplot(aes(vwc, Estimate, color = interaction (heat, water), fill = interaction (heat, water)))+
      geom_smooth(se=F)+
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .1, lwd = .1)+
  xlab("Soil Moisture")+
   ylab("Reproduction")+
  theme_bw()+
  scale_color_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  scale_fill_manual(values = c("black", "red", "dodgerblue", "purple"), labels = c("Ambient", "Heat", "Water", "Heat+Water"))+
  theme(text = element_text(size=14), legend.title = element_blank())
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
