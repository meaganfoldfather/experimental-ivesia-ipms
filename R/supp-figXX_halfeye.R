# Make half-eye plots 

#### Library & Remote Info ####
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
library(tidybayes)

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
  
  source("R/build-vital-rate-models.R")
}
####







survival_model <- readr::read_rds(file.path("data/data_output", surv_mod_fname))
growth_model <- readRDS(file.path("data/data_output", growth_mod_fname))
establishment_model <- readRDS(file.path("data/data_output", establishment_mod_fname))
hurdleRep <- readRDS(file.path("data/data_output", hurdle_mod_fname))

#### Make half-eye plots
# survival 
surv_samps <- posterior_samples(survival_model) %>% 
   dplyr::select(starts_with("b_")) %>% 
   pivot_longer(cols = starts_with("b_"), names_to = "coefficient", values_to = "estimate")
 
 ordered_surv_samps <- 
   surv_samps %>% 
   mutate(colon.count = stringi::stri_count(coefficient, fixed = ":")) %>% 
   mutate(coefficient = stringi::stri_sub(str = coefficient, from = 3, to = -1)) %>% 
   arrange(colon.count) %>% 
   mutate(coefficient = factor(coefficient, levels = rev(unique(coefficient))))
   
halfeye_survival <- ggplot(ordered_surv_samps, aes(y = coefficient, x = estimate))+
  stat_halfeye()+
  geom_vline(aes(xintercept = 0), lty = "dashed")+
  theme_minimal() +
  labs(title = "Survival")+
  xlab("Estimate")+
  ylab("Coefficient")
halfeye_survival

ggsave(halfeye_survival, filename = "figs/supp-figXX-halfeye-survival.png", width = 10, height = 12, units = "in")

# growth
grwth_samps <- posterior_samples(growth_model) %>% 
   dplyr::select(starts_with("b_")) %>% 
   pivot_longer(cols = starts_with("b_"), names_to = "coefficient", values_to = "estimate")

  ordered_grwth_samps <- 
   grwth_samps %>% 
   mutate(colon.count = stringi::stri_count(coefficient, fixed = ":")) %>% 
   mutate(coefficient = stringi::stri_sub(str = coefficient, from = 3, to = -1)) %>% 
   arrange(colon.count) %>% 
   mutate(coefficient = factor(coefficient, levels = rev(unique(coefficient))))
   
halfeye_growth<-ggplot(ordered_grwth_samps, aes(y = coefficient, x = estimate))+
  stat_halfeye()+
  geom_vline(aes(xintercept = 0), lty = "dashed")+
  theme_minimal() +
  labs(title = "Growth")+
  xlab("Estimate")+
  ylab("Coefficient")
halfeye_growth

ggsave(halfeye_growth, filename = "figs/supp-figXX-halfeye-growth.png", width = 10, height = 12, units = "in")

#reproduction
rep_samps <- posterior_samples(hurdleRep) %>% 
   dplyr::select(starts_with("b_")) %>% 
   pivot_longer(cols = starts_with("b_"), names_to = "coefficient", values_to = "estimate")

ordered_rep_samps <- 
  rep_samps %>% 
   mutate(colon.count = stringi::stri_count(coefficient, fixed = ":")) %>% 
  mutate(hu.count = stringi::stri_count(coefficient, fixed = "hu")) %>% 
   mutate(coefficient = stringi::stri_sub(str = coefficient, from = 3, to = -1)) %>% 
   arrange(hu.count, colon.count) %>% 
   mutate(coefficient = factor(coefficient, levels = rev(unique(coefficient))))

halfeye_rep<- ggplot(ordered_rep_samps, aes(y = coefficient, x = estimate))+
  stat_halfeye()+
  geom_vline(aes(xintercept = 0), lty = "dashed")+
  theme_minimal() +
  labs(title = "Reproduction")+
     xlab("Estimate")+
  ylab("Coefficient")
halfeye_rep

ggsave(halfeye_rep, filename = "figs/supp-figXX-halfeye-rep.png", width = 10, height = 12, units = "in")

#establishment
est_samps <- posterior_samples(establishment_model) %>% 
   dplyr::select(starts_with("b_")) %>% 
   pivot_longer(cols = starts_with("b_"), names_to = "coefficient", values_to = "estimate")

ordered_est_samps <- 
   est_samps %>% 
   mutate(colon.count = stringi::stri_count(coefficient, fixed = ":")) %>% 
   mutate(coefficient = stringi::stri_sub(str = coefficient, from = 3, to = -1)) %>% 
   arrange(colon.count) %>% 
   mutate(coefficient = factor(coefficient, levels = rev(unique(coefficient))))

halfeye_establishment <- ggplot(ordered_est_samps, aes(y = coefficient, x = estimate))+
  stat_halfeye()+
  geom_vline(aes(xintercept = 0), lty = "dashed")+
  theme_minimal() +
  labs(title = "Establishment")+
    xlab("Estimate")+
  ylab("Coefficient")
halfeye_establishment 

ggsave(halfeye_establishment , filename = "figs/supp-figXX-halfeye-establishment.png", width = 10, height = 12, units = "in")

