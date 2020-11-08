#Analyze site-specific IPM output and make a corresponding figure


#### Header ####
library(dplyr)
library(brms)
library(purrr)
library(furrr)
library(glue)
library(ggplot2)
library(data.table)
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

lambda_df_fname <- "site_specific_additive.csv"

overwrite <- TRUE

#### Get data from remote source ####

if(!file.exists(file.path("data/data_output", lambda_df_fname))) {
  download.file(url = glue::glue("{remote_source}/{lambda_df_fname}"), 
                destfile = file.path("data/data_output", lambda_df_fname))}

# Get site elevation
plot.chars <- read.csv(glue::glue("{remote_source}/trt.plot.data.csv"))
colnames(plot.chars) <- tolower(colnames(plot.chars))

site_chars <- 
  plot.chars %>% 
  dplyr::group_by(site) %>% 
  dplyr::summarize(elevation = mean(elevation))

site_lambdas <- 
  readr::read_csv(file.path("data/data_output", lambda_df_fname)) %>% 
  dplyr::left_join(site_chars, by = "site") %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(elevation = round(elevation,digits = 0)) %>% 
  dplyr::mutate(elevation = factor(elevation, levels = sort(unique(elevation))))


site_lambdas <-
  site_lambdas %>% 
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient"))
  
  # dplyr::mutate(trt = case_when(heat.s + heat.g + heat.f > 0 ~ "heat",
  #                               water.s + water.g + water.f > 0 ~ "water",
  #                               heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
  #                               TRUE ~ "ambient")) %>% 
  # dplyr::mutate(manipulated_vr = case_when(heat.f == 1 | water.f == 1 ~ "fecundity",
  #                                          heat.g == 1 | water.g == 1 ~ "growth",
  #                                          heat.s == 1 | water.s == 1 ~ "survivorship"))

# site-specific lambdas per elevation
ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("ambient")), aes(x = elevation, y = lambda)) +
  tidybayes::stat_halfeye(limits = c(0,2))+
  theme_classic()
  
ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("heat", "water", "ambient", "hw")), aes(x = lambda, y = elevation, fill = manipulated_vr)) +
  tidybayes::stat_halfeye()

ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("heat", "water", "ambient", "hw")), aes(x = lambda, y = elevation)) +
  tidybayes::stat_halfeye() +
  facet_wrap(facets = "manipulated_vr")

ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("heat", "water", "ambient", "hw")), aes(x = elevation, y = lambda, col = manipulated_vr)) +
  tidybayes::stat_halfeye() +
  facet_wrap(facets = "manipulated_vr")+
  theme_classic()


### change in lambdas by manipulating individual vital rates
ggplot(site_lambdas %>% dplyr::filter(!(manipulated_vr %in% c("heat", "water", "ambient", "hw"))), aes(x = lambda, y = elevation)) +
  tidybayes::stat_halfeye() +
  facet_wrap(facets = "manipulated_vr")

