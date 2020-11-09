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
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient"))

# site-specific lambdas per elevation
ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("ambient")), aes(x = elevation, y = lambda)) +
  tidybayes::stat_halfeye(limits = c(0,2))+
  theme_classic()

ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("heat", "water", "ambient", "hw")), aes(x = lambda, y = elevation, fill = manipulated_vr)) +
  tidybayes::stat_halfeye()

ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("heat", "water", "ambient", "hw")), aes(x = elevation, y = lambda, col = manipulated_vr)) +
  tidybayes::stat_halfeye() +
  facet_wrap(facets = "manipulated_vr")+
  theme_classic()

#### Contrasts ####

#### Contrasts of experimental treatments against ambient ####
contrasts <- 
  site_lambdas %>% 
  filter(manipulated_vr %in% c("ambient", "heat", "water", "hw")) %>% 
  mutate(heat = ifelse(heat.f == 1, yes = 1, no = 0), 
         water = ifelse(water.f == 1, yes = 1, no = 0)) %>% 
  dplyr::select(-heat.f:-water.s, - manipulated_vr) %>% 
  tidyr::unite(trt, heat, water) %>% 
  tidyr::nest(lambda = lambda) %>% 
  tidyr::pivot_wider(names_from = trt, values_from = lambda) %>%
  tidyr::unnest(cols = c(`0_0`, `1_0`, `0_1`, `1_1`), names_sep = "_") %>% 
  dplyr::rename(ambient = `0_0_lambda`, 
                heat = `1_0_lambda`, 
                water = `0_1_lambda`, 
                hw = `1_1_lambda`) %>%
  dplyr::mutate(heat_effect = (heat+hw)-(ambient+water), 
                water_effect = (water+hw) - (ambient + heat), 
                hw_effect = hw - (1/3)*(ambient + heat + water)) %>% 
  dplyr::select(-heat, -water, -ambient, -hw) %>% 
  tidyr::pivot_longer(names_to = "trt", values_to = "delta_lambda", -(site:elevation)) 

contrasts

contrasts_summary <-
  contrasts %>% 
  dplyr::group_by(site, elevation, degree.days, vwc, trt) %>% 
  dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
                   median_delta_lambda = median(delta_lambda),
                   lwr95 = quantile(delta_lambda, probs = 0.025),
                   upr95 = quantile(delta_lambda, probs = 0.975),
                   mass_over_under_zero = ifelse(mean_delta_lambda < 0,
                                                 yes = ecdf(delta_lambda)(0),
                                                 no = 1 - ecdf(delta_lambda)(0)),
                   sig_lambda = ifelse(mass_over_under_zero < 0.95,
                                       yes = NA, no = mean_delta_lambda), 
                   n=n())


ggplot(contrasts_summary, aes(x = elevation, y = mean_delta_lambda)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr95, ymax = upr95)) +
  facet_wrap(.~trt, nrow = 1, labeller = labeller(trt = c(heat_effect  = "HEAT", water_effect= "WATER", hw_effect = "HEAT + WATER")))+
  theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90))+
  xlab("Elevation")+
  ylab(expression(paste(Delta, lambda)))+
  geom_hline(yintercept = 0, color = "grey")


#### Contrasts of simulated treatment on individual vital rates against ambient ####


