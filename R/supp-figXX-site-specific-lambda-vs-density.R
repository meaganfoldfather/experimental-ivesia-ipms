
#### Dependencies ####
library(dplyr)
library(brms)
library(purrr)
library(furrr)
library(glue)
library(ggplot2)
library(data.table)
library(tidybayes)

dir.create("data/data_output", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/", showWarnings = FALSE)

# remote_source is where we can find the vital rate models to download
remote_source <- "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms"
# remote_target is where we put the lambda estimates from the IPM
remote_target <- "s3://earthlab-mkoontz/experimental-ivesia-ipms"

lambda_df_fname <- "site_specific_additive.csv"

#### Get data from remote source ####

if(!file.exists(file.path("data/data_output", lambda_df_fname))) {
  download.file(url = glue::glue("{remote_source}/{lambda_df_fname}"), 
                destfile = file.path("data/data_output", lambda_df_fname))}

### Get the site metadata
site_metadata <- readr::read_csv(file = "data/data_output/site_metadata.csv")

# # Get site elevation
# plot.chars <- read.csv(glue::glue("{remote_source}/trt.plot.data.csv"))
# colnames(plot.chars) <- tolower(colnames(plot.chars))
# 
# site_chars <- 
#   plot.chars %>% 
#   dplyr::group_by(site) %>% 
#   dplyr::summarize(elevation = mean(elevation))

site_lambdas <- 
  readr::read_csv(file.path("data/data_output", lambda_df_fname)) %>% 
  dplyr::left_join(site_metadata, by = "site") %>% 
  dplyr::mutate(elevation = round(elevation_raw, digits = 0)) %>% 
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient"))

site_lambdas_summary <-
  site_lambdas %>% 
  filter(manipulated_vr %in% c("ambient", "heat", "water", "hw")) %>% 
  group_by(manipulated_vr, site) %>% 
  tidybayes::median_hdci(lambda)

#### Lambda across population density ####
# calculate average density across the plots
all <- read.csv(glue::glue("{remote_source}/VitalRates_Microclimate.csv"), stringsAsFactors = FALSE)
colnames(all)[1] <- "num.id"
head(all)

density <- 
  all %>% 
  group_by(site, t1, plot) %>%
  summarize(density = n_distinct(num.id)) %>% 
  group_by(site) %>% 
  summarize(density = round(mean(density))) 
density

# merge density with 'site_lambda
site_density <-left_join(site_lambdas_summary, density)

supp_figXX <- 
  ggplot(site_density, aes(x = density, y = lambda, ymin = .lower, ymax = .upper, color = manipulated_vr)) +
  geom_point() +
  geom_errorbar() + 
  labs(x = "Site Density",
       y = expression(lambda)) +
  geom_hline(yintercept = 1, color = "grey")+
  scale_color_manual(values=c("black", "red" ,"purple", "dodgerblue"))+
  theme_classic()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90), legend.position = "none")+
  facet_wrap(facets = "manipulated_vr", labeller = labeller(manipulated_vr = c(ambient = "AMBIENT", heat = "HEAT", water = "WATER", hw = "HEAT + WATER")))

ggsave(plot = supp_figXX, filename = "figs/supp-figXX-site-specific-lambda-versus-density.png")
