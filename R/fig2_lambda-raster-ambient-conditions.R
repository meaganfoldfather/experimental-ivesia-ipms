# Figures 3; lambda across ambient microclimate gradient of degree days and VWC

#### Dependencies ####
library(dplyr)
library(ggplot2)
library(brms)
library(purrr)
library(furrr)
library(data.table)
library(glue)
library(tidyr)
library(cowplot)
library(raster)
library(broom)
library(ggrepel)
library(tidybayes)
library(rgeos)

# remote_source is where we can find the vital rate models to download
remote_source <- "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms"
# remote_target is where we put the lambda estimates from the IPM
remote_target <- "s3://earthlab-mkoontz/experimental-ivesia-ipms"

# Desired filenames for the new vital rate models
# If they are the same as previously-uploaded .rds files at the `remote_target`,
# AND the `overwrite` variable is TRUE (it is FALSE by default), then those
# vital rate models will be overwritten

# surv_mod_fname <- "surv_mod_additive.rds"
# growth_mod_fname <- "grwth_mod_additive.rds"
# establishment_mod_fname <- "recruit_mod_additive.rds"
# hurdle_mod_fname <- "hurdle_mod_additive.rds"

lambda_df_fname <- "mc_vr_effects_additive.csv"

if(!file.exists(glue::glue("data/data_output/{lambda_df_fname}"))) {
  
  download.file(url = glue::glue("{remote_source}/{lambda_df_fname}"),
                destfile = glue::glue("data/data_output/{lambda_df_fname}"))
}

mcvr <- readr::read_csv(glue::glue("data/data_output/{lambda_df_fname}"))

mcvr_lambda_summary <- 
  mcvr %>% 
  filter(manipulated_vr %in% c("ambient", "heat", "water", "hw")) %>% 
  mutate(heat = ifelse(heat.f == 1, yes = 1, no = 0), 
         water = ifelse(water.f == 1, yes = 1, no = 0)) %>% 
  group_by(degree.days, vwc, heat, water) %>% 
  tidybayes::median_hdci(lambda) %>% 
  dplyr::mutate(sig_lambda = ifelse(test = .lower >= 1 | .upper <= 1,
                                    yes = lambda,
                                    no = NA))

mcvr_lambda_summary$heat <- factor(mcvr_lambda_summary$heat)
mcvr_lambda_summary$water <- factor(mcvr_lambda_summary$water)

# make  lines around ambient cells
r <- 
  mcvr_lambda_summary %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(heat == 0 & water == 0) %>% 
  dplyr::mutate(stable = is.na(sig_lambda)) %>% 
  dplyr::select(degree.days, vwc, stable) %>% 
  raster::rasterFromXYZ()

crs(r) <- sf::st_crs(4326)

pp <- raster::rasterToPolygons(r, dissolve = TRUE)
outline <- sf::st_as_sf(pp) %>% dplyr::filter(stable == 1)

### Get the site metadata
site_metadata <- readr::read_csv(file = "data/data_output/site_metadata.csv")

# bring in site lambda values
### Get the site metadata
site_lambdas<- readr::read_csv(file = "data/data_output/site_specific_additive.csv")
site_lambdas
dim(site_lambdas)

ambient_site_lambdas <- 
  site_lambdas %>% 
  dplyr::left_join(site_metadata, by = "site") %>% 
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient")) %>% 
  dplyr::filter(manipulated_vr == 'ambient') %>% 
  group_by(site, degree.days, vwc, site.num) %>% 
  tidybayes::median_hdci(lambda) %>% 
  dplyr::mutate(sig_lambda = ifelse(test = .lower >= 1 | .upper <= 1,
                                    yes = lambda,
                                    no = NA))


ambient_site_lambdas


# Build the figure
figure2 <-
  mcvr_lambda_summary %>%
  filter(heat == 0,
         water == 0) %>%
  ggplot(aes(x = degree.days, y = vwc)) +
  geom_raster(aes(fill = lambda)) +
  scale_fill_gradient2(mid = "grey90",
                       midpoint = 1,
                       na.value = "grey90",
                       high = "darkgreen",
                       low = "darkgoldenrod1")+
  theme_classic() +
  guides(alpha = FALSE) +
  xlab("Degree-Days (scaled)") +
  ylab("Soil Moisture (scaled)") +
  theme(text = element_text(size=12)) +
  labs(fill = expression(lambda)) +
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA, lwd = 1)+
  geom_point(data = ambient_site_lambdas[is.na(ambient_site_lambdas$sig_lambda),], aes(degree.days, vwc, bg = lambda), pch =23, cex = 2) +
    geom_point(data = ambient_site_lambdas[!is.na(ambient_site_lambdas$sig_lambda),], aes(degree.days, vwc, bg = lambda), pch =25, cex = 2) +
  geom_text_repel(data = ambient_site_lambdas, aes(degree.days, vwc, label = site.num), cex = 2, point.padding
= 1) +
   scale_color_gradient2(mid = "grey90",
                       midpoint = 1,
                       na.value = "grey90",
                       high = "darkgreen",
                       low = "darkgoldenrod1")

figure2

ggsave(plot = figure2, filename = "figs/fig2-ambient-lambda-across-microclimate-conditions.tiff", dpi = 600, width = 110, units = "mm")
