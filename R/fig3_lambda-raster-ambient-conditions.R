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
  summarize(mean_lambda = mean(lambda),
            lwr95 = quantile(lambda, probs = 0.025),
            upr95 = quantile(lambda, probs = 0.975),
            mass_over_under_one = ifelse(mean_lambda < 1,
                                         yes = ecdf(lambda)(1),
                                         no = 1 - ecdf(lambda)(1)),
            sig_lambda = ifelse(mass_over_under_one < 0.90,
                                yes = NA, no = mean_lambda))

mcvr_lambda_summary$heat <- factor(mcvr_lambda_summary$heat)
mcvr_lambda_summary$water <- factor(mcvr_lambda_summary$water)

# make  lines around ambient cells
r <- 
  mcvr_lambda_summary %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(heat == 0 & water == 0) %>% 
  dplyr::mutate(expanding = ifelse(sig_lambda > 1, yes = 1, no = 0)) %>% 
  dplyr::select(degree.days, vwc, expanding) %>% 
  raster::rasterFromXYZ()

crs(r) <- sf::st_crs(4326)

pp <- raster::rasterToPolygons(r, dissolve = TRUE)
outline <- sf::st_as_sf(pp) %>% dplyr::filter(expanding == 1)

### Get the site metadata
site_metadata <- readr::read_csv(file = "data/data_output/site_metadata.csv")

# Build the figure
figure3 <- 
  mcvr_lambda_summary %>% 
  filter(heat == 0,
         water == 0) %>% 
  ggplot(aes(x = degree.days, y = vwc)) +
  geom_raster(aes(fill = sig_lambda)) +
  scale_fill_gradient2(mid = "grey90", 
                       midpoint = 1, 
                       na.value = "grey90", 
                       high = "darkgreen", 
                       low = "darkgoldenrod1")+
  theme_classic() +
  guides(alpha = FALSE) +
  xlab("Degree-Days (scaled)") +
  ylab("Soil Moisture (scaled)") +
  theme(text = element_text(size=16)) +
  labs(fill = expression(lambda)) +
  geom_point(data = site_metadata, aes(degree.days_scaled, vwc_scaled), pch = 23, cex = 8) +
  geom_text(data = site_metadata, aes(degree.days_scaled, vwc_scaled, label = site.num), cex = 5) +
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA)

figure3
