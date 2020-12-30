# Figure 5c 

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

lambda_df_fname <- "mc_vr_effects_additive.csv"

if(!file.exists(glue::glue("data/data_output/{lambda_df_fname}"))) {
  
  download.file(url = glue::glue("{remote_source}/{lambda_df_fname}"),
                destfile = glue::glue("data/data_output/{lambda_df_fname}"))
}

mcvr <- readr::read_csv(glue::glue("data/data_output/{lambda_df_fname}"))

#### Use the ambient plots to make the outline around cells with an expanding 
# population
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
make_stable_outline <- function(data, h, w) {
  r <- 
    data %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(heat == h & water == w) %>% 
    dplyr::mutate(stable = is.na(sig_lambda)) %>% 
    dplyr::select(degree.days, vwc, stable) %>% 
    raster::rasterFromXYZ()
  
  crs(r) <- sf::st_crs(4326)
  
  pp <- raster::rasterToPolygons(r, dissolve = TRUE)
  stable_outline <- sf::st_as_sf(pp) %>% dplyr::filter(stable == 1)
  
  return(stable_outline)
}

ambient_stable_outline <- make_stable_outline(data = mcvr_lambda_summary, h = 0, w = 0)

#### End code that generates outline

# make contrast
# spread trt to look at contrasts
vr_contrasts <- 
  mcvr %>% 
  dplyr::filter(manipulated_vr %in% c("fecundity", "growth", "survivorship", "ambient")) %>%
  dplyr::select(- manipulated_vr) %>% 
  tidyr::unite(trt, heat.f, water.f, heat.g, water.g, heat.s, water.s) %>% 
  tidyr::nest(lambda = lambda) %>% 
  tidyr::pivot_wider(names_from = trt, values_from = lambda) %>%
  tidyr::unnest(cols = c(`1_0_0_0_0_0`, `0_1_0_0_0_0`, `1_1_0_0_0_0`, `0_0_1_0_0_0`, 
                         `0_0_0_1_0_0`, `0_0_1_1_0_0`, `0_0_0_0_1_0`, `0_0_0_0_0_1`, 
                         `0_0_0_0_1_1`, `0_0_0_0_0_0`), names_sep = "_") %>% 
  dplyr::rename(fecundity_heat =`1_0_0_0_0_0_lambda`, 
                fecundity_water =`0_1_0_0_0_0_lambda`, 
                fecundity_hw = `1_1_0_0_0_0_lambda`, 
                growth_heat = `0_0_1_0_0_0_lambda`, 
                growth_water = `0_0_0_1_0_0_lambda`, 
                growth_hw = `0_0_1_1_0_0_lambda`, 
                survival_heat = `0_0_0_0_1_0_lambda`, 
                survival_water = `0_0_0_0_0_1_lambda`, 
                survival_hw =`0_0_0_0_1_1_lambda`, 
                ambient = `0_0_0_0_0_0_lambda`) %>%
  dplyr::mutate(fecundity_heat_effect = (fecundity_heat+fecundity_hw)-(ambient+fecundity_water), 
                fecundity_water_effect = (fecundity_water+fecundity_hw) - (ambient + fecundity_heat), 
                fecundity_hw_effect = fecundity_hw - (1/3)*(ambient + fecundity_heat + fecundity_water),
                
                growth_heat_effect = (growth_heat+growth_hw)-(ambient+growth_water), 
                growth_water_effect = (growth_water+growth_hw) - (ambient + growth_heat), 
                growth_hw_effect = growth_hw - (1/3)*(ambient + growth_heat + growth_water),
                
                survival_heat_effect = (survival_heat+survival_hw)-(ambient+survival_water), 
                survival_water_effect = (survival_water+survival_hw) - (ambient + survival_heat), 
                survival_hw_effect = survival_hw - (1/3)*(ambient + survival_heat + survival_water)) %>% 
  
  dplyr::select(-fecundity_heat, -fecundity_water, -fecundity_hw, -growth_heat, -growth_water, -growth_hw, -survival_heat, -survival_water, -survival_hw, -ambient) %>% 
  tidyr::pivot_longer(names_to = "trt", values_to = "delta_lambda", -(degree.days:vwc)) 

vr_contrasts

vr_contrasts_summary <-
  vr_contrasts %>% 
  group_by(degree.days, vwc, trt) %>% 
  tidybayes::median_hdci(delta_lambda) %>% 
  dplyr::mutate(sig_delta_lambda = ifelse(test = .lower >= 0 | .upper <= 0,
                                          yes = delta_lambda,
                                          no = NA)) %>% 
  dplyr::mutate(manipulated_vr = dplyr::case_when(stringr::str_detect(string = trt, pattern = "fecundity") ~ "FECUNDITY",
                                                  stringr::str_detect(string = trt, pattern = "growth") ~ "GROWTH",
                                                  stringr::str_detect(string = trt, pattern = "survival") ~ "SURVIVAL"),
                trt = dplyr::case_when(stringr::str_detect(string = trt, pattern = "heat") ~ "HEAT",
                                       stringr::str_detect(string = trt, pattern = "water") ~ "WATER",
                                       stringr::str_detect(string = trt, pattern = "hw") ~ "HEAT+WATER"))

vr_contrasts_summary

fig5 <- 
  ggplot(data = vr_contrasts_summary, aes(x = degree.days, y= vwc, fill = sig_delta_lambda)) +
  geom_raster() +
  facet_grid(manipulated_vr ~ trt) +
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80") +
  theme_classic() +
  guides(alpha = FALSE) +
  theme(text = element_text(size=12)) +
  theme(strip.text.y = element_text(angle = 90))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  # ggtitle("HEAT")+
  geom_sf(data = ambient_stable_outline, inherit.aes = FALSE, fill = NA)

#ggsave(plot = fig5, filename = "figs/fig5-experimental-lambda-contrasts-by-vital-rate.png")

# Add in site specific contrats
site_metadata <- readr::read_csv(file = "data/data_output/site_metadata.csv")
site_lambdas <- 
  readr::read_csv(file.path("data/data_output/site_specific_additive.csv")) %>% 
  dplyr::left_join(site_metadata, by = "site") %>% 
  dplyr::mutate(elevation = round(elevation_raw, digits = 0)) %>% 
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient"))

trt_on_fecundity <-
  site_lambdas %>%
  dplyr::filter(manipulated_vr %in% c("fecundity", "ambient")) %>%
  dplyr::select(site, elevation, manipulated_vr, degree.days, vwc, heat.f, water.f, lambda)

fecundity_contrasts <-
  trt_on_fecundity %>%
  dplyr::select(-manipulated_vr) %>%
  dplyr::group_by(site, elevation) %>%
  dplyr::group_split() %>%
  lapply(FUN = function(x) {
    x <-
      x %>%
      dplyr::rename(heat = starts_with("heat"),
                    water = starts_with("water"))

    my_df <-
      tibble::tibble(site = unique(x$site),
                     elevation = unique(x$elevation),
                     degree.days = unique(x$degree.days),
                     vwc = unique(x$vwc),
                     h = x$lambda[x$heat == 1 & x$water == 0],
                     w = x$lambda[x$heat == 0 & x$water == 1],
                     hw = x$lambda[x$heat == 1 & x$water == 1],
                     a = x$lambda[x$heat == 0 & x$water == 0],
                     heat_effect = ((h + hw) / 2) - ((w + a) / 2),
                     water_effect = ((w + hw) / 2) - ((h + a) / 2),
                     hw_effect = (hw) - ((h + w + a) / 3)) %>%
      dplyr::select(site, elevation, degree.days, vwc, dplyr::ends_with("effect")) %>%
      tidyr::pivot_longer(cols = c(heat_effect, water_effect, hw_effect),
                          names_to = "trt",
                          values_to = "delta_lambda") %>%
      dplyr::group_by(site, elevation, degree.days, vwc, trt) %>%
      #dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
       #                lwr95 = quantile(delta_lambda, probs = 0.025),
        #               upr95 = quantile(delta_lambda, probs = 0.975),
         #              .groups = "keep")
  tidybayes::median_hdci(delta_lambda) %>% 
  dplyr::mutate(sig_lambda = ifelse(test = .lower >= 0 | .upper <= 0,
                                    yes = delta_lambda,
                                    no = NA))
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(manipulated_vr = "fecundity") %>%
  dplyr::select(site, elevation, degree.days, vwc, manipulated_vr, tidyr::everything())

trt_on_survivorship <-
  site_lambdas %>%
  dplyr::filter(manipulated_vr %in% c("survivorship", "ambient")) %>%
  dplyr::select(site, elevation, manipulated_vr, degree.days, vwc, heat.s, water.s, lambda)

survivorship_contrasts <-
  trt_on_survivorship %>%
  dplyr::select(-manipulated_vr) %>%
  dplyr::group_by(site, elevation) %>%
  dplyr::group_split() %>%
  lapply(FUN = function(x) {
    x <-
      x %>%
      dplyr::rename(heat = starts_with("heat"),
                    water = starts_with("water"))

    my_df <-
      tibble::tibble(site = unique(x$site),
                     elevation = unique(x$elevation),
                     degree.days = unique(x$degree.days),
                     vwc = unique(x$vwc),
                     h = x$lambda[x$heat == 1 & x$water == 0],
                     w = x$lambda[x$heat == 0 & x$water == 1],
                     hw = x$lambda[x$heat == 1 & x$water == 1],
                     a = x$lambda[x$heat == 0 & x$water == 0],
                     heat_effect = ((h + hw) / 2) - ((w + a) / 2),
                     water_effect = ((w + hw) / 2) - ((h + a) / 2),
                     hw_effect = (hw) - ((h + w + a) / 3)) %>%
      dplyr::select(site, elevation, degree.days, vwc, dplyr::ends_with("effect")) %>%
      tidyr::pivot_longer(cols = c(heat_effect, water_effect, hw_effect),
                          names_to = "trt",
                          values_to = "delta_lambda") %>%
      dplyr::group_by(site, elevation, degree.days, vwc, trt) %>%
      #dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
       #                lwr95 = quantile(delta_lambda, probs = 0.025),
        #               upr95 = quantile(delta_lambda, probs = 0.975),
         #              .groups = "keep")
  tidybayes::median_hdci(delta_lambda) %>% 
  dplyr::mutate(sig_lambda = ifelse(test = .lower >= 0 | .upper <= 0,
                                    yes = delta_lambda,
                                    no = NA))
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(manipulated_vr = "survivorship") %>%
  dplyr::select(site, elevation, degree.days, vwc, manipulated_vr, tidyr::everything())

trt_on_growth <-
  site_lambdas %>%
  dplyr::filter(manipulated_vr %in% c("growth", "ambient")) %>%
  dplyr::select(site, elevation, manipulated_vr, degree.days, vwc, heat.g, water.g, lambda)

growth_contrasts <-
  trt_on_growth %>%
  dplyr::select(-manipulated_vr) %>%
  dplyr::group_by(site, elevation) %>%
  dplyr::group_split() %>%
  lapply(FUN = function(x) {
    x <-
      x %>%
      dplyr::rename(heat = starts_with("heat"),
                    water = starts_with("water"))

    my_df <-
      tibble::tibble(site = unique(x$site),
                     elevation = unique(x$elevation),
                     degree.days = unique(x$degree.days),
                     vwc = unique(x$vwc),
                     h = x$lambda[x$heat == 1 & x$water == 0],
                     w = x$lambda[x$heat == 0 & x$water == 1],
                     hw = x$lambda[x$heat == 1 & x$water == 1],
                     a = x$lambda[x$heat == 0 & x$water == 0],
                     heat_effect = ((h + hw) / 2) - ((w + a) / 2),
                     water_effect = ((w + hw) / 2) - ((h + a) / 2),
                     hw_effect = (hw) - ((h + w + a) / 3)) %>%
      dplyr::select(site, elevation, degree.days, vwc, dplyr::ends_with("effect")) %>%
      tidyr::pivot_longer(cols = c(heat_effect, water_effect, hw_effect),
                          names_to = "trt",
                          values_to = "delta_lambda") %>%
      dplyr::group_by(site, elevation, degree.days, vwc, trt) %>%
      #dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
       #                lwr95 = quantile(delta_lambda, probs = 0.025),
        #               upr95 = quantile(delta_lambda, probs = 0.975),
         #              .groups = "keep")
  tidybayes::median_hdci(delta_lambda) %>% 
  dplyr::mutate(sig_lambda = ifelse(test = .lower >= 0 | .upper <= 0,
                                    yes = delta_lambda,
                                    no = NA))
    
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(manipulated_vr = "growth") %>%
  dplyr::select(site, elevation, degree.days, vwc, manipulated_vr, tidyr::everything())

all_contrasts <-
  rbind(fecundity_contrasts, survivorship_contrasts, growth_contrasts) %>% 
  dplyr::mutate(
    manipulated_vr = dplyr::case_when(
      stringr::str_detect(string = manipulated_vr, pattern = "fecundity") ~ "FECUNDITY",
      stringr::str_detect(string = manipulated_vr, pattern = "growth") ~ "GROWTH",
      stringr::str_detect(string = manipulated_vr, pattern = "survivorship") ~ "SURVIVAL"),
      trt = dplyr::case_when(stringr::str_detect(string = trt, pattern = "heat") ~ "HEAT",
      stringr::str_detect(string = trt, pattern = "water") ~ "WATER",
      stringr::str_detect(string = trt, pattern = "hw") ~ "HEAT+WATER"))


# NEW FIGURE 5
fig5 <- 
  ggplot(data = vr_contrasts_summary, aes(x = degree.days, y= vwc, fill = sig_delta_lambda)) +
  geom_raster() +
  facet_grid(manipulated_vr ~ trt) +
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80") +
  theme_classic() +
  guides(alpha = FALSE) +
  theme(text = element_text(size=12)) +
  theme(strip.text.y = element_text(angle = 90))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  geom_sf(data = ambient_stable_outline, inherit.aes = FALSE, fill = NA)+
  
geom_point(data = all_contrasts[is.na(all_contrasts$sig_lambda),], aes(degree.days, y= vwc, bg = sig_lambda), pch = 21,show.legend = FALSE, cex = 3)+

geom_point(data = all_contrasts[which(all_contrasts$sig_lambda < 0),], aes(degree.days, y= vwc, bg = sig_lambda), pch = 25,show.legend = FALSE, cex = 3)+

geom_point(data = all_contrasts[which(all_contrasts$sig_lambda > 0 ),], aes(degree.days, y= vwc, bg = sig_lambda), pch = 24,show.legend = FALSE, cex = 3)

fig5 

ggsave(plot = fig5, filename = "figs/fig5-experimental-lambda-contrasts-by-vital-rate-with -sites.png")




