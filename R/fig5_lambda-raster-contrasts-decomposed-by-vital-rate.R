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
  summarize(mean_delta_lambda = mean(delta_lambda),
            median_delta_lambda = median(delta_lambda),
            lwr95 = quantile(delta_lambda, probs = 0.025),
            upr95 = quantile(delta_lambda, probs = 0.975),
            mass_over_under_zero = ifelse(mean_delta_lambda < 0,
                                          yes = ecdf(delta_lambda)(0),
                                          no = 1 - ecdf(delta_lambda)(0)),
            sig_lambda = ifelse(mass_over_under_zero < 0.95,
                                yes = NA, no = mean_delta_lambda))
vr_contrasts_summary

# heat plot
vr_contrasts_summary_heat <- vr_contrasts_summary[grepl(vr_contrasts_summary$trt, pattern = "_heat_"),]

heat_plot1 <- 
  ggplot(data = vr_contrasts_summary_heat, aes(x = degree.days, y= vwc, fill = sig_lambda)) +
  geom_raster()+
  facet_grid(trt~., labeller = labeller(trt = c(survival_heat_effect = "SURVIVAL",growth_heat_effect = "GROWTH", fecundity_heat_effect = "FECUNDITY")))+
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80",limits = c(-.4, .1))+
  theme_classic() +
  guides(alpha = F)+
  theme(text = element_text(size=12))+
  theme(strip.text.y = element_text(angle = 90))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  ggtitle("HEAT")+
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA)

heat_plot1

#water plots
vr_contrasts_summary_water <- vr_contrasts_summary[grepl(vr_contrasts_summary$trt, pattern = "_water_"),]

water_plot1 <- 
  ggplot(data = vr_contrasts_summary_water, aes(x = degree.days, y= vwc, fill = sig_lambda)) +
  geom_raster()+
  facet_grid(trt~., labeller = labeller(trt = c(survival_water_effect = "SURVIVAL",growth_water_effect = "GROWTH", fecundity_water_effect = "FECUNDITY")))+
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80",limits = c(-.4, .1))+
  theme_classic() +
  guides(alpha = F)+
  theme(text = element_text(size=12))+
  theme(strip.text.y = element_text(angle = 90))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  ggtitle("WATER")+
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA)

water_plot1

# hw plot
vr_contrasts_summary_hw <- vr_contrasts_summary[grepl(vr_contrasts_summary$trt, pattern = "_hw_"),]

hw_plot1 <- 
  ggplot(data = vr_contrasts_summary_hw, aes(x = degree.days, y= vwc, fill = sig_lambda)) +
  geom_raster()+
  facet_grid(trt~., labeller = labeller(trt = c(survival_hw_effect = "SURVIVAL",growth_hw_effect = "GROWTH", fecundity_hw_effect = "FECUNDITY")))+
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80",  limits = c(-.4, .1))+
  
  theme_classic() +
  guides(alpha = F)+
  theme(text = element_text(size=12))+
  theme(strip.text.y = element_text(angle = 90))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  ggtitle("HEAT+WATER")+
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA)

hw_plot1

### all plot
#combine plots
prow <- 
  plot_grid(
    heat_plot1 + theme(legend.position="none"),
    hw_plot1 + theme(legend.position="none"),
    water_plot1 + theme(legend.position="none"), nrow = 1, labels = c("C.i", "C.ii", "C .iii"))

legend_b <- get_legend(heat_plot1)

plot_grid(prow, legend_b, ncol = 2, rel_widths =c(1, .2))
