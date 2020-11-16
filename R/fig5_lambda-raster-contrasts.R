# Figure 5 showing 

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

mcvr_lambda_summary$trt <- "ambient"
mcvr_lambda_summary[mcvr_lambda_summary$heat == 1 & mcvr_lambda_summary$water == 0, "trt"] <- "heat"
mcvr_lambda_summary[mcvr_lambda_summary$heat == 1 & mcvr_lambda_summary$water == 1, "trt"] <- "hw"
mcvr_lambda_summary[mcvr_lambda_summary$heat == 0 & mcvr_lambda_summary$water == 1, "trt"] <- "water"

Fig5A <- mcvr_lambda_summary %>% 
  filter(heat == 1 | water == 1) %>% 
  ggplot(aes(x = degree.days, y = vwc, fill = sig_lambda)) +
  geom_raster()+
  scale_fill_gradient2(mid = "grey90", midpoint = 1, na.value = "grey90", high = "darkgreen", low = "darkgoldenrod1")+
  theme_classic() +
  guides(alpha = F)+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  theme(text = element_text(size=16))+
  labs(fill = expression(lambda))+
  facet_grid(.~trt, labeller = labeller(trt = c(heat = "HEAT", water = "WATER", hw = "HEAT + WATER")))+
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA)

Fig5A


#panel B
contrasts <- 
  mcvr %>% 
  filter(manipulated_vr %in% c("ambient", "heat", "water", "hw")) %>% 
  mutate(heat = ifelse(heat.f == 1, yes = 1, no = 0), 
         water = ifelse(water.f == 1, yes = 1, no = 0)) %>% 
  dplyr::select(-heat.f:-water.s, - manipulated_vr) %>% 
  unite(trt, heat, water) %>% 
  nest(lambda = lambda) %>% 
  pivot_wider(names_from = trt, values_from = lambda) %>%
  unnest(cols = c(`0_0`, `1_0`, `0_1`, `1_1`),names_sep = "_") %>% 
  rename(ambient = `0_0_lambda`, heat = `1_0_lambda`, water = `0_1_lambda`, hw = `1_1_lambda`) %>%
  mutate(heat_effect = (heat+hw)-(ambient+water), water_effect = (water+hw) - (ambient + heat), hw_effect = hw - (1/3)*(ambient + heat + water)) %>% 
  dplyr::select(-heat, -water, -ambient, -hw) %>% 
  tidyr::pivot_longer(names_to = "trt", values_to = "delta_lambda", -(degree.days:vwc)) 

contrasts

contrasts_summary <-
  contrasts %>% 
  group_by(degree.days, vwc, trt) %>% 
  summarize(mean_delta_lambda = mean(delta_lambda),
            median_delta_lambda = median(delta_lambda),
            lwr95 = quantile(delta_lambda, probs = 0.025),
            upr95 = quantile(delta_lambda, probs = 0.975),
            mass_over_under_zero = ifelse(mean_delta_lambda < 0,
                                          yes = ecdf(delta_lambda)(0),
                                          no = 1 - ecdf(delta_lambda)(0)),
            sig_lambda = ifelse(mass_over_under_zero < 0.95,
                                yes = NA, no = mean_delta_lambda), 
            n=n())

head(contrasts_summary)
range(contrasts_summary$mean_delta_lambda)

Fig5B <- 
  ggplot(contrasts_summary, aes(x = degree.days, y= vwc, fill = sig_lambda)) +
  geom_raster()+
  facet_grid(.~trt, labeller = labeller(trt = c(heat_effect = "HEAT", hw_effect = "HEAT+WATER", water_effect = "WATER"))) +
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80")+
  theme_classic() +
  guides(alpha = F)+
  theme(text = element_text(size=16))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  geom_sf(data = outline, inherit.aes = FALSE, fill = NA)
Fig5B

plot_grid(Fig5A, Fig5B, NULL, nrow = 3, labels = c("A", "B", "")) 
