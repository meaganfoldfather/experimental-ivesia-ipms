# Figure 5a showing the effect of experimental treatments on lambda
# Figure 5b showing the contrasts (delta lambda) between experimental treatments
# and ambient

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
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

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

### Get the site metadata
site_metadata <- readr::read_csv(file = "data/data_output/site_metadata.csv")

#panel A
contrasts <- 
  mcvr %>% 
  dplyr::filter(manipulated_vr %in% c("ambient", "heat", "water", "hw")) %>% 
  dplyr::mutate(heat = ifelse(heat.f == 1, yes = 1, no = 0), 
                water = ifelse(water.f == 1, yes = 1, no = 0)) %>% 
  dplyr::select(-heat.f:-water.s, - manipulated_vr) %>% 
  tidyr::unite(trt, heat, water) %>% 
  tidyr::nest(lambda = lambda) %>% 
  tidyr::pivot_wider(names_from = trt, values_from = lambda) %>%
  unnest(cols = c(`0_0`, `1_0`, `0_1`, `1_1`),names_sep = "_") %>% 
  rename(ambient = `0_0_lambda`, heat = `1_0_lambda`, water = `0_1_lambda`, hw = `1_1_lambda`) %>%
  mutate(heat_effect = (heat+hw)-(ambient+water), water_effect = (water+hw) - (ambient + heat), hw_effect = hw - (1/3)*(ambient + heat + water)) %>% 
  dplyr::select(-heat, -water, -ambient, -hw) %>% 
  tidyr::pivot_longer(names_to = "trt", values_to = "delta_lambda", -(degree.days:vwc)) 

contrasts

contrasts_summary <-
  contrasts %>% 
  group_by(degree.days, vwc, trt) %>% 
  tidybayes::median_hdci(delta_lambda) %>% 
  dplyr::mutate(sig_delta_lambda = ifelse(test = .lower >= 0 | .upper <= 0,
                                          yes = delta_lambda,
                                          no = NA))


head(contrasts_summary)
range(contrasts_summary$delta_lambda)

Fig4A <- 
  ggplot(contrasts_summary, aes(x = degree.days, y= vwc, fill = sig_delta_lambda)) +
  geom_raster()+
  facet_grid(.~trt, labeller = labeller(trt = c(heat_effect = "HEAT", hw_effect = "HEAT+WATER", water_effect = "WATER"))) +
  scale_fill_gradient2(midpoint = 0, mid = "grey80", na.value = "grey80")+
  theme_classic() +
  guides(alpha = F)+
  theme(text = element_text(size=16))+
  labs(fill = expression(paste(Delta, lambda)))+
  xlab("Degree-Days")+
  ylab("Soil Moisture")+
  geom_sf(data = ambient_stable_outline, inherit.aes = FALSE, fill = NA)

Fig4A

# Panel B

heat_stable_outline <- 
  make_stable_outline(data = mcvr_lambda_summary, h = 1, w = 0) %>% 
  dplyr::mutate(trt = "heat")
water_stable_outline <- 
  make_stable_outline(data = mcvr_lambda_summary, h = 0, w = 1) %>% 
  dplyr::mutate(trt = "water")
hw_stable_outline <- 
  make_stable_outline(data = mcvr_lambda_summary, h = 1, w = 1) %>% 
  dplyr::mutate(trt = "hw")

trt_stable_outline <-
  rbind(heat_stable_outline, water_stable_outline, hw_stable_outline)

mcvr_lambda_summary$trt <- "ambient"
mcvr_lambda_summary[mcvr_lambda_summary$heat == 1 & mcvr_lambda_summary$water == 0, "trt"] <- "heat"
mcvr_lambda_summary[mcvr_lambda_summary$heat == 1 & mcvr_lambda_summary$water == 1, "trt"] <- "hw"
mcvr_lambda_summary[mcvr_lambda_summary$heat == 0 & mcvr_lambda_summary$water == 1, "trt"] <- "water"

fig4b_data <-
  mcvr_lambda_summary %>% 
  filter(heat == 1 | water == 1)

Fig4B <-     
  ggplot(fig5b_data, aes(x = degree.days, y = vwc, fill = lambda)) +
  geom_raster() +
  scale_fill_gradient2(mid = "grey90", 
                       midpoint = 1, 
                       na.value = "grey90", 
                       high = "darkgreen", 
                       low = "darkgoldenrod1")+
  theme_classic() +
  guides(alpha = F) +
  xlab("Degree-Days") +
  ylab("Soil Moisture") +
  theme(text = element_text(size=16)) +
  labs(fill = expression(lambda)) +
  geom_sf(data = ambient_stable_outline, inherit.aes = FALSE, fill = NA, lty = 1) +
  ggpattern::geom_sf_pattern(data = trt_stable_outline, 
                             pattern = "circle",
                             inherit.aes = FALSE, 
                             fill = NA, 
                             lty = 1,
                             color = NA) +
  facet_grid(.~trt, labeller = labeller(trt = c(heat = "HEAT", water = "WATER", hw = "HEAT + WATER")))

Fig4B

fig4 <- plot_grid(Fig5A, Fig5B,nrow = 2, labels = c("A", "B")) 

ggsave(plot = fig4, filename = "figs/fig4-experimental-lambda-contrasts-and-lambda-across-microclimates.png")
