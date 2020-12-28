#Analyze site-specific IPM output and make a corresponding figure

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

# # site-specific lambdas per elevation
# ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("ambient")), aes(x = elevation, y = lambda)) +
#   tidybayes::stat_halfeye(limits = c(0,2))+
#   theme_classic()
# 
# ggplot(site_lambdas %>% dplyr::filter(manipulated_vr %in% c("heat", "water", "ambient", "hw")), aes(x = elevation, y = lambda, col = manipulated_vr)) +
#   tidybayes::stat_halfeye() +
#   facet_wrap(facets = "manipulated_vr")+
#   theme_classic()

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
  tidyr::pivot_longer(names_to = "trt", values_to = "delta_lambda", -(site:elevation)) %>% 
  dplyr::select(site, elevation, vwc, degree.days, delta_lambda, site.num, everything())

contrasts

contrasts_summary <-
  contrasts %>% 
  dplyr::group_by(site, elevation, degree.days, vwc, trt) %>% 
  tidybayes::median_hdci(delta_lambda) %>% 
  dplyr::mutate(sig_delta_lambda = ifelse(test = .lower >= 0 | .upper <= 0,
                                          yes = delta_lambda,
                                          no = NA))

fig3 <-
  ggplot(contrasts_summary, aes(x = elevation, y = delta_lambda)) +
  geom_point() +
  geom_errorbar(aes(ymin = .lower, ymax = .upper)) +
  facet_wrap(.~trt, ncol = 1, labeller = labeller(trt = c(heat_effect  = "HEAT", water_effect= "WATER", hw_effect = "HEAT + WATER")))+
  theme_classic() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90))+
  xlab("Elevation")+
  ylab(expression(paste(Delta, lambda)))+
  geom_hline(yintercept = 0, color = "grey")

fig3

ggsave(filename = "figs/fig3_site-specific-lambda-contrasts-experimental-treatments.png", plot = fig4, dpi = 600)

#### Contrasts of simulated treatment on individual vital rates against ambientt
# fecundity ---------------------------------------------------------------

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
      dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
                       lwr95 = quantile(delta_lambda, probs = 0.025),
                       upr95 = quantile(delta_lambda, probs = 0.975),
                       .groups = "keep")
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(manipulated_vr = "fecundity") %>%
  dplyr::select(site, elevation, degree.days, vwc, manipulated_vr, tidyr::everything())


ggplot(fecundity_contrasts, aes(x = elevation, y = mean_delta_lambda)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr95, ymax = upr95)) +
  facet_wrap(.~trt, nrow = 1, labeller = labeller(trt = c(heat_effect  = "HEAT", water_effect= "WATER", hw_effect = "HEAT + WATER")))+
  theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90))+
  xlab("Elevation")+
  ylab(expression(paste(Delta, lambda)))+
  geom_hline(yintercept = 0, color = "grey")
#
#
# # survivorship --------------------------------------------------------
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
      dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
                       lwr95 = quantile(delta_lambda, probs = 0.025),
                       upr95 = quantile(delta_lambda, probs = 0.975),
                       .groups = "keep")
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(manipulated_vr = "survivorship") %>%
  dplyr::select(site, elevation, degree.days, vwc, manipulated_vr, tidyr::everything())


ggplot(survivorship_contrasts, aes(x = elevation, y = mean_delta_lambda)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr95, ymax = upr95)) +
  facet_wrap(.~trt, nrow = 1, labeller = labeller(trt = c(heat_effect  = "HEAT", water_effect= "WATER", hw_effect = "HEAT + WATER")))+
  theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90))+
  xlab("Elevation")+
  ylab(expression(paste(Delta, lambda)))+
  geom_hline(yintercept = 0, color = "grey")
#
# # growth -------------------------------------------------------------

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
      dplyr::summarize(mean_delta_lambda = mean(delta_lambda),
                       lwr95 = quantile(delta_lambda, probs = 0.025),
                       upr95 = quantile(delta_lambda, probs = 0.975),
                       .groups = "keep")
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(manipulated_vr = "growth") %>%
  dplyr::select(site, elevation, degree.days, vwc, manipulated_vr, tidyr::everything())


ggplot(growth_contrasts, aes(x = elevation, y = mean_delta_lambda)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr95, ymax = upr95)) +
  facet_wrap(.~trt, nrow = 1, labeller = labeller(trt = c(heat_effect  = "HEAT", water_effect= "WATER", hw_effect = "HEAT + WATER")))+
  theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90))+
  xlab("Elevation")+
  ylab(expression(paste(Delta, lambda)))+
  geom_hline(yintercept = 0, color = "grey")
#
#
# # all contrasts -------------------------------------------------------
#
all_contrasts <-
  rbind(fecundity_contrasts, survivorship_contrasts, growth_contrasts)

fig4b <-
  ggplot(all_contrasts, aes(x = elevation, y = mean_delta_lambda)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr95, ymax = upr95)) +
  facet_grid(rows = vars(manipulated_vr), cols = vars(trt), labeller = labeller(trt = c(heat_effect  = "HEAT", water_effect= "WATER", hw_effect = "HEAT + WATER")))+
  theme_bw() +
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90))+
  xlab("Elevation")+
  ylab(expression(paste(Delta, lambda)))+
  geom_hline(yintercept = 0, color = "grey")+
  ggtitle("B")
# 
ggsave(filename = "figs/fig3b_site-specific-lambda-contrasts-decomposed-by-vital-rate.png", plot = fig4b, dpi = 600)
