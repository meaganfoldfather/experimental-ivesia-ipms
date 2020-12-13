
#### Library ####
library(dplyr)
library(ggplot2)
library(brms)
library(purrr)
library(furrr)
library(data.table)
library(glue)
library(tidyr)
library(cowplot)

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

lambda_df_fname <- "mc_vr_effects_additive.csv"

overwrite <- FALSE

##### Get vital rate models from remote_source ####

if(!file.exists(file.path("data/data_output", surv_mod_fname))) {
  download.file(url = glue::glue("{remote_source}/{surv_mod_fname}"), 
                destfile = file.path("data/data_output", surv_mod_fname))
}

if(!file.exists(file.path("data/data_output", growth_mod_fname))) {
  download.file(url = glue::glue("{remote_source}/{growth_mod_fname}"),
                destfile = file.path("data/data_output", growth_mod_fname))
}

if(!file.exists(file.path("data/data_output", establishment_mod_fname))) {
  
  download.file(url = glue::glue("{remote_source}/{establishment_mod_fname}"),
                destfile = file.path("data/data_output", establishment_mod_fname))
}

if(!file.exists(file.path("data/data_output", hurdle_mod_fname))) {
  download.file(url = glue::glue("{remote_source}/{hurdle_mod_fname}"),
                destfile = file.path("data/data_output", hurdle_mod_fname))
}

# If the files STILL don't exist on disk, then they weren't available from the
# remote source, so we have to run the script to build the vital rate models
# using the filenames declared above
if(!file.exists(file.path("data/data_output", surv_mod_fname)) |
   !file.exists(file.path("data/data_output", growth_mod_fname)) |
   !file.exists(file.path("data/data_output", establishment_mod_fname)) |
   !file.exists(file.path("data/data_output", hurdle_mod_fname))) {
  
  source("R/build-vital-rate-models.R")
}

survival_model <- readr::read_rds(file.path("data/data_output", surv_mod_fname))
growth_model <- readRDS(file.path("data/data_output", growth_mod_fname))
establishment_model <- readRDS(file.path("data/data_output", establishment_mod_fname))
hurdleRep <- readRDS(file.path("data/data_output", hurdle_mod_fname))

# Colorblind friendly palette from from https://davidmathlogic.com/colorblind/#%23FFC20A-%230C7BDC-%23eb1416-%23140b2c

four_colors <- rev(c("#FFC20A", "#0C7BDC", "#EB1416", "#140B2C"))

traceplot <- function(model, n_pars = 28, par_filter = NA_character_, negate = FALSE) {
  samps <- 
    brms::posterior_samples(model, add_chain = TRUE) %>% 
    tidyr::pivot_longer(cols = c(-iter, -chain), names_to = "parameter", values_to = "estimate") 
  if (is.na(par_filter)) {
    samps <- 
      samps %>% 
      dplyr::filter(parameter %in% unique(parameter)[1:n_pars]) %>% 
      dplyr::mutate(parameter = factor(parameter, levels = unique(parameter)))
  } else {
    samps <-
      samps %>% 
      dplyr::filter(stringr::str_detect(string = parameter, pattern = {{ par_filter }}, negate = negate)) %>% 
      dplyr::filter(parameter %in% unique(parameter)[1:n_pars]) %>% 
      dplyr::mutate(parameter = factor(parameter, levels = unique(parameter)))
  }
  
  the_traceplot <- 
    ggplot(samps, aes(x = iter, y = estimate, color = chain)) +
    geom_line(alpha = 0.25) +
    facet_wrap(facets = vars(parameter), ncol = 4, scales = "free_y") +
    theme_bw() +
    labs(x = "iteration") +
    scale_color_manual(values = four_colors)
  
}
#### Survivorship model trace plot
survival_traceplot <- traceplot(survival_model)
ggsave(survival_traceplot, filename = "figs/supp-figXX-survival-model-traceplot.pdf", width = 10, height = 12)

#### Growth model traceplot
growth_traceplot <- traceplot(growth_model)
ggsave(growth_traceplot, filename = "figs/supp-figXX-growth-model-traceplot.pdf", width = 10, height = 12)

#### Establishment model
establishment_traceplot <- traceplot(establishment_model, n_pars = 15)
ggsave(establishment_traceplot, filename = "figs/supp-figXX-establishment-model-traceplot.pdf", width = 10, height = 8)

### recruitment model
hurdleRep_traceplot_hu <- traceplot(hurdleRep, par_filter = "hu", n_pars = 26)
ggsave(hurdleRep_traceplot_hu, filename = "figs/supp-figXX-hurdleRep-model-hu-traceplot.pdf", width = 10, height = 12)

hurdleRep_traceplot_not_hu <- traceplot(hurdleRep, par_filter = "hu", n_pars = 28, negate = TRUE)
ggsave(hurdleRep_traceplot_not_hu, filename = "figs/supp-figXX-hurdleRep-model-not-hu-traceplot.pdf", width = 10, height = 12)

