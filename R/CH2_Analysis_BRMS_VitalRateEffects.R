#Intent: Analysis for CH2 - Analysis to look at Specific Responses of Each Vital Rate to the Experimental Manipulations
#Date Last Edited: 20190915

#### Library ####
library(dplyr)
library(brms)
library(purrr)
library(furrr)
library(data.table)
library(glue)

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
   !file.exists(file.path("data/data_output", hurdle_mod_fname)) |
   !file.exists(file.path("data/data_output", hurdle_mod_fname))) {
  
  source("R/CH2_VRModels_Revised.R")
}

survival_model <- readr::read_rds(file.path("data/data_output", surv_mod_fname))
growth_model <- readRDS(file.path("data/data_output", growth_mod_fname))
establishment_model <- readRDS(file.path("data/data_output", establishment_mod_fname))
hurdleRep <- readRDS(file.path("data/data_output", hurdle_mod_fname))

# R squared of vital rate models
bayes_R2(survival_model) 
bayes_R2(growth_model)
bayes_R2(hurdleRep)
bayes_R2(establishment_model)

plot1 <- pp_check(survival_model, nsamples = 100); plot1
plot2 <- pp_check(growth_model, nsamples = 100); plot2
plot3 <- pp_check(hurdleRep, nsamples = 100); plot3
plot4 <- pp_check(establishment_model, nsamples = 100); plot4

library(cowplot)
plot_grid(plot1, plot2, plot3, plot4, labels = c("A", "B", "C", "D"))

#### Bring in vital rate data --> I.df ####
vr.mc<-read.csv(glue::glue("{remote_source}/VitalRates_Microclimate.csv"), stringsAsFactors = F); head(vr.mc)

#Add in experimental metdata
plot.chars<-read.csv(glue::glue("{remote_source}/trt.plot.data.csv"))
colnames(plot.chars)<-c("plot","zone","site","trt","type", "elevation");plot.chars$plot<-gsub("IL","",plot.chars$plot)
head(plot.chars)

vr.mc<-merge(vr.mc, plot.chars, all.x=T, by=c("plot","site")); head(vr.mc)

# Add in delta column for difference between size and sizeNext
vr.mc$delta<-vr.mc$sizeNext - vr.mc$size

# Add in Heat and Water columns
vr.mc$heat<-0; vr.mc$water<-0; vr.mc[vr.mc$trt == "H" | vr.mc$tr == "HW","heat"]<-1; vr.mc[vr.mc$trt == "W" | vr.mc$trt == "HW","water"]<-1
head(vr.mc); tail(vr.mc); dim(vr.mc)

# Change df to I.df 
I.df<-vr.mc; head(I.df); dim(I.df)

# Add column that is probability of flowering & subsequently make the flwr magnitude zero --> NA
I.df$fprob<-NA
I.df$fprob[which(I.df$fec>0)] <- 1
I.df$fprob[which(I.df$fec==0)] <- 0
I.df$fec[which(I.df$fec==0)]<- NA

I.df$fprobNext<-NA 
I.df$fprobNext[which(I.df$fecNext>0)] <- 1
I.df$fprobNext[which(I.df$fecNext==0)] <- 0
I.df$fecNext[which(I.df$fecNext==0)]<- NA

# Add a column that the number of seeds produced in each year
I.df$seed<-I.df$fec*5
I.df$seedNext<-I.df$fecNext*5

#change t to a factor
I.df$t<-as.factor(I.df$t)
I.df$t1<-as.factor(I.df$t1)

# subset to only experimental sites
I.df<-subset(I.df, type == "EXPERIMENTAL")

# final data df 
I.df[1:20,]; dim (I.df) # 8785 individuals, 9 sites 


#### Scales size and included ambient climate vectors --> df ####
# create  dataframe called df
df <- I.df
head(df)

# get rid of current climate data 
df <- df [, -c(12:17)]

# bring in climate data averaged across all ambient plots and years
microclimate<-read.csv(glue::glue("{remote_source}/Microclimate.csv"))
head(microclimate)
head(plot.chars)
mc.trt<-merge(microclimate, plot.chars)
head(mc.trt)
mc<-mc.trt %>%
  filter(type == "EXPERIMENTAL") %>% 
  filter(trt == "A") %>%
  group_by(site) %>% 
  summarise(vwc = round(mean(vwc, na.rm=T), 0), degree.days = round(mean(degree.days, na.rm=T),0), snow.days = round(mean(snow.days, na.rm=T,0)))
mc

df <- merge(df, mc, by = c("site"), all.x =T)
head(df); dim(df)

#scaling size and sizeNext by the sizes of all individuals measured across all timepoints, sites and treatments
df$site <- as.factor(df$site)
df$plot <- as.factor(df$plot)

#filter out rows with 0 for size or size next
df <- df[-which(df$size == 0),]
df <- df[-which(df$sizeNext == 0),]

size_vec <- c(df$size, df$sizeNext)
df$size.s <- (df$size - mean(size_vec, na.rm=T))/sd(size_vec, na.rm=T)
df$sizeNext.s <- (df$sizeNext - mean(size_vec, na.rm=T))/sd(size_vec, na.rm=T)

# scale climate vectors
df$vwc<- scale(df$vwc)
df$degree.days <- scale(df$degree.days)
df$snow.days<- scale(df$snow.days)

# scale elevation vectors
df$elevation.s <- scale(df$elevation)

# kernel construction
min.size=min(df[,c("size.s", "sizeNext.s")], na.rm=T)
max.size=max(df[, c("size.s", "sizeNext.s")], na.rm=T)
n=100 # number of cells in the matrix
b=min.size+c(0:n)*(max.size-min.size)/n # boundary points
y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points
h=y[2]-y[1] # step size

#1. survival probability function
s.x=function(x, degree.days = 0, vwc = 0, heat = 0, water = 0) {
  new.data <- data.frame(size.s = x, degree.days = degree.days, vwc = vwc, heat = heat, water = water, t1 = NA) # use t1 = NA so that the fitted() function marginalizes across all levels of the t1 categorical predictor, using the grand mean (because we use sum coding of that factor)
  # https://github.com/paul-buerkner/brms/issues/66
  survivorship <- fitted(survival_model, newdata = new.data, re_formula = NA, scale = "response", summary = FALSE)
  return(survivorship)
}

# 2. growth function
g.yx=function(xp,x,degree.days = 0, vwc = 0, heat = 0, water = 0) {
  new.data <- data.frame(size.s = x, degree.days = degree.days, vwc = vwc, heat = heat, water = water, t1 = NA)
  sizeNext.s <- predict(growth_model, newdata = new.data, re_formula = NA, scale = "response", summary = FALSE) 
  
  
  growth <- apply(sizeNext.s, MARGIN = 2, FUN = function(sizeNext.s_column) {approxfun(density(sizeNext.s_column))(xp)})
  # if the probability is so small that it doesn't exist in our empirical
  # probability distribution (and thus the approxfun(density(x)) function 
  # returns NA), just make the probability 0
  growth[is.na(growth)] <- 0 
  
  return(growth)
}

# 3. reproduction function
f.yx <- function(x, xp ,degree.days = 0, vwc = 0, heat = 0, water = 0) {
  new.data <- data.frame(sizeNext.s = x, degree.days = degree.days, vwc = vwc, heat = heat, water = water, t1 = NA)
  
  reproduction <- fitted(hurdleRep, newdata = new.data, re_formula = NA, scale = "response", summary = FALSE)
  reproduction
  dim(reproduction)
  
  establishment_probability <- fitted(establishment_model, newdata = data.frame(new.data, seeds.avaliable = 1), re_formula = NA, scale = "response", summary = T)[1, "Estimate"]
  
  recruit.size.probabilty <- approxfun(density(df[which(df$new == 1), "sizeNext.s"], na.rm = T))(xp)
  # code is based on https://stats.stackexchange.com/questions/78711/how-to-find-estimate-probability-density-function-from-density-function-in-r
  # if the probability is so small that it doesn't exist in our empirical
  # probability distribution (and thus the approxfun(density(x)) function 
  # returns NA), just make the probability 0
  recruit.size.probabilty[is.na(recruit.size.probabilty)] <- 0 
  length(recruit.size.probabilty) #100
  recruit.size.probabilty
  
  fecundity <- (reproduction * establishment_probability) # expected number of new recruits for a individual of a certain size
  
  # Our goal is, for each size, determine the transition rate to *all* possible sizes attributed to fecundity. This would look like a 100 x 100 matrix (with each size_t as a different column and each size_(t+1) as a different row). Because we have a bunch of samples from our model, we want a 3 dimensional array such that one of the above-mentioned 100 x 100 matrices will be created for each sample (row) in the original fecundity matrix (which is 100 x 2000 [num of size classes x num of samples])
  fecundity_3d <- lapply(data.frame(t(fecundity)), FUN = function(f) recruit.size.probabilty %o% f)
  
  return(fecundity_3d)
}

#### Making Whole Construction Function ####
kernel_construction <- function(x, degree.days = 0, vwc = 0, heat.s = 0, water.s = 0, heat.g = 0, water.g = 0, heat.f = 0, water.f = 0){
  
  F.mat = h*simplify2array(f.yx(x, x, degree.days, vwc, heat = heat.f, water = water.f))
  S = s.x(x, degree.days, vwc, heat = heat.s, water = water.s)
  S = array(t(S), dim = c(1, 100, 2000))
  G = h * g.yx(x, x, degree.days, vwc, heat = heat.g, water = water.g)
  G = array(G, dim = c(100, 100, 1))
  
  P = array(dim=c(100,100,2000))
  for(k in 1:n) P[,k,]=G[,k,]*S[,k,] # growth/survival matrix
  # full matrix
  K=P+F.mat
  
  lambda <- apply(K, MARGIN = 3, FUN = function(k) Re(eigen(k)$values[1]))
  return(lambda)
}

#Treatment Effects Across Climatic Gradients
# Running on the Alien with 10 workers on 2020-05-06; runtime = 

plan(multiprocess, workers = 10)
start_time <- Sys.time()
start_time

vital_effects <- 
  data.frame(
    heat.f = c(1,0,1, rep(0,7)),
    water.f = c(0,1,1, rep(0,7)), 
    heat.g = c(0,0,0,1,0,1,0,0,0,0), 
    water.g = c(0,0,0,0,1,1,0,0,0,0), 
    heat.s = c(rep(0,6),1,0,1,0), 
    water.s = c(rep(0,7),1,1,0)) %>% 
  rbind(c(1, 0, 1, 0, 1, 0)) %>% # heat every vital rate
  rbind(c(0, 1, 0, 1, 0, 1)) %>% # water every vital rate
  rbind(c(1, 1, 1, 1, 1, 1)) # heat and water every vital rate

microclimate_effects <-
  tidyr::crossing(degree.days = seq(-1.5, 2, by = 0.5),
                  vwc = seq(-1.5, 1, by = 0.5))
mc_vr_effects <-
  microclimate_effects %>% 
  dplyr::mutate(vital_effects = list(vital_effects)) %>% 
  tidyr::unnest(vital_effects) %>% 
  mutate(lambda = furrr::future_pmap(., kernel_construction, y)) %>%
  tidyr::unnest()

mcvr <- 
  mc_vr_effects %>% 
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient"))

if(overwrite | !file.exists(glue::glue("data/data_output/{lambda_df_fname}"))) {
  data.table::fwrite(x = mcvr, file = here::here("data", "data_output", lambda_df_fname))
  
  system2(command = "aws", args = glue::glue("s3 cp data/data_output/{lambda_df_fname} {remote_target}/{lambda_df_fname}"))
}
(end_time <- Sys.time())

(base::difftime(end_time, start_time, units = "mins"))


if(!file.exists(glue::glue("data/data_output/{lambda_df_fname}"))) {
  
  download.file(url = glue::glue("{remote_source}/{lambda_df_fname}"),
                destfile = "data/data_output/{lambda_df_fname}")
}
mcvr <- readr::read_csv(glue::glue("data/data_output/{lambda_df_fname}"))


# fecundity ---------------------------------------------------------------

trt_on_fecundity <- 
  mcvr %>% 
  dplyr::filter(manipulated_vr %in% c("fecundity", "ambient")) %>% 
  dplyr::select(manipulated_vr, degree.days, vwc, snow.days, heat.f, water.f, lambda)

fecundity_contrasts <- 
  trt_on_fecundity %>% 
  dplyr::select(-manipulated_vr) %>% 
  base::split(.[, c("degree.days", "vwc", "snow.days")]) %>% 
  lapply(FUN = function(x) {
    x <- 
      x %>% 
      dplyr::rename(heat = starts_with("heat"),
                    water = starts_with("water"))
    
    my_df <- 
      tibble(degree.days = unique(x$degree.days), 
             vwc = unique(x$vwc), 
             snow.days = unique(x$snow.days), 
             h = x$lambda[x$heat == 1 & x$water == 0],
             w = x$lambda[x$heat == 0 & x$water == 1],
             hw = x$lambda[x$heat == 1 & x$water == 1],
             a = x$lambda[x$heat == 0 & x$water == 0],
             contrast_h = ((h + hw) / 2) - ((w + a) / 2),
             contrast_w = ((w + hw) / 2) - ((h + a) / 2),
             contrast_hw = (hw) - ((h + w + a) / 3)) %>% 
      dplyr::select(degree.days, vwc, snow.days, starts_with("contrast")) %>% 
      pivot_longer(cols = c(contrast_h, contrast_w, contrast_hw),
                   names_to = "trt_contrast",
                   values_to = "delta_lambda") %>% 
      dplyr::group_by(degree.days, vwc, snow.days, trt_contrast) %>% 
      dplyr::summarize(delta_lambda_mean = mean(delta_lambda),
                       delta_lambda_lwr = quantile(delta_lambda, probs = 0.025),
                       delta_lambda_upr = quantile(delta_lambda, probs = 0.975))
  }) %>% 
  bind_rows()


# growth ------------------------------------------------------------------

trt_on_growth <- 
  mcvr %>% 
  filter(manipulated_vr %in% c("growth", "ambient")) %>% 
  dplyr::select(manipulated_vr, degree.days, vwc, snow.days, heat.g, water.g, lambda)


growth_contrasts <- 
  trt_on_growth %>% 
  dplyr::select(-manipulated_vr) %>% 
  base::split(.[, c("degree.days", "vwc", "snow.days")]) %>% 
  lapply(FUN = function(x) {
    x <- 
      x %>% 
      dplyr::rename(heat = starts_with("heat"),
                    water = starts_with("water"))
    
    my_df <- 
      tibble(degree.days = unique(x$degree.days), 
             vwc = unique(x$vwc), 
             snow.days = unique(x$snow.days), 
             h = x$lambda[x$heat == 1 & x$water == 0],
             w = x$lambda[x$heat == 0 & x$water == 1],
             hw = x$lambda[x$heat == 1 & x$water == 1],
             a = x$lambda[x$heat == 0 & x$water == 0],
             contrast_h = ((h + hw) / 2) - ((w + a) / 2),
             contrast_w = ((w + hw) / 2) - ((h + a) / 2),
             contrast_hw = (hw) - ((h + w + a) / 3)) %>% 
      dplyr::select(degree.days, vwc, snow.days, starts_with("contrast")) %>% 
      pivot_longer(cols = c(contrast_h, contrast_w, contrast_hw),
                   names_to = "trt_contrast",
                   values_to = "delta_lambda") %>% 
      dplyr::group_by(degree.days, vwc, snow.days, trt_contrast) %>% 
      dplyr::summarize(delta_lambda_mean = mean(delta_lambda),
                       delta_lambda_lwr = quantile(delta_lambda, probs = 0.025),
                       delta_lambda_upr = quantile(delta_lambda, probs = 0.975))
  }) %>% 
  bind_rows()


# survivorship ------------------------------------------------------------

trt_on_survivorship <- 
  mcvr %>% 
  filter(manipulated_vr %in% c("survivorship", "ambient")) %>% 
  dplyr::select(manipulated_vr, degree.days, vwc, snow.days, heat.s, water.s, lambda)


survivorship_contrasts <- 
  trt_on_survivorship %>% 
  dplyr::select(-manipulated_vr) %>% 
  base::split(.[, c("degree.days", "vwc", "snow.days")]) %>% 
  lapply(FUN = function(x) {
    x <- 
      x %>% 
      dplyr::rename(heat = starts_with("heat"),
                    water = starts_with("water"))
    
    my_df <- 
      tibble(degree.days = unique(x$degree.days), 
             vwc = unique(x$vwc), 
             snow.days = unique(x$snow.days), 
             h = x$lambda[x$heat == 1 & x$water == 0],
             w = x$lambda[x$heat == 0 & x$water == 1],
             hw = x$lambda[x$heat == 1 & x$water == 1],
             a = x$lambda[x$heat == 0 & x$water == 0],
             contrast_h = ((h + hw) / 2) - ((w + a) / 2),
             contrast_w = ((w + hw) / 2) - ((h + a) / 2),
             contrast_hw = (hw) - ((h + w + a) / 3)) %>% 
      dplyr::select(degree.days, vwc, snow.days, starts_with("contrast")) %>% 
      pivot_longer(cols = c(contrast_h, contrast_w, contrast_hw),
                   names_to = "trt_contrast",
                   values_to = "delta_lambda") %>% 
      dplyr::group_by(degree.days, vwc, snow.days, trt_contrast) %>% 
      dplyr::summarize(delta_lambda_mean = mean(delta_lambda),
                       delta_lambda_lwr = quantile(delta_lambda, probs = 0.025),
                       delta_lambda_upr = quantile(delta_lambda, probs = 0.975))
  }) %>% 
  bind_rows()


# figures -----------------------------------------------------------------


# fecundity manipulation; degree days, vwc --------------------------------

ggplot(fecundity_contrasts %>% filter(snow.days == 0), aes(x = degree.days, y = vwc, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# fecundity manipulation; snow days, vwc --------------------------------

ggplot(fecundity_contrasts %>% filter(degree.days == 0), aes(x = snow.days, y = vwc, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# fecundity manipulation; snow days, degree.days --------------------------------

ggplot(fecundity_contrasts %>% filter(vwc == 0), aes(x = degree.days, y = snow.days, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# growth manipulation; degree days, vwc --------------------------------

ggplot(growth_contrasts %>% filter(snow.days == 0), aes(x = degree.days, y = vwc, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# growth manipulation; snow days, vwc --------------------------------

ggplot(growth_contrasts %>% filter(degree.days == 0), aes(x = snow.days, y = vwc, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# growth manipulation; snow days, degree.days --------------------------------

ggplot(growth_contrasts %>% filter(vwc == 0), aes(x = degree.days, y = snow.days, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# survivorship manipulation; degree days, vwc --------------------------------

ggplot(survivorship_contrasts %>% filter(snow.days == 0), aes(x = degree.days, y = vwc, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# survivorship manipulation; snow days, vwc --------------------------------

ggplot(survivorship_contrasts %>% filter(degree.days == 0), aes(x = snow.days, y = vwc, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()

# survivorship manipulation; snow days, degree.days --------------------------------

ggplot(survivorship_contrasts %>% filter(vwc == 0), aes(x = degree.days, y = snow.days, fill = delta_lambda_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(trt_contrast)) +
  coord_fixed() +
  theme_minimal()
