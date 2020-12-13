# Get estimated lambda for each site (using the site-specific microclimate variables)

#### Library ####
library(dplyr)
library(brms)
library(purrr)
library(furrr)
library(glue)
library(ggplot2)
library(data.table)

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

lambda_df_fname <- "site_specific_additive.csv"

overwrite <- TRUE

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



#### Begin prepare the plant measurements
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

# filter out rows with 0 for size or size next
# Just a couple plants had no leaves when sampled, but were still alive.
# unclear what this really means biologically (perhaps just an artifact of
# when sampling occurred during the season) so we drop them
df <- df[-which(df$size == 0),]
df <- df[-which(df$sizeNext == 0),]

# One plant (2118.OS.15.0) was recorded as having a sizeNext of 0.3 in t1=2015
# which is almost certainly a typo (should almost certainly be 3, since previous
# size was 3). We will drop this row, and the entry for this plant in the next
# year, when the error propagates and the "size" column then becomes 0.3.
df <- df[-which(df$size == 0.3), ]
df <- df[-which(df$sizeNext == 0.3), ]

size_vec <- c(df$size, df$sizeNext)
mean_size <- mean(size_vec, na.rm = TRUE)
sd_size <- sd(size_vec, na.rm = TRUE)

df$size.s <- (df$size - mean_size) / sd_size
df$sizeNext.s <- (df$sizeNext - mean_size)/ sd_size

# scale climate vectors
df$vwc<- scale(df$vwc)
df$degree.days <- scale(df$degree.days)
df$snow.days<- scale(df$snow.days)

# scale elevation vectors
df$elevation.s <- scale(df$elevation)

#### End prepare the plant measurements



#### Begin prepare the IPM infrastructure
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
  sizeNext <- predict(growth_model, newdata = new.data, re_formula = NA, scale = "response", summary = FALSE) 
  
  sizeNext.s <- (sizeNext - mean_size) / sd_size
  
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


#### end IPM infrastructure setup

#### get site-level microclimate
site_mc_scaled <-
  df %>%
  group_by(site) %>%
  summarize(vwc = unique(vwc)[1], degree.days = unique(degree.days)[1])

#Treatment Effects Across Climatic Gradients
# Running on the Alien with 10 workers on 2020-05-06; runtime = 

plan(multiprocess, workers = 6)
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

# Now create a dataframe that includes the site-specific microclimate as well
# as the experimental effects we want to implement on lambda (as defined by
# the vital_effects matrix above)
# Then, we iterate over each row using future_pmap() [using dplyr::select(., -site)
# to remove the site variable from the values in the row being passed to the
# kernel_construction() function] and build the IPM to get lambda for that
# particular combination of microclimate and experimental treatment (either 
# applied to all vital rates at once, or one at a time). The result is a very
# large dataframe with 2000 values of lambda per combination of microclimate
# variables (of which there is one combo per site) and per experimental manipulation

site_specific_lambda <- 
  site_mc_scaled %>%
  dplyr::mutate(vital_effects = list(vital_effects)) %>% 
  tidyr::unnest(vital_effects) %>%
  mutate(lambda = furrr::future_pmap(.l = dplyr::select(., -site), kernel_construction, y, .options = furrr::furrr_options(seed = TRUE))) %>%
  tidyr::unnest(cols = lambda)

site_specific_lambda <- 
  site_specific_lambda %>% 
  dplyr::mutate(manipulated_vr = case_when(heat.f == 1 & heat.g == 1 & heat.s == 1 & water.f == 1 & water.g == 1 & water.s == 1 ~ "hw",
                                           heat.f == 1 & heat.g == 1 & heat.s == 1 ~ "heat",
                                           water.f == 1 & water.g == 1 & water.s == 1 ~ "water",
                                           heat.f == 1 | water.f == 1 ~ "fecundity",
                                           heat.g == 1 | water.g == 1 ~ "growth",
                                           heat.s == 1 | water.s == 1 ~ "survivorship",
                                           TRUE ~ "ambient"))

if(overwrite | !file.exists(glue::glue("data/data_output/{lambda_df_fname}"))) {
  data.table::fwrite(x = site_specific_lambda, file = here::here("data", "data_output", lambda_df_fname))
  
  system2(command = "aws", args = glue::glue("s3 cp data/data_output/{lambda_df_fname} {remote_target}/{lambda_df_fname} --acl public-read"))
}
(end_time <- Sys.time())

(base::difftime(end_time, start_time, units = "mins"))
