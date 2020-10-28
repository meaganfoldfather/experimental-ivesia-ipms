# CH2 Revision - Figuring out model structure in response to reviewers
# 20201022

# libraries
library(lme4)
library(lmerTest)
library(ggplot2)
library(car)
library(popdemo)
library(dplyr)
library(tidyr)
library(viridis)
library(MASS)
library(MuMIn)
library(rstan)
library(brms)
library(purrr)
library(broom)
library(devtools)
library(furrr)
library(tidyverse)
library(scales)
library(cowplot)
library(tictoc)

#### Bring in vital rate data --> I.df ####
vr.mc<-read.csv("https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/VitalRates_Microclimate.csv", stringsAsFactors = F); head(vr.mc)

#Add in experimental metadata
plot.chars<-read.csv("https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/trt.plot.data.csv")
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
microclimate<-read.csv("https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/Microclimate.csv")

head(microclimate)
head(plot.chars)
mc.trt<-merge(microclimate, plot.chars)
head(mc.trt)
mc<-mc.trt %>%
  filter(type == "EXPERIMENTAL") %>% 
  filter(trt == "A") %>%
#  group_by(site, year) %>% 
  group_by(site) %>% 
  summarise(vwc = round(mean(vwc, na.rm=T), 0), degree.days = round(mean(degree.days, na.rm=T),0), snow.days = round(mean(snow.days, na.rm=T,0)))
#mutate(year = as.factor(year))
mc

# join vr and microclimate data
#df <- left_join(df, mc, by = c(site = "site", t1 = "year"))
df <- left_join(df, mc, by = c(site = "site"))
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

#### Change to Sum coding for year ####
contrasts(df$t1) = c("contr.sum", "contr.poly")
#### Survival ####
tic()
mod = brms::brm(
  surv ~ size.s*degree.days*heat*water*vwc + t1 + (1|site/plot),
  data = df,
  family = 'bernoulli',
  prior = set_prior('normal(0, 3)'),
  iter = 2000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
toc()

summary(mod)
plot(mod)
bayes_R2(mod) 
pp_check(mod,nsamples = 50)

#saveRDS(object = mod, file = "surv_mod.rds")

# Trying T1 (year) as a fixed effect
# m1 = surv ~ size.s*degree.days*heat*water + size.s*vwc*heat*water + t1 + (1|site/plot)
#adapt_delta = .995, max_treedepth = 15
# No divergences or treedepth issue! Took 32 mins to run! 

# m2 = surv ~ size.s*degree.days*heat*water*vwc + t1 + (1|site/plot)
#adapt_delta = .995, max_treedepth = 15
# No divergences or treedepth issue! Took 91 mins to run! 


# Average climate modes # (Did these models with dd* vwc interaction for first 4)
# model with random effects of site/plot + 1/year, no snow, had 24 divergent transitions and 2 transition that exceed treedepth (alpha = .8, tree depth specified to 10) took 30 mins to run
# model with random effects of site/plot + 1/year, no snow, had 25 divergent transitions and 0 transition that exceed treedepth (alpha = .9, tree depth specified to 15) took 38 mins to run
# model with random effects of site/plot + 1/year, no snow, had 3 divergent transitions and 0 transition that exceed treedepth (alpha = .99, tree depth specified to 15) took 11.4 hrs to run
# model with random effects of 1/plot + 1/year, no snow, had 35 divergent transitions and 0 transition that exceed treedepth (alpha = .9, tree depth specified to 15) took 24 mins to run


# no interactions btw climate vars for these models #
# model with random effects of 1/plot + 1/year, no snow, had 5 divergent transitions and 0 transition that exceed treedepth (alpha = .99, tree depth specified to 15) took 59 mins to run
# model with random effects of 1/plot + 1/year, no snow, had 5 divergent transitions and 0 transition that exceed treedepth (alpha = .99, tree depth specified to 15) took 43 mins to run
# model with random effects of 1/plot + 1/year, no snow, had 5 divergent transitions and 0 transition that exceed treedepth (alpha = .995, tree depth specified to 15) took 44 mins to run


#Year-specific climate #
# model with random effects of site/plot + 1/year had 39 divergent transitions, alpha = .9 and took 22 mins to run 
# model with random effects of site/plot + 1/year had 37 divergent transitions, alpha = .95 and took 29 mins to run 
# model with random effects of site/plot + 1/year had 9 divergent transitions, alpha = .99, 1238 transitions exceeded max treedepth and took 1hr to run 
# model with random effects of 1/plot + 1/year had 1 divergent transitions, alpha = .99, 25 transitions exceeded max treedepth and took 1hr and 50 mins to run 
# model with random effects of 1/plot + 1/year had 5 divergent transitions, alpha = .995,tree depth specified to 15 & no transitions exceeded max treedepth and took 10.1 hrs to run 
# model with random effects of 1/site/plot + 1/year, but with snow had 85 divergent transitions (alpha = .8,tree depth specified to 10) and took 1.5 hrs to run 
# model with random effects of site/plot + 1/year, no snow included, but dd*vwc had 61 divergent transitions (alpha = .9, tree depth specified to 10) took 17 mins to run
# model with random effects of site/plot + 1/year, no snow included, had 70 divergent transitions (alpha = .9, tree depth specified to 10) took 16 mins to run
# model with random effects of site/plot + 1/year, no snow included, had 1 divergent transitions and 1 transition that exceed treedepth (alpha = .95, tree depth specified to 10) took 52 mins to run
# model with random effects of site/plot + 1/year, no snow included, had 7 divergent transitions and 247 transitions that exceed treedepth (alpha = .99, tree depth specified to 10) took 51 mins to run
# model with random effects of site/plot + 1/year, no snow included, had 5 divergent transitions and 0 transitions that exceed treedepth (alpha = .99, tree depth specified to 15) took 1.3 hours to run


#saveRDS(object = mod, file = "surv_mod_add_ave.rds")

#### Hurdle Reproduction ####
df[which(df$fprobNext == 0), "seedNext"] <- 0
tic()
model.formula <- bf(seedNext ~ sizeNext.s*degree.days*heat*water*vwc + t1 + (1|site/plot), hu ~ sizeNext.s*degree.days*heat*water*vwc + t1 + (1|site/plot))
hurdleRep = brms::brm(formula = model.formula,
  data = df,
  family = "hurdle_poisson",
  prior = set_prior('normal(0, 3)'),
  iter = 1000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = .85, max_treedepth = 15)
)
toc()

summary(hurdleRep)
plot(hurdleRep)
bayes_R2(hurdleRep) #.67!
pp_check(hurdleRep,nsamples = 50)

saveRDS(object = hurdleRep, file = "hurdle_mod_add_ave_sum.rds")

# model with fixed effect of year, random effects of site/plot, no snow, had 0 divergent transitions but 753 transition that exceed treedepth (alpha = .8, tree depth specified to 10) and  took 45 mins to run --> treedpeth issues are more a efficiency concern

#### Growth ####
tic()
df$sizeNext.s <- as.integer(df$sizeNext.s)
mod.grwth = brms::brm(
  sizeNext.s ~ size.s*degree.days*heat*water*vwc + t1 + (1|site/plot),
  data = df,
  family = 'gaussian',
  prior = set_prior('normal(0, 3)'),
  iter = 1500,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
toc()

summary(mod.grwth)
plot(mod.grwth)
bayes_R2(mod.grwth) 
pp_check(mod.grwth,nsamples = 50)

# model with fixed effect of year, random effects of site/plot, no snow, had 5 divergent transitions and 0  transition that exceed treedepth (alpha = .8, tree depth specified to 10) and  took 12 mins to run;also had low bulk and tail effective sample size, indicating the chains need more iters
# model with fixed effect of year, random effects of site/plot, no snow, had 0 divergent transitions and 5 transition that exceed treedepth (alpha = .9, tree depth specified to 10) and  took 14 mins to run;also had low bulk effective sample size, indicating the chains need more iters
# model with fixed effect of year, random effects of site/plot, no snow, had 1 divergent transitions and 0 transition that exceed treedepth (alpha = .9, tree depth specified to 15) and  took 16 mins to run;also had low bulk effective sample size, indicating the chains need more iters
# model with fixed effect of year, random effects of site/plot, no snow, had 0 divergent transitions and 0 transition that exceed treedepth (alpha = .99, tree depth specified to 15) and  took 42 mins to run;but had low bulkand tail effective sample size, indicating the chains need more iters
#Final model with 1500 iters took 42 mins to run

saveRDS(object = mod.grwth, file = "grwth_mod_add_ave_sum.rds")

#### Establishment ####
establishment.df<- 
  df %>% 
  mutate(vwc = as.numeric(vwc), degree.days = as.numeric(degree.days), snow.days = as.numeric(snow.days)) %>% 
  group_by(t1, plot, site, heat, water, vwc, degree.days, snow.days) %>% 
  summarise(recruit = sum(new), seeds.avaliable = sum(seed, na.rm = T)) %>%
  filter(seeds.avaliable > 0) %>% 
  mutate(recruit = ifelse(recruit > seeds.avaliable, yes = seeds.avaliable, no = recruit)) 
establishment.df

recruit.mod = brms::brm(recruit | trials(seeds.avaliable) ~ degree.days*heat*water*vwc + t1 + (1|site), data = establishment.df, family = 'binomial', prior = set_prior('normal(0, 3)'), iter = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = .95, max_treedepth = 10)) 

summary(recruit.mod) 
plot(recruit.mod)
bayes_R2(recruit.mod) 
pp_check(recruit.mod,nsamples = 50)

saveRDS(object = recruit.mod, file = "recruit_mod_add_ave_sum.rds")

#not_summed_establishment_model<-readRDS("Recruit_mod_add_ave.rds")
#summary(not_summed_establishment_model)
