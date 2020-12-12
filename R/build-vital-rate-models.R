# CH2 Revision - Figuring out model structure in response to reviewers
# 20201022

# libraries
library(dplyr)
library(tidyr)
library(brms)
library(purrr)
library(tictoc)
library(glue)
library(ggplot2)

overwrite <- FALSE

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

#### Bring in vital rate data --> I.df ####
vr.mc<-read.csv(glue::glue("{remote_source}/VitalRates_Microclimate.csv"), stringsAsFactors = F); head(vr.mc)

#Add in experimental metadata
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
survival_model = brms::brm(
  surv ~ size.s*heat*water*degree.days + size.s*heat*water*vwc + t1 + (1|site/plot),
  data = df,
  family = 'bernoulli',
  prior = set_prior('normal(0, 3)'),
  iter = 2000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = .99, max_treedepth = 15)
)
toc() # 1126.566 seconds on Alien

summary(survival_model)
plot(survival_model)
bayes_R2(survival_model) 
pp_check(survival_model, nsamples = 50)

if(overwrite | !file.exists(glue::glue("data/data_output/{surv_mod_fname}"))) {
  saveRDS(object = survival_model, file = glue::glue("data/data_output/{surv_mod_fname}"))
  
  system2(command = "aws", args = glue::glue("s3 cp data/data_output/{surv_mod_fname} {remote_target}/{surv_mod_fname} --acl public-read"))
}

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

#### Hurdle Reproduction ####
set.seed(1541)
df[which(df$fprobNext == 0), "seedNext"] <- 0
(start <- Sys.time())
model.formula <- bf(seedNext ~ sizeNext.s*heat*water*degree.days + sizeNext.s*heat*water*vwc + t1 + (1|site/plot), 
                    hu ~ sizeNext.s*heat*water*degree.days + sizeNext.s*heat*water*vwc + t1 + (1|site/plot))
hurdleRep = brms::brm(formula = model.formula,
                      data = df,
                      family = "hurdle_negbinomial",
                      prior = set_prior('normal(0, 3)'),
                      iter = 1000,
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.9)
)

(difftime(time1 = Sys.time(), time2 = start, units = "mins"))

summary(hurdleRep)
plot(hurdleRep)
bayes_R2(hurdleRep) #.67!
pp_check(hurdleRep, nsamples = 50) + coord_cartesian(xlim = c(0, 50))

if(overwrite | !file.exists(glue::glue("data/data_output/{hurdle_mod_fname}"))) {
  saveRDS(object = hurdleRep, file = glue::glue("data/data_output/{hurdle_mod_fname}"))
  
  system2(command = "aws", args = glue::glue("s3 cp data/data_output/{hurdle_mod_fname} {remote_target}/{hurdle_mod_fname} --acl public-read"))
}


# model with fixed effect of year, random effects of site/plot, no snow, had 0 divergent transitions but 753 transition that exceed treedepth (alpha = .8, tree depth specified to 10) and  took 45 mins to run --> treedpeth issues are more a efficiency concern

#### Growth ####
set.seed(822)
(start <- Sys.time())
growth_model = brms::brm(
  sizeNext ~ size.s*heat*water*degree.days + size.s*heat*water*vwc + t1 + (1|site/plot),
  data = df,
  family = 'negbinomial',
  prior = set_prior('normal(0, 3)'),
  iter = 1000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.9)
)
difftime(time1 = Sys.time(), time2 = start, units = "mins")

summary(growth_model)
plot(growth_model)
bayes_R2(growth_model) 
growth_pp_check <- pp_check(growth_model)
growth_pp_check + coord_cartesian(xlim = c(0, 100))

# model with fixed effect of year, random effects of site/plot, no snow, had 5 divergent transitions and 0  transition that exceed treedepth (alpha = .8, tree depth specified to 10) and  took 12 mins to run;also had low bulk and tail effective sample size, indicating the chains need more iters
# model with fixed effect of year, random effects of site/plot, no snow, had 0 divergent transitions and 5 transition that exceed treedepth (alpha = .9, tree depth specified to 10) and  took 14 mins to run;also had low bulk effective sample size, indicating the chains need more iters
# model with fixed effect of year, random effects of site/plot, no snow, had 1 divergent transitions and 0 transition that exceed treedepth (alpha = .9, tree depth specified to 15) and  took 16 mins to run;also had low bulk effective sample size, indicating the chains need more iters
# model with fixed effect of year, random effects of site/plot, no snow, had 0 divergent transitions and 0 transition that exceed treedepth (alpha = .99, tree depth specified to 15) and  took 42 mins to run;but had low bulkand tail effective sample size, indicating the chains need more iters
#Final model with 1500 iters took 42 mins to run

if(overwrite | !file.exists(glue::glue("data/data_output/{growth_mod_fname}"))) {
  saveRDS(object = growth_model, file = glue::glue("data/data_output/{growth_mod_fname}"))
  
  system2(command = "aws", args = glue::glue("s3 cp data/data_output/{growth_mod_fname} {remote_target}/{growth_mod_fname} --acl public-read"))
}

#### Establishment ####
establishment.df<- 
  df %>% 
  mutate(vwc = as.numeric(vwc), degree.days = as.numeric(degree.days), snow.days = as.numeric(snow.days)) %>% 
  group_by(t1, plot, site, heat, water, vwc, degree.days, snow.days) %>% 
  summarise(recruit = sum(new), seeds.avaliable = sum(seed, na.rm = T)) %>%
  filter(seeds.avaliable > 0) %>% 
  mutate(recruit = ifelse(recruit > seeds.avaliable, yes = seeds.avaliable, no = recruit)) 
establishment.df

establishment_model = brms::brm(recruit | trials(seeds.avaliable) ~ heat*water*degree.days + heat*water*vwc + t1 + (1|site), data = establishment.df, family = 'binomial', prior = set_prior('normal(0, 3)'), iter = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = .95, max_treedepth = 10)) 

summary(establishment_model) 
plot(establishment_model)
bayes_R2(establishment_model) 
pp_check(establishment_model, nsamples = 50)

if(overwrite | !file.exists(glue::glue("data/data_output/{establishment_mod_fname}"))) {
  saveRDS(object = establishment_model, file = glue::glue("data/data_output/{establishment_mod_fname}"))
  
  system2(command = "aws", args = glue::glue("s3 cp data/data_output/{establishment_mod_fname} {remote_target}/{establishment_mod_fname} --acl public-read"))
}

survival_model
hurdleRep
growth_model
establishment_model
