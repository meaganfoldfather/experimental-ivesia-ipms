#Intent: Analysis for CH2 - Analysis to look at Specific Responses of Each Vital Rate to the Experimental Manipulations
#Date Last Edited: 20190915

#### Library ####
library(tidyverse)
library(rstan)
library(brms)
library(purrr)
library(devtools)
library(furrr)

#### Bring in vital rate data --> I.df ####
vr.mc<-read.csv("https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/VitalRates_Microclimate.csv", stringsAsFactors = F); head(vr.mc)

#Add in experimental metdata
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

##### Define Functions of vital rates ####
#load in models
dir.create("data/data_output", recursive = TRUE)
download.file(url = "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/surv_mod_add_ave.rds", 
              destfile = "data/data_output/surv_mod_add_ave.rds")
survival_model <- readRDS("data/data_output/surv_mod_add_ave.rds")

download.file(url = "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/grwth_mod_add_ave.rds",
              destfile = "data/data_output/grwth_mod_add_ave.rds")
growth_model <- readRDS("data/data_output/grwth_mod_add_ave.rds")

download.file(url = "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/recruit_mod_add_ave.rds",
              destfile = "data/data_output/recruit_mod_add_ave.rds")
establishment_model <- readRDS("data/data_output/recruit_mod_add_ave.rds")
establishment_model$data <- establishment_model$data %>% mutate(seeds.avaliable = ifelse(seeds.avaliable == 0, yes = NA, no = seeds.avaliable)) # removed zeros from $data from the model object because error trapping for number of trials became more strict in brms 2.10

download.file(url = "https://earthlab-mkoontz.s3-us-west-2.amazonaws.com/experimental-ivesia-ipms/hurdle_mod_add_ave.rds",
              destfile = "data/data_output/hurdle_mod_add_ave.rds")
hurdleRep <- readRDS("data/data_output/hurdle_mod_add_ave.rds")

# kernel construction
min.size=min(df[,c("size.s", "sizeNext.s")], na.rm=T)
max.size=max(df[, c("size.s", "sizeNext.s")], na.rm=T)
n=100 # number of cells in the matrix
b=min.size+c(0:n)*(max.size-min.size)/n # boundary points
y=0.5*(b[1:n]+b[2:(n+1)]) # mesh points
h=y[2]-y[1] # step size

#1. survival probability function
s.x=function(x, degree.days = 0, vwc = 0, snow.days = 0, heat = 0, water = 0) {
  new.data <- data.frame(size.s = x, degree.days = degree.days, vwc = vwc, snow.days = snow.days, heat = heat, water = water)
  survivorship <- fitted(survival_model, newdata = new.data, re_formula = NA, scale = "response", summary = FALSE)
  return(survivorship)
}

# 2. growth function
g.yx=function(xp,x,degree.days = 0, vwc = 0, snow.days = 0, heat = 0, water = 0) {
  new.data <- data.frame(size.s = x, degree.days = degree.days, vwc = vwc, snow.days = snow.days, heat = heat, water = water)
  sizeNext.s <- predict(growth_model, newdata = new.data, re_formula = NA, scale = "response", summary = FALSE) 
  
  
  growth <- apply(sizeNext.s, MARGIN = 2, FUN = function(sizeNext.s_column) {approxfun(density(sizeNext.s_column))(xp)})
  # if the probability is so small that it doesn't exist in our empirical
  # probability distribution (and thus the approxfun(density(x)) function 
  # returns NA), just make the probability 0
  growth[is.na(growth)] <- 0 
  
  return(growth)
}

# 3. reproduction function
f.yx <- function(x, xp ,degree.days = 0, vwc = 0, snow.days = 0, heat = 0, water = 0) {
  new.data <- data.frame(sizeNext.s = x, degree.days = degree.days, vwc = vwc, snow.days = snow.days, heat = heat, water = water)
  
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
kernel_construction <- function(x,degree.days = 0, vwc = 0, snow.days = 0, heat.s = 0, water.s = 0,heat.g = 0, water.g = 0, heat.f = 0, water.f = 0){
  
  F.mat = h*simplify2array(f.yx(x, x, degree.days, vwc, snow.days, heat = heat.f, water = water.f))
  S = s.x(x, degree.days, vwc, snow.days, heat = heat.s, water = water.s)
  S = array(t(S), dim = c(1, 100, 2000))
  G = h * g.yx(x, x, degree.days, vwc, snow.days, heat = heat.g, water = water.g)
  G = array(G, dim = c(100, 100, 1))
  
  P = array(dim=c(100,100,2000))
  for(k in 1:n) P[,k,]=G[,k,]*S[,k,] # growth/survival matrix
  # full matrix
  K=P+F.mat
  
  lambda <- apply(K, MARGIN = 3, FUN = function(k) Re(eigen(k)$values[1]))
  return(lambda)
}

#Treatment Effects Across Climatic Gradients
plan(multiprocess)
start_time <- Sys.time()
start_time

vital_effects <- 
  data.frame(
    heat.f = c(1,0,1, rep(0,7)),
    water.f = c(0,1,1, rep(0,7)), 
    heat.g = c(0,0,0,1,0,1,0,0,0,0), 
    water.g = c(0,0,0,0,1,1,0,0,0,0), 
    heat.s = c(rep(0,6),1,0,1,0), 
    water.s = c(rep(0,7),1,1,0))

length.out <- 5
microclimate_effects <-data.frame(expand.grid(degree.days = seq(-2, 2, length.out = length.out), vwc = seq(-2, 2, length.out = length.out), snow.days =  seq(-2, 2, length.out = length.out)))

mc_vr_effects <-
  microclimate_effects %>% 
  dplyr::mutate(vital_effects = list(vital_effects)) %>% 
  unnest(vital_effects) %>% 
  mutate(lambda = future_pmap(., kernel_construction, y)) %>%
  unnest()

end_time <- Sys.time()
end_time - start_time

#write_csv() 

