#pollen analysis
#this is for the ragweed pollen sampling that was conducted near patches of ragweed on 9/16/15, 9/17/15, and 9/22/15
#the purpose is to see how ragweed concentrations diminish as a function of distance from source and from other covariates
#several other pilots had been conducted, but this is the final one

#set up work environment
rm(list = ls())

library(ggplot2)
library(dplyr)

#setwd("Q:/Ibanez Lab/Dan Katz/post doc search/NIH F32/ragweed air sampling")

#load in pollen count data
p <- read.csv("Pollen Grain Counting Record_151019.csv")
p <- subset(p, !is.na(p$Count.by.Dan))  #taking out the unmeasured stuff for now
head(p)

#getting distances for each transect on the filter
p$filter_y
p$filter_y_start <- as.numeric(substr(as.character(p$filter_y), 1,2))
p$filter_y_end <- as.numeric(substr(as.character(p$filter_y), 6,7))
p$filter_dist <- abs(p$filter_y_start - p$filter_y_end)

p$p_per_mm2 <- p$grains_dk / (p$filter_dist * 0.5)
hist(p$p_per_mm2)

#mean and sd of each filter
s_summary <- tbl_df(p) %>% 
  group_by(round, sampler) %>% 
  summarize(p_per_mm2_mean = mean(p_per_mm2, na.rm = TRUE),
            p_per_mm2_sd = sd(p_per_mm2, na.rm = TRUE))

#coversion to cubic meters - will do individually for air samplers
#it might be that 1/4 of area is a bad estimate - should do the basic geometry
#a = pie *37mm (?) ^2 
s_summary$p_per_m2_mean <- (s_summary$p_per_mm2_mean * 4 * 4 ) /  #1/4 of slide has transects, 1/4 of that area is viewed
                            ((2/1000)*120)  # 2 L per min (= o.002 m3), 120 min,

s_summary$p_per_m2_sd <- (s_summary$p_per_mm2_sd * 4 * 4 ) /  #1/4 of slide has transects, 1/4 of that area is viewed
  ((2/1000)*120)  # 2 L per min (= o.002 m3), 120 min,


#load in collector data and link it to s_summary
c <- read.csv("samplinglocations151019.csv")
s_summary <- left_join(s_summary, c, by = c("round", "sampler"))

#load in wind data  
w <- read.csv("winddata151019.csv")


#####################  
#pollen visualization
#####################
s_summary$p_per_m2_mean

#pollen as a function of distance
ggplot(s_summary, aes(x= distance_measured, y = p_per_m2_mean, ymin = p_per_m2_mean - p_per_m2_sd,
                      ymax = p_per_m2_mean + p_per_m2_sd)) + geom_point() + 
  ylab("distance from patch") + xlab("pollen grains per m3") + theme_bw() +
  facet_wrap( ~ round) +
  geom_errorbar()

#wind roses for each site
head(w)
#speed over time
ggplot(w, aes(x = min_since_start, y = speed)) + geom_point() + facet_wrap( ~ round)

#direction 
ggplot(w, aes(x = direction, y = speed)) + geom_point() + facet_wrap( ~ round) + 
  coord_polar() + theme_bw()


#pollen as a function of wind direction and speed

#pollen grains as a function of interior vs. exterior of slide
ggplot(p, aes(x = as.factor(transect), y = p_per_mm2)) + geom_boxplot() 
ggplot(p, aes(x = as.factor(transect), y = p_per_mm2)) + geom_jitter() 




ggplot(p, aes(x = as.factor(Sample), y = p)) + geom_boxplot()

ggplot(pc, aes(x = distance_measured, y = p, color = as.factor(sampler))) + geom_point() + theme_bw()
ggplot(pc, aes(x = distance_measured, y = p, color = Transect)) + geom_point()


ggplot(pc, aes(x = distance_measured))


ggplot(pc, aes(x=distance_measured, y = p_mean.avg)) + geom_boxplot()


#wind visualization
head(w)


ggplot(w, aes(x=min_since_start, y = speed)) + geom_point() + facet_wrap( ~ round)
ggplot(w, aes(x=direction)) + geom_histogram() + facet_wrap( ~ round) + coord_polar()


#effects of wind on pollen concentrations
