#pollen analysis
#this is for the ragweed pollen sampling that was conducted near patches of ragweed on 9/16/15, 9/17/15, and 9/22/15
#the purpose is to see how ragweed concentrations diminish as a function of distance from source and from other covariates
#several other pilots had been conducted, but this is the final one

#note: need to get rid of sd, due to unequal transect length

#set up work environment
rm(list = ls())

library(ggplot2)
library(dplyr)

setwd("Q:/Ibanez Lab/Dan Katz/post doc search/NIH F32/ragweed air sampling/ragweed-pollen-fall-2015")

#load in pollen count data
p <- read.csv("Pollen Grain Counting Record_151026.csv")
#p <- subset(p, !is.na(p$grains_dk))  #taking out the unmeasured stuff for now
head(p)

#getting distances for each transect on the filter
p$filter_y
p$filter_y_start <- as.numeric(substr(as.character(p$filter_y), 1,2))
p$filter_y_end <- as.numeric(substr(as.character(p$filter_y), 6,7))
p$filter_dist <- abs(p$filter_y_start - p$filter_y_end)

p$p_per_mm2 <- p$pollen / (p$filter_dist * 0.5)
hist(p$p_per_mm2)

#mean and sd of each filter
s_summary <- tbl_df(p) %>% 
  group_by(round, sampler) %>% 
  summarize(sum_p = sum(pollen, na.rm = TRUE),
            sum_area = sum(filter_dist * 0.5))

#coversion to cubic meters 
a <- 3.14159*((37/2)^2)  #filter area (mm2), diameter is 37mm
s_summary$p_per_m3 <- (s_summary$sum_p * (a/s_summary$sum_area)) /  #estimated number of pollen grains per slide
                            ((3/1000)*120)  # 3 L per min (= o.003 m3), 120 min,


#load in collector data and link it to s_summary
c <- read.csv("samplinglocations151019.csv")
s_summary <- left_join(s_summary, c, by = c("round", "sampler"))

#load in wind data  
w <- read.csv("winddata151019.csv")
w_summary <- tbl_df(w) %>% filter(ss == 0) %>% group_by(round) %>%
  summarise(direction_mean = mean(direction), direction_sd = sd(direction),
            speed_mean = mean(speed), speed_sd = sd(speed))

#####################  
#pollen visualization
#####################
s_summary$p_per_m2_mean

#pollen as a function of distance
ggplot(s_summary, aes(x= distance_measured, y = p_per_m3)) + geom_point() + 
  xlab("distance from patch") + ylab("pollen grains per m3") + theme_bw() +
  facet_wrap( ~ round) 

#wind roses for each site
head(w)
#speed over time
ggplot(w, aes(x = min_since_start, y = speed)) + geom_point() + facet_wrap( ~ round)

#direction 
ggplot(w, aes(x = direction, y = speed)) + geom_point(alpha = 0.2, size = 3) + facet_wrap( ~ round) + 
  coord_polar() + theme_bw()

#pollen collectors displayed spatially
ggplot(s_summary, aes(x = heading, y = distance_measured, color = log10(p_per_m3))) + 
  geom_point(size = 4) + facet_wrap( ~ round) + 
  coord_polar() + theme_bw() +
  xlab("distance (m)") + ylab("distance (m)") +
  scale_color_continuous(na.value = "gray80", high = "gray10", low = "gray80",
                           name = expression(paste("log pollen per m"^"3"))) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

#wind rose
install.packages("clifro")
library(clifro)
windrose(speed = w$speed, direction = w$direction, facet = as.factor(w$round), 
         speed_cuts = c(0,0.5,1,1.5,2,2.5,4), n_directions = 20,
         ggtheme = "bw", col_pal = "YlOrRd", n_col = 3, legend_title = "Wind Speed\n(m/s)")


#pollen collectors AND wind
ggplot(data = subset(s_summary, heading != 82), aes(x = heading, y = distance_measured)) + 
  geom_point(aes(color = log(p_per_m3)), size = 5)+
  geom_point(data = w, aes(x = direction, y = speed * 5), alpha = 0.5)+
   coord_polar() + theme_bw() +  facet_wrap( ~ round)  





#pollen grains as a function of interior vs. exterior of slide
ggplot(p, aes(x = as.factor(transect), y = p_per_mm2)) + geom_boxplot() 
ggplot(p, aes(x = as.factor(transect), y = p_per_mm2, color = as.factor(sampler))) + geom_jitter(size = 5) + theme_bw() 

