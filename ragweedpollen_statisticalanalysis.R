#pollen statistical analysis
#this is the statistical analysis of the data visualized in ragweedpollenanalysis.R

#set up work environment

library(ggplot2)
library(dplyr)
library(rjags)

##########
rm(list = ls())


##########
#load in data
##########

#load in pollen count data
p <- read.csv("Pollen Grain Counting Record_151019.csv")
#p <- subset(p, !is.na(p$grains_dk))  #taking out the unmeasured stuff for now
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
w_summary <- tbl_df(w) %>% filter(ss == 0) %>% group_by(round) %>%
  summarise(direction_mean = mean(direction), direction_sd = sd(direction),
            speed_mean = mean(speed), speed_sd = sd(speed))



###########
#statistical analysis
###########






niterations <- 100


#############
#log normal hurdle model 
#############

sink("model.txt")
cat("   
    
    model{  
    
    #######model discrete data
    for(i in 1:nsamplers){
    y[i] ~ dnorm(mu[i], tau) 

    mu <- exp(beta_0 + beta_1 * dist[i]) + beta_2*w[i]
    w[i] <- 
    }

    tau <- 1/(sigma*sigma)
    sigma ~ dunif( 0.01, 10)



    
    logit(probzero[i]) <- #kappa_conba10_d*conspecific_ba_d[i] + kappa_relbasite_d*conspecific_relba_site_d[i] +
    #    kappa_congen10_d *  congeneric_ba_d[i] +
    multiscaled_intercept_d[multiscaled_d[i]] + #conspecific_dist_slope_d * conspecific_dist_d[i] +
    delta_prep_d*prepheight_s_d[i] + delta_age_days_d*age_days_d[i] + delta_nleaves_d * nleaves_d[i] +
    beta_gsf_d * gsf_d[i] + beta_sm_d * sm_d[i] +     
    #alpha_d[region_d[i]] +
    gamma_seedling_d[seedling_d[i]]  + gamma_year_d[year_d[i]] 
    #+ gamma_plot_d[plot_d[i]] + gamma_plot_year_d[plot_d[i],year_d[i]]
    }
    
    
    #           kappa_conba10_d ~ dnorm(0, 0.0001)
    #           kappa_relbasite_d ~ dnorm(0, 0.0001)
    #           kappa_congen10_d ~ dnorm(0, 0.0001)
    delta_prep_d ~ dnorm(0, 0.0001)
    delta_age_days_d ~ dnorm(0, 0.0001)
    delta_nleaves_d ~ dnorm(0, 0.0001)
    
    for(i in 1:nmultiscaled_d){
    multiscaled_intercept_d[i]~ dnorm(0, 0.0001)
    }
    
    #          conspecific_dist_slope_d ~ dnorm(0, 0.0001)
    
    beta_gsf_d ~ dnorm(0, 0.0001)
    beta_sm_d ~ dnorm(0, 0.0001)
    
    # #             for(i in 1:nregion_d){
    # #               beta_gsf_d[i] ~ dnorm(0, 0.0001)
    # #               beta_sm_d[i] ~ dnorm(0, 0.0001)
    #            #   alpha_d[i] ~ dnorm(0, 0.0001)
    #             }
    
    ## for(i in 1:nmigstat_d){
    #  alpha_d[i] ~ dnorm(0, 0.0001)
    #}
    
    #  for(i in 1:nplot_d){
    #          gamma_plot_d[i] ~ dnorm(0, tau.plot_d)
    ##          tau.plot_year_d[i] <- pow(sigma.plot_year_d[i], -2)
    #          sigma.plot_year_d[i] ~ dunif(0,10)
    
    #      for(j in 1:nyear_d){
    #          gamma_plot_year_d[i,j] ~ dnorm(gamma_plot_d[i], tau.plot_year_d[i])
    #      }                
    #   }
    #   tau.plot_d <- pow(sigma.plot_d, -2)
    #   sigma.plot_d ~ dunif(0,10)
    
    for(i in 1:nyear_d){
    gamma_year_d[i] ~ dnorm(0, tau.year_d)
    }
    tau.year_d <- pow(sigma.year_d, -2)
    sigma.year_d ~ dunif(0,10)
    
    #    for(i in 1:nplot_d){
    #          gamma_plot_d[i] ~ dnorm(0, tau.plot_d)
    #      }
    #      tau.plot_d <- pow(sigma.plot_d, -2)
    #      sigma.plot_d ~ dunif(0,10)
    
    
    for(i in 1:nseedling_d){
    gamma_seedling_d[i] ~ dnorm(0, tau.seedling_d)
    }
    tau.seedling_d <- pow(sigma.seedling_d, -2)
    sigma.seedling_d ~ dunif(0,10)
    
    ##multiscaled comparison
    for(i in 1:nmultiscaled_d){
    multiscaled_dif_d[i] <- multiscaled_intercept_d[1] - multiscaled_intercept_d[i]
    }
    
    #significance between multiscaled
    for(x in 1:nmultiscaled_d){
    for(y in 1:nmultiscaled_d){
    m_sig_d[x,y] <-  multiscaled_intercept_d[x] - multiscaled_intercept_d[y] 
    }
    }
    
    
    
    
    #########log normal model of continuous data
    for(i in 1:nobs_c){
    y_c[i] ~ dlnorm( meanval[i], tau.total_c ) 
    
    meanval[i] <- multiscaled_intercept_c[multiscaled_c[i]] + #conspecific_dist_slope_c * conspecific_dist_c[i] +
    #                         alpha_c[region_c[i]] +                        
    #                         kappa_conba10_c * conspecific_ba_c[i] + kappa_relbasite_c * conspecific_relba_site_c[i] + kappa_congen10_c * 
    #                             congeneric_ba_c[i] +
    delta_prep_c * prepheight_s_c[i] + delta_age_days_c * age_days_c[i] + delta_nleaves_c * nleaves_c[i] +
    beta_gsf_c * gsf_c[i] + beta_sm_c * sm_c[i] + 
    gamma_seedling_c[seedling_c[i]]  + gamma_year_c[year_c[i]]   #+ gamma_plot_c[plot_c[i]]
    # gamma_plot_year_c[plot_c[i],year_c[i]]  +
    
    
    #predicted values for when herbiv occurs
    y_non_zeroes_pred[i] ~ dlnorm( meanval[i], tau.total_c)
    }
    
    #priors
    #for(i in 1:nmigstat_c){
    #  alpha_c[i] ~ dnorm(0, 0.0001)
    #}
    
    for(i in 1:nmultiscaled_c){
    multiscaled_intercept_c[i]~ dnorm(0, 0.0001)
    }
    
    #          conspecific_dist_slope_c ~ dnorm(0, 0.0001)
    #             kappa_relbasite_c ~ dnorm(0, 0.0001)
    #             kappa_conba10_c ~ dnorm(0, 0.0001)
    #             kappa_congen10_c ~ dnorm(0, 0.0001)
    
    delta_prep_c ~ dnorm(0, 0.0001)
    delta_age_days_c ~ dnorm(0, 0.0001)
    delta_nleaves_c ~ dnorm(0, 0.0001)
    
    beta_gsf_c ~ dnorm(0, 0.0001)
    beta_sm_c ~ dnorm(0, 0.0001)
    
    #             for(i in 1:nregion_c){
    #               beta_gsf_c[i] ~ dnorm(0, 0.0001)
    #               beta_sm_c[i] ~ dnorm(0, 0.0001)
    # #               alpha_c[i] ~ dnorm(0, 0.0001)
    #             }
    
    for(i in 1:nyear_c){
    gamma_year_c[i] ~ dnorm(0, tau.year_c)
    }
    tau.year_c <- pow(sigma.year_c, -2)
    sigma.year_c ~ dunif(0,10)
    
    
    #    for(i in 1:nplot_c){
    #        gamma_plot_c[i] ~ dnorm(0, tau.plot_c)
    #    }
    #    tau.plot_c <- pow(sigma.plot_c, -2)
    #    sigma.plot_c ~ dunif(0,10)
    
    for(i in 1:nseedling_c){
    gamma_seedling_c[i] ~ dnorm(0, tau.seedling_c)
    }
    tau.seedling_c <- pow(sigma.seedling_c, -2)
    sigma.seedling_c ~ dunif(0,10)
    
    tau.total_c <- pow(sigma_c, -2)
    sigma_c ~ dunif(0,10)
    
    #multiscaled comparison
    for(i in 1:nmultiscaled_c){
    multiscaled_dif_c[i] <- multiscaled_intercept_c[1] - multiscaled_intercept_c[i]
    }
    
    #significance between multiscaled
    for(x in 1:nmultiscaled_c){
    for(y in 1:nmultiscaled_c){
    m_sig_c[x,y] <-  multiscaled_intercept_c[x] - multiscaled_intercept_c[y] 
    }
    }
    
    }
    ",fill=TRUE)
sink() 

jags <- jags.model('model.txt',
                   data = list(
                     
                     #########  discrete     
                     #indices and response vars
                     'nobs_d' = nrow(survexp_sub),
                     #                      'region_d' = survexp_sub$region_seq,
                     #                     #'site_d' = survexp_sub$site_seq,
                     #                      'nregion_d' = max(survexp_sub$region_seq),
                     #'nsite_d' = max(survexp_sub$site_seq),
                     # 'plot_d' = survexp_sub$plot_seq,
                     # 'nplot_d' = max(survexp_sub$plot_seq),
                     'seedling_d' = survexp_sub$seedling,
                     'nseedling_d' = max(survexp_sub$seedling),
                     #  'migstat_d' = survexp_sub$migstat_seq,   
                     #  'nmigstat_d' = max(survexp_sub$migstat_seq), 
                     'year_d' = survexp_sub$year_seq,
                     'nyear_d' = max(survexp_sub$year_seq),
                     
                     #seedling data
                     'y_d' = ifelse(survexp_sub$p > 0, 1, 0),
                     'prepheight_s_d' = survexp_sub$prepheight_s,
                     'age_days_d' = survexp_sub$age_days_s, 
                     'nleaves_d'= survexp_sub$nleaves_s,
                     
                     #environmental variables
                     'gsf_mean_d' = survexp_sub$gsf_mean,
                     'gsf_sd_d' = survexp_sub$gsf_sd,
                     'sm_mean_d' = survexp_sub$sm_mean,
                     'sm_sd_d' = survexp_sub$sm_sd,
                     'sm_totalmean_mean_d' = sm_totalmean_mean,
                     'sm_totalmean_sd_d' = sm_totalmean_sd,
                     'sm_totalsd_mean_d' = sm_totalsd_mean,
                     'sm_totalsd_sd_d' = sm_totalsd_sd,
                     
                     #covariates
                     'multiscaled_d' = survexp_sub$multiscaled2,
                     'nmultiscaled_d' = max(survexp_sub$multiscaled2),
                     #                   'conspecific_dist_d' = survexp_sub$conspecific_d,
                     #                      'conspecific_ba_d' = survexp_sub$conspecific_d,  #conspecific BA within 10m
                     #                                         #  'conspecific_relba_site_d' = survexp_sub$conspecific_relba_site,  #conspecific BA @ site                   
                     #                      'conspecific_relba_site_d' =      #conspecific BA @ site - only present at enough sites for acru, quru 
                     #                        if(sp_subset == "acru" | sp_subset == "quru"){survexp_sub$conspecific_relba_site}else{
                     #                        survexp_sub$conspecific_relba_site * 0},                   
                     #                      'congeneric_ba_d' = survexp_sub$congeneric_ba,  #congeneric BA within 10m
                     #                      'conspecific_ba_horner_mean_d' = 
                     #                         ifelse(sp_subset == "beth" |sp_subset == "ceor" |sp_subset == "elum" |sp_subset == "cagl", 
                     #                                0,  #no invasives planted at horner
                     #                              mean(survexp_sub$conspecific_ba[survexp_sub$site== "horner"], na.rm =TRUE)),  
                     #                      'conspecific_ba_horner_sd_d' = 
                     #                         ifelse(sp_subset == "beth" |sp_subset == "ceor" |sp_subset == "elum"|sp_subset == "cagl", 1000,  #no invasives planted at horner
                     #                              sd(survexp_sub$conspecific_ba[survexp_sub$site== "horner"], na.rm =TRUE)),  
                     
                     #########  continuous     
                     #indices and response vars
                     'nobs_c' = nrow(survexp_sub_cont),                    
                     #                        'region_c' = survexp_sub_cont$region_seq,
                     #                        #'site_c' = survexp_sub_cont$site_seq,
                     #                        'nregion_c' = max(survexp_sub_cont$region_seq),
                     #                        #'nsite_c' = max(survexp_sub_cont$site_seq),
                     #'plot_c' = survexp_sub_cont$plot_seq,
                     # 'nplot_c' = max(survexp_sub_cont$plot_seq),
                     'seedling_c' = survexp_sub_cont$seedling,
                     'nseedling_c' = max(survexp_sub_cont$seedling),                   
                     # 'migstat_c' = survexp_sub_cont$migstat_seq,   
                     # 'nmigstat_c' = max(survexp_sub_cont$migstat_seq), 
                     'year_c' = survexp_sub_cont$year_seq,
                     'nyear_c' = max(survexp_sub_cont$year_seq),
                     
                     #seedling data
                     'y_c' = survexp_sub_cont$p,
                     'prepheight_s_c' = survexp_sub_cont$prepheight_s,
                     'age_days_c' = survexp_sub_cont$age_days_s,   
                     'nleaves_c' = survexp_sub_cont$nleaves_s,
                     
                     #environmental variables
                     'gsf_mean_c' = survexp_sub_cont$gsf_mean,
                     'gsf_sd_c' = survexp_sub_cont$gsf_sd,
                     'sm_mean_c' = survexp_sub_cont$sm_mean,
                     'sm_sd_c' = survexp_sub_cont$sm_sd,
                     'sm_totalmean_mean_c' = sm_totalmean_mean,
                     'sm_totalmean_sd_c' = sm_totalmean_sd,
                     'sm_totalsd_mean_c' = sm_totalsd_mean,
                     'sm_totalsd_sd_c' = sm_totalsd_sd,
                     
                     #covariates
                     'multiscaled_c' = survexp_sub_cont$multiscaled2,
                     'nmultiscaled_c' = max(survexp_sub_cont$multiscaled2)
                     #                   'conspecific_dist_c' = survexp_sub_cont$conspecific_d
                     #                        'conspecific_ba_c' = survexp_sub_cont$conspecific_ba,  #conspecific BA within 10m
                     #                        'conspecific_relba_site_c' =      #conspecific BA @ site - only present at enough sites for acru, quru 
                     #                              if(sp_subset == "acru" | sp_subset == "quru"){survexp_sub_cont$conspecific_relba_site}else{
                     #                                survexp_sub_cont$conspecific_relba_site * 0},                          
                     #                        'congeneric_ba_c' = survexp_sub_cont$congeneric_ba,  #congeneric BA within 10m
                     #                        'conspecific_ba_horner_mean_c' = 
                     #                             ifelse(sp_subset == "beth" |sp_subset == "ceor" |sp_subset == "elum" | sp_subset == "qual"|sp_subset == "cagl"
                     #                                    |sp_subset == "rops"
                     #                                    , 0,  #no invasives planted at horner
                     #                             mean(survexp_sub_cont$conspecific_ba[survexp_sub_cont$site== "horner"], na.rm =TRUE)),  
                     #                        'conspecific_ba_horner_sd_c' = 
                     #                             ifelse(sp_subset == "beth" |sp_subset == "ceor" |sp_subset == "elum" | sp_subset == "cagl"|
                     #                                      sp_subset == "qual"|sp_subset == "rops", 1000,  #no invasives planted at horner
                     #                             sd(survexp_sub_cont$conspecific_ba[survexp_sub_cont$site== "horner"], na.rm =TRUE))  
                     
                   ),
                   n.chains = 3,
                   n.adapt = 100)

update(jags,n.iter= niterations) #update(jags,n.iter=100) 
mcmc_samples <- coda.samples(jags, variable.names=c( #"alpha_d", "alpha_c",
  #"gamma_plot_year_c", "gamma_plot_c", "gamma_plot_year_c",
  #"kappa_conba10_d", "kappa_conba10_c", "kappa_congen10_d", "kappa_congen10_c",
  #"kappa_relbasite_d", "kappa_relbasite_c",
  "multiscaled_intercept_c", "multiscaled_intercept_d",
  "multiscaled_dif_d","m_sig_d",
  "multiscaled_dif_c","m_sig_c",  
  # "conspecific_dist_slope_d", "conspecific_dist_slope_c",
  "delta_prep_d", "delta_prep_c", "delta_age_days_d", "delta_age_days_c", 
  "delta_nleaves_d", "delta_nleaves_c",
  "beta_gsf_d", "beta_gsf_c", "beta_sm_d", "beta_sm_c"      
), n.iter= niterations*3)
#plot(mcmc_samples)  #plot(samples[,'a'])
#gelman.diag(mcmc_samples)
#gelman.plot(mcmc_samples)
# dic.jags <- (dic.samples(jags, 1000, "pD"))
#dicjags2 <- (dic.samples(jags, 1000, "pD"))
#diffdic(dic.jags, dic.jags2)

plotname<-paste(sp_subset,"_chainhistory.pdf",sep="")
pdf(plotname);plot(mcmc_samples); dev.off()    #saving chain history into a pdf

#saving parameter estimates
resultsall<-summary(mcmc_samples)    #takes a while for residuals
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param<-substr(results$parameter,1,2)
results$param2<-substr(results$parameter,1,5)
results$param3<-substr(results$parameter,1,1)
results$index<-1:length(results$Mean)
results$sp <- sp_subset  

results$timestamp <- Sys.time()
resultsfilename <- paste("results_",sp_subset,".csv",sep="")
write.csv(results,resultsfilename)

#residuals
mcmc_samples <- coda.samples(jags, variable.names=c("y_non_zeroes_pred"),  n.iter= niterations)

resultsall<-summary(mcmc_samples)    #takes a while for residuals
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param<-substr(results$parameter,1,2)
results$param2<-substr(results$parameter,1,5)
results$param3<-substr(results$parameter,1,1)
results$index<-1:length(results$Mean)


residplots <- data.frame(results$Mean, results$X2.5., results$X97.5., survexp_sub_cont$p, survexp_sub_cont$ba_allsp)
fit <- summary(lm(residplots$results.Mean ~ residplots$survexp_sub_cont.p))
plottitle <-  paste(sp_subset," R2 =",round(fit$r.squared,3))

ggplot(residplots, aes(x= survexp_sub_cont.p, y = results.Mean)) +
  geom_point() + theme_bw(base_size = 16) + geom_smooth(method = "lm") + ggtitle(as.character(plottitle)) +
  xlab("observed p") + ylab("predicted p") 
#  geom_errorbar(aes(x = survexp_sub_cont.p, ymax = results.X97.5., ymin = results.X2.5., alpha = 0.05))   
plotname <- paste(sp_subset,"_modelfit.pdf",sep="")
ggsave(plotname)

