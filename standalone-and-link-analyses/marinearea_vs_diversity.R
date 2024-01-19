run.nr=1

library(layeranalyzer)


# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Read marina area
load("marinearea.Rdata")

# Marine area structure:
marineareast<-layer.series.structure(marinearea,numlayers=1,
  lin.time=T,prior=pr)





# Read diversification

load("pradel_results.RData")
load("pradel_stdev.RData")
load("stages.Rdata")

origts <- layer.data.series(time.points=(-stages$bottom[6:95]),
                            value.points=pradel_results[[1]][3:92],
                            std.dev=pradel_stdev[[1]][3:92],
                            name="origts")

extts <- layer.data.series(time.points=(-stages$top[4:93]),
                           value.points=pradel_results[[2]][1:90],
                           std.dev=pradel_stdev[[2]][1:90],
                           name="extts")

sampts <- layer.data.series(time.points=(-stages$mid[5:94]),
                            value.points=pradel_results[[3]][2:91],
                            std.dev=pradel_stdev[[3]][2:91],
                            name="sampts")


# Diversification structure: 2-layered linear time trend:

origst<-layer.series.structure(origts,numlayers=2,
  lin.time=T,prior=pr)

extst<-layer.series.structure(extts,numlayers=2,
  lin.time=T,prior=pr)

sampst<-layer.series.structure(sampts,numlayers=2,
  lin.time=T,prior=pr)



#####################################
# Marine area vs origination:
#####################################

t1=Sys.time()
mo0=layer.analyzer(marineareast, origst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)
t2=Sys.time(); t2-t1
ml2=mo0$ML


mo1caus=layer.analyzer(marineareast, origst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


mo2caus=layer.analyzer(marineareast, origst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


mo1corr=layer.analyzer(marineareast, origst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


mo2corr=layer.analyzer(marineareast, origst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)





compare.layered(mo0, mo1caus, mo2caus, mo1corr, mo2corr,
   ML.IC="AICc")     
# Majority result (3/4)  
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         922.0937        0.87807
#Model   2         925.3092       21.87770
#Model   3         925.1734       19.09975
#Model   4         926.1476       50.59764
#Model   5         924.2448        7.54684
# Correlation to first layer of diversification

# Minority result (1/4):
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         922.4785        1.84937
#Model   2         924.7384       17.72121
#Model   3         924.9993       23.00357
#Model   4         925.6787       45.37728
#Model   5         924.3526       12.04858
# Also correlation to first layer of diversification



summary(mo1corr)
#Coefficients:
#                              ML estimate Bayesian Lower 95%
#mu_log10_surface_area              6.565165        1.966001
#lin_t_log10_surface_area          -0.000055       -0.002377
#dt_log10_surface_area_1          120.356412       61.995573
#sigma_log10_surface_area_1         0.035384        0.033412
#obs_sd_log10_surface_area          0.000123        0.000001
#mu_origts                         -3.386582       -3.555130
#lin_t_origts                      -0.002073       -0.003193
#dt_origts_1                        0.026632        0.016460
#sigma_origts_1                     5.819578        0.000001
#dt_origts_2                        0.662095        0.099654
#sigma_origts_2                     0.000001        0.000000
#corr_log10_surface_area,1_origts,1 0.947717       -0.649947
#                                   Bayesian Upper 95%
#mu_log10_surface_area                        7.380726
#lin_t_log10_surface_area                     0.001642
#dt_log10_surface_area_1                   2766.763729
#sigma_log10_surface_area_1                   0.037741
#obs_sd_log10_surface_area                    0.004421
#mu_origts                                   -3.234436
#lin_t_origts                                -0.001146
#dt_origts_1                                  2.558845
#sigma_origts_1                               7.664866
#dt_origts_2                                 16.613116
#sigma_origts_2                               1.315480
#corr_log10_surface_area,1_origts,1           0.983931
#
#Log-likelihood:   938.274
#AIC : -1852.548
#AICc: -1852.295
#BIC : -1799.180


save(mo1corr, file="origbest.Rdata")



# Hypothesis tests:

anova(mo0, mo1caus)$Pr[2]
# 0.01447855

anova(mo0, mo1caus)$Pr[2]
# 0.01447855

anova(mo0, mo1corr)$Pr[2]
# 0.001445325

anova(mo0, mo2corr)$Pr[2]
# 0.01179452

save.image(sprintf("marinearea_origdone_run%d", run.nr))


#####################################
# Marine area vs extinction:
#####################################

me0=layer.analyzer(marineareast, extst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


me1caus=layer.analyzer(marineareast, extst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


me2caus=layer.analyzer(marineareast, extst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


me1corr=layer.analyzer(marineareast, extst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


me2corr=layer.analyzer(marineareast, extst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)





compare.layered(me0, me1caus, me2caus, me1corr, me2corr,
   ML.IC="AICc") 
# Majority result (3/4)
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         902.5115        3.07655
#Model   2         902.4701        2.95179
#Model   3         902.6385        3.49296
#Model   4         904.0329       14.08618
#Model   5         905.7236       76.39252
# Correlation to second layer of diversification


# Minorty result (1/4)
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         902.2753        2.38719
#Model   2         902.5237        3.06021
#Model   3         902.1797        2.16961
#Model   4         904.0966       14.75262
#Model   5         905.7572       77.63038
# Also correlation to second layer of diversification

summary(me2corr)
#Coefficients:
#                         ML estimate Bayesian Lower 95%
#mu_log10_surface_area            6.349604       4.064909
#lin_t_log10_surface_area        -0.000253      -0.002075
#dt_log10_surface_area_1        181.007560      58.978968
#sigma_log10_surface_area_1       0.035173       0.033507
#obs_sd_log10_surface_area        0.001417       0.000000
#mu_extts                        -3.622451      -3.775371
#lin_t_extts                     -0.002313      -0.003410
#dt_extts_1                       0.496307       0.057866
#sigma_extts_1                    1.488298       0.000000
#dt_extts_2                       0.624240       0.370283
#sigma_extts_2                    0.654941       0.000818
#corr_log10_surface_area,1_extts,2 -0.812276    -0.943510
#                                  Bayesian Upper 95%
#mu_log10_surface_area                       7.478777
#lin_t_log10_surface_area                    0.001738
#dt_log10_surface_area_1                  1975.311290
#sigma_log10_surface_area_1                  0.037609
#obs_sd_log10_surface_area                   0.004460
#mu_extts                                   -3.355616
#lin_t_extts                                -0.000899
#dt_extts_1                                  0.956922
#sigma_extts_1                               2.699421
#dt_extts_2                                  8.577298
#sigma_extts_2                               2.658534
#corr_log10_surface_area,1_extts,2           0.108723
#
#Log-likelihood:   917.883
#AIC : -1811.767
#AICc: -1811.514
#BIC : -1758.399



save(me2corr, file="extbest.Rdata")



anova(me0, me1caus)$Pr[2]
# 0.3759899

anova(me0, me1caus)$Pr[2]
# 0.3759899

anova(me0, me1corr)$Pr[2]
# 0.02417614

anova(me0, me2corr)$Pr[2]
# 0.003623773

save.image(sprintf("marinearea_extdone_run%d", run.nr))


#####################################
# Marine area vs sampling:
#####################################

ms0=layer.analyzer(marineareast, sampst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ms1caus=layer.analyzer(marineareast, sampst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ms2caus=layer.analyzer(marineareast, sampst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ms1corr=layer.analyzer(marineareast, sampst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ms2corr=layer.analyzer(marineareast, sampst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)





compare.layered(ms0, ms1caus, ms2caus, ms1corr, ms2corr,
   ML.IC="AICc")        
# Majority result (3/4):
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         949.5476       15.97950
#Model   2         950.3497       35.63890
#Model   3         950.0027       25.18998
#Model   4         949.3980       13.75866
#Model   5         949.0205        9.43295
# causal link to first layer of sampling

# minority result (1/4)
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         949.6600       18.67612
#Model   2         950.2699       34.36700
#Model   3         950.0130       26.58220
#Model   4         948.9850        9.50936
#Model   5         949.1183       10.86532
# Also causal link to first layer of sampling




summary(ms1caus)
#Coefficients:
#                             ML estimate Bayesian Lower 95%
#mu_log10_surface_area             6.365310       4.917617
#lin_t_log10_surface_area         -0.005969      -0.009727
#dt_log10_surface_area_1         149.680171      60.100724
#sigma_log10_surface_area_1        0.035417       0.033522
#obs_sd_log10_surface_area         0.002208       0.000001
#mu_sampts                        -1.886929      -2.530402
#lin_t_sampts                     -0.000730      -0.002774
#dt_sampts_1                       0.693684       0.098455
#sigma_sampts_1                    0.673507       0.359763
#dt_sampts_2                      24.964481       7.775981
#sigma_sampts_2                    0.152940       0.004975
#beta_log10_surface_area,1_sampts,1-0.847055     -1.114175
#                                   Bayesian Upper 95%
#mu_log10_surface_area                        6.855242
#lin_t_log10_surface_area                     0.009844
#dt_log10_surface_area_1                   2805.437731
#sigma_log10_surface_area_1                   0.037628
#obs_sd_log10_surface_area                    0.004722
#mu_sampts                                   -0.865482
#lin_t_sampts                                 0.002496
#dt_sampts_1                                  5.350627
#sigma_sampts_1                               1.860264
#dt_sampts_2                                140.867024
#sigma_sampts_2                               0.248754
#beta_log10_surface_area,1_sampts,1           0.311999
#
#Log-likelihood:   962.476
#AIC : -1900.952
#AICc: -1900.699
#BIC : -1847.584


save(ms1caus, file="sampbest.Rdata")



anova(ms0, ms1caus)$Pr[2]
# 0.1617471

anova(ms0, ms2caus)$Pr[2]
# 0.2288407

anova(ms0, ms1corr)$Pr[2]
# 0.187151

anova(ms0, ms2corr)$Pr[2]
# 0.3209712



save.image(sprintf("marinearea_sampdone_run%d", run.nr))


