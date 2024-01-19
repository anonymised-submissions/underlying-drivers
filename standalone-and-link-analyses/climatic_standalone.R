library(layeranalyzer)

# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Marine area:

load("marinearea.Rdata")

res=traverse.standalone.layered(marinearea,
  max.layers=2, talkative=TRUE, 
  allow.one.feedback.loop = TRUE, 
  just.stationary = FALSE, no.rw = FALSE, 
  allow.deterministic.layers = TRUE, 
  num.MCMC = 100, spacing = 10, burnin = 1000, num.temp = 4, 
  do.maximum.likelihood = TRUE, 
  maximum.likelihood.numstart = 200,
  prior = pr) 

res.nolin=list(res[[1]],res[[2]],res[[4]],res[[5]],
  res[[6]],res[[7]],res[[8]],res[[9]])

compare.layered(res.nolin,ML.IC="AICc")
# 2 times:
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         1044.642        0.10223
#Model   2         1043.168        0.02341
#Model   3         1051.488       96.07977
#Model   4         1039.558        0.00063
#Model   5         1041.396        0.00398
#Model   6         1043.105        0.02197
#Model   7         1042.874        0.01744
#Model   8         1046.321        0.54798
#Model   9         1045.956        0.38020
#Model  10         1041.349        0.00380
#Model  11         1040.649        0.00188
#Model  12         1047.612        1.99259
#Model  13         1046.729        0.82413

# 1 time:
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         1044.174        3.72450
#Model   2         1041.743        0.32743
#Model   3         1046.699       46.50927
#Model   4         1045.413       12.85840
#Model   5         1043.430        1.76937
#Model   6         1043.711        2.34518
#Model   7         1041.150        0.18112
#Model   8         1042.677        0.83338
#Model   9         1041.931        0.39522
#Model  10         1044.236        3.96139
#Model  11         1041.350        0.22105
#Model  12         1045.143        9.81639
#Model  13         1045.696       17.05730

# Anyway, 1-layered OU with linear trend

# Take out best:
ma1=c(1044.642, 1043.743, 1051.488, 1045.413, 1043.430, 1043.711, 1042.874,
  1046.321, 1045.956, 1044.236, 1041.350, 1047.612, 1046.729)
a1=-2*ma1; a1
exp((ma1-1044))/sum(exp((ma1-1044)))*100
# [1]  0.101845274  0.041448627 95.746058797  0.220182059  0.030309305
# [6]  0.040143268  0.017382340  0.545910337  0.378969128  0.067860542
#[11]  0.003786548  1.985159902  0.820943875







summary(res[[3]])


save(res,file="marea_standalone.Rdata")




# Recreate:
marea.st=layer.series.structure(marinearea,no.pull=F,
 lin.time=T,prior=pr,init.0=T)
marbest=layer.analyzer(marea.st, 
      num.MCMC=10000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,maximum.likelihood.numstart=5000,
      silent.mode=FALSE,talkative.burnin=TRUE,talkative.likelihood=FALSE)
summary(marbest)
#Coefficients:
#                                 ML estimate Bayesian Lower 95%
#mu_log10_marinearea_area            6.611600           4.709222
#lin_t_log10_marinearea_area        -0.001192          -0.007457
#dt_log10_marinearea_area_1         41.828936          26.755578
#sigma_log10_marinearea_area_1       0.034590           0.033229
#obs_sd_log10_marinearea_area        0.000000           0.000000
#init_log10_marinearea_area_l1_s0    6.092191           6.088516
#                                 Bayesian Upper 95%
#mu_log10_marinearea_area                   7.078587
#lin_t_log10_marinearea_area               -0.000242
#dt_log10_marinearea_area_1              2120.707972
#sigma_log10_marinearea_area_1              0.037560
#obs_sd_log10_marinearea_area               0.004970
#init_log10_marinearea_area_l1_s0           6.095691
#
#Log-likelihood:  1058.925
#AIC : -2105.850
#AICc: -2105.771
#BIC : -2080.090




save(marbest,file="marinearea_best.Rdata")








# sealevel:

load("sealevel.Rdata")

res2=traverse.standalone.layered(sealevel,
  max.layers=2, talkative=TRUE, 
  allow.one.feedback.loop = TRUE, 
  just.stationary = FALSE, no.rw = FALSE, 
  allow.deterministic.layers = TRUE, 
  num.MCMC = 1000, spacing = 10, burnin = 10000, num.temp = 4, 
  do.maximum.likelihood = TRUE, 
  maximum.likelihood.numstart = 1000,
  prior = pr) 

compare.layered(res2,ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -1875.189        0.11592
#Model   2        -1873.937        0.40540
#Model   3        -1877.414        0.01253
#Model   4        -1878.131        0.00612
#Model   5        -1880.315        0.00069
#Model   6        -1880.642        0.00050
#Model   7        -1877.741        0.00903
#Model   8        -1876.814        0.02283
#Model   9        -1868.435       99.41928
#Model  10        -1881.195        0.00029
#Model  11        -1879.618        0.00138
#Model  12        -1879.645        0.00135
#Model  13        -1878.397        0.00469

# 3 times same result?

ma2=c(-1875.189, -1873.937, -1877.414, -1878.131, -1880.315, -1880.642,
  -1877.741, -1876.814, -1868.435, -1881.195, -1879.618, -1879.645, -1878.397)
a2=-2*ma2; a2





summary(res2[[9]])
#Coefficients:
#                    ML estimate Bayesian Lower 95% Bayesian Upper 95%
#mu_sealevel           -6.376800          -9.721924           9.363685
#dt_sealevel_1          2.881514           0.228070           3.429604
#sigma_sealevel_2       4.819201           3.534162           6.273011
#obs_sd_sealevel        9.651841           0.000062          22.321637
#init_sealevel_l1_s0   49.483677           7.020318          27.350942
#init_sealevel_l2_s0   75.989416          10.675213          28.712325
#
#Log-likelihood: -1862.396
#AIC :  3736.791
#AICc:  3736.870
#BIC :  3762.552

# 2 layers, RW at bottom, no feedback, 
# deterministic first layer. 
# (Prime interpretation:
# first layer is caused by post-process
# smoothing or interpolation).)

save.image("sealevel_standalone1.Rdata")



# global temperature

load("global_temperature.Rdata")

res3=traverse.standalone.layered(global_temperature,
  max.layers=2, talkative=TRUE, 
  allow.one.feedback.loop = TRUE, 
  just.stationary = FALSE, no.rw = FALSE, 
  allow.deterministic.layers = TRUE, 
  num.MCMC = 100, spacing = 10, burnin = 1000, num.temp = 4, 
  do.maximum.likelihood = TRUE, 
  maximum.likelihood.numstart = 200,
  prior = pr) 

compare.layered(res3,ML.IC="AICc")

# Large run (num.MCMC=1000, burnin=10000,
#   maximum.likelihood.numstart = 2400)
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -297.5925       65.28006
#Model   2        -300.0422        5.63488
#Model   3        -299.2197       12.82601
#Model   4        -303.3325        0.20987
#Model   5        -303.8154        0.12948
#Model   6        -300.5673        3.33288
#Model   7        -300.2892        4.40154
#Model   8        -300.1493        5.06265
#Model   9        -301.5264        1.27735
#Model  10        -302.6266        0.42512
#Model  11        -303.1526        0.25121
#Model  12        -303.8532        0.12468
#Model  13        -301.7279        1.04427
# OU, simply

# Small run, (num.MCMC=100, burnin=1000,
#   maximum.likelihood.numstart = 200)
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -297.9746       73.49246
#Model   2        -300.0426        9.29293
#Model   3        -301.0262        3.47497
#Model   4        -303.5509        0.27829
#Model   5        -304.4715        0.11084
#Model   6        -306.8245        0.01054
#Model   7        -300.2582        7.49016
#Model   8        -301.0471        3.40319
#Model   9        -301.8683        1.49707
#Model  10        -303.7592        0.22596
#Model  11        -303.8548        0.20535
#Model  12        -303.1482        0.41630
#Model  13        -304.5551        0.10194

#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -298.8595       60.45101
#Model   2        -300.0424       18.52105
#Model   3        -300.9718        7.31247
#Model   4        -303.4263        0.62820
#Model   5        -303.9952        0.35563
#Model   6        -305.5427        0.07567
#Model   7        -301.5291        4.18832
#Model   8        -302.1554        2.23877
#Model   9        -301.8764        2.95925
#Model  10        -302.9917        0.97011
#Model  11        -304.0788        0.32712
#Model  12        -306.9369        0.01877
#Model  13        -302.2917        1.95364

#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -297.8648       79.93649
#Model   2        -300.0439        9.04357
#Model   3        -301.0363        3.35242
#Model   4        -303.2796        0.35571
#Model   5        -303.2254        0.37554
#Model   6        -306.3092        0.01719
#Model   7        -300.9613        3.61343
#Model   8        -302.1087        1.14719
#Model   9        -302.4607        0.80678
#Model  10        -302.8990        0.52048
#Model  11        -303.9692        0.17849
#Model  12        -306.9818        0.00878
#Model  13        -302.6861        0.64394

# There's  a bit of variation, but the main outcome 
# does not chance: single layer OU is the best.

ma3=c(-297.5925, -300.0422, -299.2197, -301.0262, -303.2254,
  -300.5673, -300.2582, -300.1493, -301.5264, -302.6266,
  -303.1526, -303.1482, -301.7279)
a3=-2*ma3; a3

exp((ma3+300))/sum(exp((ma3+300)))*100
# [1] 63.8327467  5.5100094 12.5417908  2.0597123  0.2284053  3.2591497
# [7]  4.4396091  4.9503898  1.2490289  0.4156824  0.2456534  0.2467367
#[13]  1.0210856


summary(res3[[1]])

save(res3,file="temp_standalone.Rdata")


# Recreate:
temp.st=layer.series.structure(global_temperature,numlayers=1,prior=pr,init.0=T)
tempbest=layer.analyzer(temp.st, 
      num.MCMC=10000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,maximum.likelihood.numstart=5000,
      silent.mode=FALSE,talkative.burnin=TRUE,talkative.likelihood=FALSE)
summary(tempbest)
#Coefficients:
#                                    ML estimate Bayesian Lower 95%
#mu_global_temperature.Rdata           19.356206          -8.198254
#dt_global_temperature.Rdata_1         34.266699          92.346286
#sigma_global_temperature.Rdata_1       0.897930           0.581793
#obs_sd_global_temperature.Rdata        1.992332           1.722999
#init_global_temperature.Rdata_l1_s0   26.699873          20.057173
#                                    Bayesian Upper 95%
#mu_global_temperature.Rdata                  14.654750
#dt_global_temperature.Rdata_1              3978.122842
#sigma_global_temperature.Rdata_1              1.056634
#obs_sd_global_temperature.Rdata               2.594148
#init_global_temperature.Rdata_l1_s0          27.146174
#
#Log-likelihood:  -292.385
#AIC :   594.770
#AICc:   595.041
#BIC :   608.581


save(tempbest, file="temp_best.Rdata")

