run.nr=1

library(layeranalyzer)


# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Read sealevel
load("sealevel.Rdata")


# Sealevel structure:
seast<-layer.series.structure(sealevel,numlayers=2,
  no.sigma=c(1),
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
# Sealevel vs origination:
#####################################

t1=Sys.time()
so0=layer.analyzer(seast, origst,
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)
ml2=so0$ML

t2=Sys.time(); t2-t1



s1o1caus=layer.analyzer(seast, origst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s1o2caus=layer.analyzer(seast, origst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)

s2o1caus=layer.analyzer(seast, origst,
       causal=cbind(c(1,2,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s2o2caus=layer.analyzer(seast, origst,
       causal=cbind(c(1,2,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s2o1corr=layer.analyzer(seast, origst,
       corr=cbind(c(1,2,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s2o2corr=layer.analyzer(seast, origst,
       corr=cbind(c(1,2,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)



compare.layered(so0, s1o1caus, s1o2caus, 
   s2o1caus, s2o2caus, s2o1corr, s2o2corr,
   ML.IC="AICc")     
#  Windows
#                   weight=-0.5*AICc Post. Prob.(%)
#Model   1         545.5422       17.78740
#Model   2         544.3317        5.30161
#Model   3         544.2434        4.85352
#Model   4         545.5652       18.20101
#Model   5         545.4168       15.69084
#Model   6         544.9237        9.58318
#Model   7         546.0165       28.58244



# Hypothesis tests:

anova(so0, s1o1caus)$Pr[2]
# 1 W

anova(so0, s1o2caus)$Pr[2]
# 1 W

anova(so0, s2o1caus)$Pr[2]
# 0.3525452 W

anova(so0, s2o2caus)$Pr[2]
# 0.4089442 W

anova(so0, s2o1corr)$Pr[2]
# 0.3704312 W

anova(so0, s2o2corr)$Pr[2]
#  0.08389488


save.image(sprintf("sealevel_vs_orig%d.Rdata",run.nr))





#####################################
# Sealevel vs extinction:
#####################################

se0=layer.analyzer(seast, extst,
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s1e1caus=layer.analyzer(seast, extst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s1e2caus=layer.analyzer(seast, extst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)

s2e1caus=layer.analyzer(seast, extst,
       causal=cbind(c(1,2,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)

s2e2caus=layer.analyzer(seast, extst,
       causal=cbind(c(1,2,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)

s2e1corr=layer.analyzer(seast, extst,
       corr=cbind(c(1,2,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)

s2e2corr=layer.analyzer(seast, extst,
       corr=cbind(c(1,2,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)






compare.layered(se0, s1e1caus, s1e2caus, 
   s2e1caus, s2e2caus, s2e1corr, s2e2corr,
   ML.IC="AICc") 
# 

summary(s2e2corr)
#


anova(se0, s1e1caus)$Pr[2]
# 

anova(se0, s1e1caus)$Pr[2]
# 

anova(se0, s2e1caus)$Pr[2]
# 

anova(se0, s2e1caus)$Pr[2]
# 

anova(se0, s2e1corr)$Pr[2]
# 

anova(se0, s2e2corr)$Pr[2]
# 


save.image(sprintf("sealevel_vs_ext%d.Rdata",run.nr))



#####################################
# Sealevel vs sampling:
#####################################

ss0=layer.analyzer(seast, sampst,
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s1s1caus=layer.analyzer(seast, sampst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s1s2caus=layer.analyzer(seast, sampst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)

s2s1caus=layer.analyzer(seast, sampst,
       causal=cbind(c(1,2,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s2s2caus=layer.analyzer(seast, sampst,
       causal=cbind(c(1,2,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s2s1corr=layer.analyzer(seast, sampst,
       corr=cbind(c(1,2,2,1)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


s2s2corr=layer.analyzer(seast, sampst,
       corr=cbind(c(1,2,2,2)),
       num.MCMC=1000,spacing=5,burnin=4000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=1000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)





compare.layered(ss0, s1s1caus, s1s2caus, 
   s2s1caus, s2s2caus, s2s1corr, s2s2corr,
   ML.IC="AICc")        
#


summary(s2s1caus) # depends on which was actually best
#


save.image(file=sprintf("sealevel_vs_div%d.Rdata",run.nr))


anova(ss0, s1s1caus)$Pr[2]
# 

anova(ss0, s1s2caus)$Pr[2]
# 

anova(ss0, s2s1caus)$Pr[2]
# 

anova(ss0, s2s2caus)$Pr[2]
# 

anova(ss0, s2s1corr)$Pr[2]
# 

anova(ss0, s2s2corr)$Pr[2]
# 


save.image(sprintf("sealevel_vs_samp%d.Rdata",run.nr))




