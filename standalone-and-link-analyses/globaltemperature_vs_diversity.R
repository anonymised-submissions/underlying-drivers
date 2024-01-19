run.nr=1

library(layeranalyzer)


# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Read global temperature
load("global_temperature.Rdata")


# Global temperature structure:
tempst<-layer.series.structure(global_temperature,numlayers=1,
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
# Global temperature vs origination:
#####################################

t1=Sys.time()
to0=layer.analyzer(tempst, origst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=30000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)
ml2=to0$ML


to1caus=layer.analyzer(tempst, origst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=30000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


to2caus=layer.analyzer(tempst, origst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=30000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


to1corr=layer.analyzer(tempst, origst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=30000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


to2corr=layer.analyzer(tempst, origst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=30000,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)
t2=Sys.time(); t2-t1




compare.layered(to0, to1caus, to2caus, to1corr, to2corr,
   ML.IC="AICc")
#

summary(to1corr)
#

# Hypothesis tests:

anova(to0, to1caus)$Pr[2]
# 

anova(to0, to2caus)$Pr[2]
# 

anova(to0, to1corr)$Pr[2]
#  

anova(to0, to2corr)$Pr[2]
#  



save.image(file=sprintf("globaltemp_origdone_run%1d.Rdata",run.nr))



#####################################
# Global temeprature vs extinction:
#####################################

te0=layer.analyzer(tempst, extst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


te1caus=layer.analyzer(tempst, extst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


te2caus=layer.analyzer(tempst, extst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


te1corr=layer.analyzer(tempst, extst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


te2corr=layer.analyzer(tempst, extst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)





compare.layered(te0, te1caus, te2caus, te1corr, te2corr,
   ML.IC="AICc") 
# 

summary(te2corr)
#




anova(te0, te1caus)$Pr[2]
# 

anova(te0, te1caus)$Pr[2]
# 

anova(te0, te1corr)$Pr[2]
# 

anova(te0, te2corr)$Pr[2]
# 

save.image(file=sprintf("globaltemp_extdone_run%1d.Rdata",run.nr))


#####################################
# Marine area vs sampling:
#####################################

ts0=layer.analyzer(tempst, sampst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ts1caus=layer.analyzer(tempst, sampst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ts2caus=layer.analyzer(tempst, sampst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ts1corr=layer.analyzer(tempst, sampst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


ts2corr=layer.analyzer(tempst, sampst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=8,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)





compare.layered(ts0, ts1caus, ts2caus, ts1corr, ts2corr,
   ML.IC="AICc")        
# 



summary(ts1caus)
#


anova(ts0, ts1caus)$Pr[2]
# 

anova(ts0, ts2caus)$Pr[2]
# 

anova(ts0, ts1corr)$Pr[2]
# 

anova(ts0, ts2corr)$Pr[2]
# 



save.image(file=sprintf("globaltemp_sampdone_run%1d.Rdata",run.nr))



