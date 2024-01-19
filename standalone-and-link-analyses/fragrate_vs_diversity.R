
run.nr=4

library(layeranalyzer)


# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Read fragefentation rate
fragts=read.layer.data.series("smoothed_frag_rate.csv",sep=";",
  header=T, column.type=c("time","value"), name="fragrate")

#subsample it. This is too much!
i=seq(1,length(fragts$time),10)
fragts2=layer.data.series(time.points=fragts$time[i], 
        value.points=fragts$value[i], name="fragrate")


# fragmentation rate structure:
fragst<-layer.series.structure(fragts2,numlayers=1,
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
# Fragmentation rate vs origination:
#####################################

t1=Sys.time()
fo0=layer.analyzer(fragst, origst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)
t2=Sys.time(); t2-t1
ml2=fo0$ML


fo1caus=layer.analyzer(fragst, origst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fo2caus=layer.analyzer(fragst, origst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fo1corr=layer.analyzer(fragst, origst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fo2corr=layer.analyzer(fragst, origst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)




# Hypothesis tests:

anova(fo0, fo1caus)$Pr[2]
# 0.01447855

anova(fo0, fo1caus)$Pr[2]
# 0.01447855

anova(fo0, fo1corr)$Pr[2]
# 0.001445325

anova(fo0, fo2corr)$Pr[2]
# 0.01179452


save.image(file=sprintf("fragorig%d.Rdata",run.nr))




#####################################
# Marine area vs extinction:
#####################################

fe0=layer.analyzer(fragst, extst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fe1caus=layer.analyzer(fragst, extst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fe2caus=layer.analyzer(fragst, extst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fe1corr=layer.analyzer(fragst, extst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fe2corr=layer.analyzer(fragst, extst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)



anova(fe0, fe1caus)$Pr[2]
# 0.3759899

anova(fe0, fe1caus)$Pr[2]
# 0.3759899

anova(fe0, fe1corr)$Pr[2]
# 0.02417614

anova(fe0, fe2corr)$Pr[2]
# 0.003623773


save.image(file=sprintf("fragext%d.Rdata",run.nr))



#####################################
# Fragmentation rate vs sampling:
#####################################

fs0=layer.analyzer(fragst, sampst,
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fs1caus=layer.analyzer(fragst, sampst,
       causal=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fs2caus=layer.analyzer(fragst, sampst,
       causal=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fs1corr=layer.analyzer(fragst, sampst,
       corr=cbind(c(1,1,2,1)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


fs2corr=layer.analyzer(fragst, sampst,
       corr=cbind(c(1,1,2,2)),
       num.MCMC=1000,spacing=10,burnin=10000,num.temp=4,
       do.maximum.likelihood=TRUE,
       maximum.likelihood.numstart=600,
       silent.mode=FALSE,talkative.burnin=TRUE,
       talkative.likelihood=FALSE)


anova(fs0, fs1caus)$Pr[2]
# 0.1617471

anova(fs0, fs2caus)$Pr[2]
# 0.2288407

anova(fs0, fs1corr)$Pr[2]
# 0.187151

anova(fs0, fs2corr)$Pr[2]
# 0.3209712


save.image(file=sprintf("fragsamp%d.Rdata",run.nr))




