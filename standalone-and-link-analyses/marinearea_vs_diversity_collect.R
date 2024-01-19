run.nr.start=1
num.run=24

library(layeranalyzer)


# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Read marine area
load("marinearea.Rdata")


# Global temperature structure:
marst<-layer.series.structure(marinearea,numlayers=1,
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

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
mo=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("marinearea_sampdone_run%d",r))

  if(is.finite(mo0$ML))
  if(mo0$ML<1e+10)
  if(mo0$ML>ml.best[1])
  {
     mo[[1]]=mo0
     ml.best[1]=mo0$ML
  }

  if(is.finite(mo1caus$ML))
  if(mo1caus$ML<1e+10)
  if(mo1caus$ML>ml.best[2])
  {
     mo[[2]]=mo1caus
     ml.best[2]=mo1caus$ML
  }

  if(is.finite(mo2caus$ML))
  if(mo2caus$ML<1e+10)
  if(mo2caus$ML>ml.best[3])
  {
     mo[[3]]=mo2caus
     ml.best[3]=mo2caus$ML
  }
  
  if(is.finite(mo1corr$ML))
  if(mo1corr$ML<1e+10)
  if(mo1corr$ML>ml.best[4])
  {
     mo[[4]]=mo1corr
     ml.best[4]=mo1corr$ML
  }

  if(is.finite(mo2corr$ML))
  if(mo2corr$ML<1e+10)
  if(mo2corr$ML>ml.best[5])
  {
     mo[[5]]=mo2corr
     ml.best[5]=mo2corr$ML
  }
}


compare.layered(mo, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -411.8567        7.21168
#Model   2        -410.2384       36.38094
#Model   3        -410.6568       23.94328
#Model   4        -411.1113       15.19715
#Model   5        -410.9837       17.26694



# Hypothesis tests:

anova(mo[[1]], mo[[2]])$Pr[2]
# 0.0128

anova(mo[[1]], mo[[3]])$Pr[2]
# 0.0140

anova(mo[[1]], mo[[4]])$Pr[2]
# 0.00104

anova(mo[[1]], mo[[5]])$Pr[2]
# 0.00189



#####################################
# Marine area vs extinction:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
me=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("marinearea_sampdone_run%d",r))

  if(is.finite(me0$ML))
  if(me0$ML<1e+10)
  if(me0$ML>ml.best[1])
  {
     me[[1]]=me0
     ml.best[1]=me0$ML
  }

  if(is.finite(me1caus$ML))
  if(me1caus$ML<1e+10)
  if(me1caus$ML>ml.best[2])
  {
     me[[2]]=me1caus
     ml.best[2]=me1caus$ML
  }

  if(is.finite(me2caus$ML))
  if(me2caus$ML<1e+10)
  if(me2caus$ML>ml.best[3])
  {
     me[[3]]=me2caus
     ml.best[3]=me2caus$ML
  }
  
  if(is.finite(me1corr$ML))
  if(me1corr$ML<1e+10)
  if(me1corr$ML>ml.best[4])
  {
     me[[4]]=me1corr
     ml.best[4]=me1corr$ML
  }

  if(is.finite(me2corr$ML))
  if(me2corr$ML<1e+10)
  if(me2corr$ML>ml.best[5])
  {
     me[[5]]=me2corr
     ml.best[5]=me2corr$ML
  }
}


compare.layered(me, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -431.5616       34.49750
#Model   2        -432.8338        9.66719
#Model   3        -432.9832        8.32500
#Model   4        -432.0484       21.20302
#Model   5        -431.8327       26.30729

summary(me1corr)
#

# Hypothesis mests:

anova(me[[1]], me[[2]])$Pr[2]
#  0.327

anova(me[[1]], me[[3]])$Pr[2]
# 0.339

anova(me[[1]], me[[4]])$Pr[2]
# 0.020

anova(me[[1]], me[[5]])$Pr[2]
# 0.0015

#####################################
# Marine area vs sampling:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
ms=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("marinearea_sampdone_run%d",r))

  if(is.finite(ms0$ML))
  if(ms0$ML<1e+10)
  if(ms0$ML>ml.best[1])
  {
     ms[[1]]=ms0
     ml.best[1]=ms0$ML
  }

  if(is.finite(ms1caus$ML))
  if(ms1caus$ML<1e+10)
  if(ms1caus$ML>ml.best[2])
  {
     ms[[2]]=ms1caus
     ml.best[2]=ms1caus$ML
  }

  if(is.finite(ms2caus$ML))
  if(ms2caus$ML<1e+10)
  if(ms2caus$ML>ml.best[3])
  {
     ms[[3]]=ms2caus
     ml.best[3]=ms2caus$ML
  }
  
  if(is.finite(ms1corr$ML))
  if(ms1corr$ML<1e+10)
  if(ms1corr$ML>ml.best[4])
  {
     ms[[4]]=ms1corr
     ml.best[4]=ms1corr$ML
  }

  if(is.finite(ms2corr$ML))
  if(ms2corr$ML<1e+10)
  if(ms2corr$ML>ml.best[5])
  {
     ms[[5]]=ms2corr
     ml.best[5]=ms2corr$ML
  }
}


compare.layered(ms, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -384.5678       22.21318
#Model   2        -384.2840       29.50513
#Model   3        -384.3782       26.85088
#Model   4        -385.7292        6.95436
#Model   5        -384.9960       14.47645

summary(ms1corr)
#

# Hypothesis tests:

anova(ms[[1]], ms[[2]])$Pr[2]
# 0.227

anova(ms[[1]], ms[[3]])$Pr[2]
# 0.357

anova(ms[[1]], ms[[4]])$Pr[2]
# 0.444

anova(ms[[1]], ms[[5]])$Pr[2]
# 0.180

save.image("marinearea_vs_diversity_collected.Rdata")
