run.nr.start=1
num.run=32

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

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
to=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("globaltemp_sampdone_run%d.Rdata",r))

  if(is.finite(to0$ML))
  if(to0$ML<1e+10)
  if(to0$ML>ml.best[1])
  {
     to[[1]]=to0
     ml.best[1]=to0$ML
  }

  if(is.finite(to1caus$ML))
  if(to1caus$ML<1e+10)
  if(to1caus$ML>ml.best[2])
  {
     to[[2]]=to1caus
     ml.best[2]=to1caus$ML
  }

  if(is.finite(to2caus$ML))
  if(to2caus$ML<1e+10)
  if(to2caus$ML>ml.best[3])
  {
     to[[3]]=to2caus
     ml.best[3]=to2caus$ML
  }
  
  if(is.finite(to1corr$ML))
  if(to1corr$ML<1e+10)
  if(to1corr$ML>ml.best[4])
  {
     to[[4]]=to1corr
     ml.best[4]=to1corr$ML
  }

  if(is.finite(to2corr$ML))
  if(to2corr$ML<1e+10)
  if(to2corr$ML>ml.best[5])
  {
     to[[5]]=to2corr
     ml.best[5]=to2corr$ML
  }
}


compare.layered(to, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -411.8567        7.21168
#Model   2        -410.2384       36.38094
#Model   3        -410.6568       23.94328
#Model   4        -411.1113       15.19715
#Model   5        -410.9837       17.26694



# Hypothesis tests:

anova(to[[1]], to[[2]])$Pr[2]
# 0.0684

anova(to[[1]], to[[3]])$Pr[2]
# 0.1040

anova(to[[1]], to[[4]])$Pr[2]
# 0.0572

anova(to[[1]], to[[5]])$Pr[2]
# 0.04906



#####################################
# Global temperature vs extinction:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
te=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("globaltemp_sampdone_run%d.Rdata",r))

  if(is.finite(te0$ML))
  if(te0$ML<1e+10)
  if(te0$ML>ml.best[1])
  {
     te[[1]]=te0
     ml.best[1]=te0$ML
  }

  if(is.finite(te1caus$ML))
  if(te1caus$ML<1e+10)
  if(te1caus$ML>ml.best[2])
  {
     te[[2]]=te1caus
     ml.best[2]=te1caus$ML
  }

  if(is.finite(te2caus$ML))
  if(te2caus$ML<1e+10)
  if(te2caus$ML>ml.best[3])
  {
     te[[3]]=te2caus
     ml.best[3]=te2caus$ML
  }
  
  if(is.finite(te1corr$ML))
  if(te1corr$ML<1e+10)
  if(te1corr$ML>ml.best[4])
  {
     te[[4]]=te1corr
     ml.best[4]=te1corr$ML
  }

  if(is.finite(te2corr$ML))
  if(te2corr$ML<1e+10)
  if(te2corr$ML>ml.best[5])
  {
     te[[5]]=te2corr
     ml.best[5]=te2corr$ML
  }
}


compare.layered(te, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -431.5616       34.49750
#Model   2        -432.8338        9.66719
#Model   3        -432.9832        8.32500
#Model   4        -432.0484       21.20302
#Model   5        -431.8327       26.30729

summary(te1corr)
#

# Hypothesis tests:

anova(te[[1]], te[[2]])$Pr[2]
#  1

anova(te[[1]], te[[3]])$Pr[2]
# 1

anova(te[[1]], te[[4]])$Pr[2]
# 0.2828

anova(te[[1]], te[[5]])$Pr[2]
# 0.2080

#####################################
# Global temeprature vs sampling:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
ts=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("globaltemp_sampdone_run%d.Rdata",r))

  if(is.finite(ts0$ML))
  if(ts0$ML<1e+10)
  if(ts0$ML>ml.best[1])
  {
     ts[[1]]=ts0
     ml.best[1]=ts0$ML
  }

  if(is.finite(ts1caus$ML))
  if(ts1caus$ML<1e+10)
  if(ts1caus$ML>ml.best[2])
  {
     ts[[2]]=ts1caus
     ml.best[2]=ts1caus$ML
  }

  if(is.finite(ts2caus$ML))
  if(ts2caus$ML<1e+10)
  if(ts2caus$ML>ml.best[3])
  {
     ts[[3]]=ts2caus
     ml.best[3]=ts2caus$ML
  }
  
  if(is.finite(ts1corr$ML))
  if(ts1corr$ML<1e+10)
  if(ts1corr$ML>ml.best[4])
  {
     ts[[4]]=ts1corr
     ml.best[4]=ts1corr$ML
  }

  if(is.finite(ts2corr$ML))
  if(ts2corr$ML<1e+10)
  if(ts2corr$ML>ml.best[5])
  {
     ts[[5]]=ts2corr
     ml.best[5]=ts2corr$ML
  }
}


compare.layered(ts, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1        -384.5678       22.21318
#Model   2        -384.2840       29.50513
#Model   3        -384.3782       26.85088
#Model   4        -385.7292        6.95436
#Model   5        -384.9960       14.47645

summary(ts1corr)
#

# Hypothesis tests:

anova(ts[[1]], ts[[2]])$Pr[2]
# 0.2599

anova(ts[[1]], ts[[3]])$Pr[2]
# 0.2856

anova(ts[[1]], ts[[4]])$Pr[2]
# 1

anova(ts[[1]], ts[[5]])$Pr[2]
# 0.2596


save.image("globaltemperature_vs_diversity_collected.Rdata")
