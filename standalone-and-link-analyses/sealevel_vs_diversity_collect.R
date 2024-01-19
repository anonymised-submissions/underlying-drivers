run.nr.start=1
num.run=32

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

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
so=list()
ml.best=rep(-1e+200,7)
ml0=rep(0,num.run)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("sealevel_vs_samp%d.Rdata",r))

  ml0[r]=so0$ML

  if(is.finite(so0$ML))
  if(so0$ML<1e+10)
  if(so0$ML>ml.best[1])
  {
     so[[1]]=so0
     ml.best[1]=so0$ML
  }

  if(is.finite(s1o1caus$ML))
  if(s1o1caus$ML<1e+10)
  if(s1o1caus$ML>ml.best[2])
  {
     so[[2]]=s1o1caus
     ml.best[2]=s1o1caus$ML
  }

  if(is.finite(s2o1caus$ML))
  if(s2o1caus$ML<1e+10)
  if(s2o1caus$ML>ml.best[3])
  {
     so[[3]]=s2o1caus
     ml.best[3]=s2o1caus$ML
  }

  if(is.finite(s1o2caus$ML))
  if(s1o2caus$ML<1e+10)
  if(s1o2caus$ML>ml.best[4])
  {
     so[[4]]=s1o2caus
     ml.best[4]=s1o2caus$ML
  }

  if(is.finite(s2o2caus$ML))
  if(s2o2caus$ML<1e+10)
  if(s2o2caus$ML>ml.best[5])
  {
     so[[5]]=s2o2caus
     ml.best[5]=s2o2caus$ML
  }

  if(is.finite(s2o1corr$ML))
  if(s2o1corr$ML<1e+10)
  if(s2o1corr$ML>ml.best[6])
  {
     so[[6]]=s2o1corr
     ml.best[6]=s2o1corr$ML
  }

  if(is.finite(s2o2corr$ML))
  if(s2o2corr$ML<1e+10)
  if(s2o2corr$ML>ml.best[7])
  {
     so[[7]]=s2o2corr
     ml.best[7]=s2o2corr$ML
  }
}


compare.layered(so, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         545.7748       15.68982
#Model   2         544.7151        5.43745
#Model   3         546.1832       23.60439
#Model   4         544.5493        4.60648
#Model   5         545.8848       17.51317
#Model   6         545.2857        9.62000
#Model   7         546.1800       23.52867

summary(s2o1caus)
#

# Hypothesis tests:

anova(so[[1]], so[[2]])$Pr[2]
# 0.158

anova(so[[1]], so[[3]])$Pr[2]
# 0.195

anova(so[[1]], so[[4]])$Pr[2]
# 0.399

anova(so[[1]], so[[5]])$Pr[2]
# 0.170

anova(so[[1]], so[[6]])$Pr[2]
# 0.547

anova(so[[1]], so[[7]])$Pr[2]
# 0.119




#####################################
# Sealevel vs extinction:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
se=list()
ml.best=rep(-1e+200,7)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("sealevel_vs_samp%d.Rdata",r))

  if(is.finite(se0$ML))
  if(se0$ML<1e+10)
  if(se0$ML>ml.best[1])
  {
     se[[1]]=se0
     ml.best[1]=se0$ML
  }

  if(is.finite(s1e1caus$ML))
  if(s1e1caus$ML<1e+10)
  if(s1e1caus$ML>ml.best[2])
  {
     se[[2]]=s1e1caus
     ml.best[2]=s1e1caus$ML
  }

  if(is.finite(s1e2caus$ML))
  if(s1e2caus$ML<1e+10)
  if(s1e2caus$ML>ml.best[3])
  {
     se[[3]]=s1e2caus
     ml.best[3]=s1e2caus$ML
  }

  if(is.finite(s2e1caus$ML))
  if(s2e1caus$ML<1e+10)
  if(s2e1caus$ML>ml.best[4])
  {
     se[[4]]=s2e1caus
     ml.best[4]=s2e1caus$ML
  }

  if(is.finite(s2e2caus$ML))
  if(s2e2caus$ML<1e+10)
  if(s2e2caus$ML>ml.best[5])
  {
     se[[5]]=s2e2caus
     ml.best[5]=s2e2caus$ML
  }

  if(is.finite(s2e1corr$ML))
  if(s2e1corr$ML<1e+10)
  if(s2e1corr$ML>ml.best[6])
  {
     se[[6]]=s2e1corr
     ml.best[6]=s2e1corr$ML
  }

  if(is.finite(s2e2corr$ML))
  if(s2e2corr$ML<1e+10)
  if(s2e2corr$ML>ml.best[7])
  {
     se[[7]]=s2e2corr
     ml.best[7]=s2e2corr$ML
  }
}


compare.layered(se, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         525.7631       19.82797
#Model   2         525.2474       11.83928
#Model   3         524.8528        7.97932
#Model   4         525.6470       17.65513
#Model   5         525.9715       24.42214
#Model   6         525.0943       10.15808
#Model   7         524.8701        8.11807

summary(s2e2caus)
#

# Hypothesis tests:

anova(se[[1]], se[[2]])$Pr[2]
# 1

anova(se[[1]], se[[3]])$Pr[2]
# 1

anova(se[[1]], se[[4]])$Pr[2]
# 1

anova(se[[1]], se[[5]])$Pr[2]
# 1

anova(se[[1]], se[[6]])$Pr[2]
# 0.376

anova(se[[1]], se[[7]])$Pr[2]
# 0.489


#####################################
# Sealevel vs sampling:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
ss=list()
ml.best=rep(-1e+200,7)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("sealevel_vs_samp%d.Rdata",r))

  if(is.finite(se0$ML))
  if(se0$ML<1e+10)
  if(ss0$ML>ml.best[1])
  {
     ss[[1]]=ss0
     ml.best[1]=ss0$ML
  }

  if(is.finite(s1s1caus$ML))
  if(s1s1caus$ML<1e+10)
  if(s1s1caus$ML>ml.best[2])
  {
     ss[[2]]=s1s1caus
     ml.best[2]=s1s1caus$ML
  }

  if(is.finite(s1s2caus$ML))
  if(s1s2caus$ML<1e+10)
  if(s1s2caus$ML>ml.best[3])
  {
     ss[[3]]=s1s2caus
     ml.best[3]=s1s2caus$ML
  }

  if(is.finite(s2s1caus$ML))
  if(s2s1caus$ML<1e+10)
  if(s2s1caus$ML>ml.best[4])
  {
     ss[[4]]=s2s1caus
     ml.best[4]=s2s1caus$ML
  }

  if(is.finite(s2s2caus$ML))
  if(s2s2caus$ML<1e+10)
  if(s2s2caus$ML>ml.best[5])
  {
     ss[[5]]=s2s2caus
     ml.best[5]=s2s2caus$ML
  }

  if(is.finite(s2s1corr$ML))
  if(s2s1corr$ML<1e+10)
  if(s2s1corr$ML>ml.best[6])
  {
     ss[[6]]=s2s1corr
     ml.best[6]=s2s1corr$ML
  }

  if(is.finite(s2s2corr$ML))
  if(s2s2corr$ML<1e+10)
  if(s2s2corr$ML>ml.best[7])
  {
     ss[[7]]=s2s2corr
     ml.best[7]=s2s2corr$ML
  }
}


compare.layered(ss, ML.IC="AICc")
#          weight=-0.5*AICc Post. Prob.(%)
#Model   1         573.2722       21.92316
#Model   2         572.1049        6.82283
#Model   3         571.8384        5.22661
#Model   4         572.3904        9.07768
#Model   5         572.1857        7.39681
#Model   6         573.9651       43.83508
#Model   7         571.9282        5.71784

summary(s2s1corr)
#

# Hypothesis tests:

anova(ss[[1]], ss[[2]])$Pr[2]
# 0.962

anova(ss[[1]], ss[[3]])$Pr[2]
# 0.251

anova(ss[[1]], ss[[4]])$Pr[2]
# 1

anova(ss[[1]], ss[[5]])$Pr[2]
# 0.187

anova(ss[[1]], ss[[6]])$Pr[2]
# 1

anova(ss[[1]], ss[[7]])$Pr[2]
# 0.421

save.image("sealevel_vs_diversity_collected.Rdata")

