run.nr.start=1
num.run=4

library(layeranalyzer)


# Prior:
pr=layer.prior(mu=c(-10,10),dt=c(0.1,500), 
  sigma=c(1e-9,1000), init=c(-10,10),
  lin=c(-0.01,0.01), beta=c(-1,1), obs=c(0.0001,2))
  

# Read fragmentation rate
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

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr 
fo=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("fragsamp%d.Rdata",r))

  if(is.finite(fo0$ML))
  if(fo0$ML<1e+10)
  if(fo0$ML>ml.best[1])
  {
     fo[[1]]=fo0
     ml.best[1]=fo0$ML
  }

  if(is.finite(fo1caus$ML))
  if(fo1caus$ML<1e+10)
  if(fo1caus$ML>ml.best[2])
  {
     fo[[2]]=fo1caus
     ml.best[2]=fo1caus$ML
  }

  if(is.finite(fo2caus$ML))
  if(fo2caus$ML<1e+10)
  if(fo2caus$ML>ml.best[3])
  {
     fo[[3]]=fo2caus
     ml.best[3]=fo2caus$ML
  }
  
  if(is.finite(fo1corr$ML))
  if(fo1corr$ML<1e+10)
  if(fo1corr$ML>ml.best[4])
  {
     fo[[4]]=fo1corr
     ml.best[4]=fo1corr$ML
  }

  if(is.finite(fo2corr$ML))
  if(fo2corr$ML<1e+10)
  if(fo2corr$ML>ml.best[5])
  {
     fo[[5]]=fo2corr
     ml.best[5]=fo2corr$ML
  }
}


compare.layered(fo, ML.IC="AICc")



# Hypothesis tests:

anova(fo[[1]], fo[[2]])$Pr[2]
# 1

anova(fo[[1]], fo[[3]])$Pr[2]
# 1

anova(fo[[1]], fo[[4]])$Pr[2]
# 1

anova(fo[[1]], fo[[5]])$Pr[2]
# 1



#####################################
# Fragmentation rate vs extinction:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
fe=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("fragsamp%d.Rdata",r))

  if(is.finite(fe0$ML))
  if(fe0$ML<1e+10)
  if(fe0$ML>ml.best[1])
  {
     fe[[1]]=fe0
     ml.best[1]=fe0$ML
  }

  if(is.finite(fe1caus$ML))
  if(fe1caus$ML<1e+10)
  if(fe1caus$ML>ml.best[2])
  {
     fe[[2]]=fe1caus
     ml.best[2]=fe1caus$ML
  }

  if(is.finite(fe2caus$ML))
  if(fe2caus$ML<1e+10)
  if(fe2caus$ML>ml.best[3])
  {
     fe[[3]]=fe2caus
     ml.best[3]=fe2caus$ML
  }
  
  if(is.finite(fe1corr$ML))
  if(fe1corr$ML<1e+10)
  if(fe1corr$ML>ml.best[4])
  {
     fe[[4]]=fe1corr
     ml.best[4]=fe1corr$ML
  }

  if(is.finite(fe2corr$ML))
  if(fe2corr$ML<1e+10)
  if(fe2corr$ML>ml.best[5])
  {
     fe[[5]]=fe2corr
     ml.best[5]=fe2corr$ML
  }
}


compare.layered(fe, ML.IC="AICc")

summary(fe1corr)
#

# Hypothesis tests:

anova(fe[[1]], fe[[2]])$Pr[2]
# 0.856

anova(fe[[1]], fe[[3]])$Pr[2]
# 0.831

anova(fe[[1]], fe[[4]])$Pr[2]
# 0.353

anova(fe[[1]], fe[[5]])$Pr[2]
# 0.614

#####################################
# Fragmentation rate vs sampling:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr
fs=list()
ml.best=rep(-1e+200,5)

for(r in run.nr.start-1+(1:num.run))
{
  load(sprintf("fragsamp%d.Rdata",r))

  if(is.finite(fs0$ML))
  if(fs0$ML<1e+10)
  if(fs0$ML>ml.best[1])
  {
     fs[[1]]=fs0
     ml.best[1]=fs0$ML
  }

  if(is.finite(fs1caus$ML))
  if(fs1caus$ML<1e+10)
  if(fs1caus$ML>ml.best[2])
  {
     fs[[2]]=fs1caus
     ml.best[2]=fs1caus$ML
  }

  if(is.finite(fs2caus$ML))
  if(fs2caus$ML<1e+10)
  if(fs2caus$ML>ml.best[3])
  {
     fs[[3]]=fs2caus
     ml.best[3]=fs2caus$ML
  }
  
  if(is.finite(fs1corr$ML))
  if(fs1corr$ML<1e+10)
  if(fs1corr$ML>ml.best[4])
  {
     fs[[4]]=fs1corr
     ml.best[4]=fs1corr$ML
  }

  if(is.finite(fs2corr$ML))
  if(fs2corr$ML<1e+10)
  if(fs2corr$ML>ml.best[5])
  {
     fs[[5]]=fs2corr
     ml.best[5]=fs2corr$ML
  }
}


compare.layered(fs, ML.IC="AICc")



# Hypothesis tesfs:

anova(fs[[1]], fs[[2]])$Pr[2]
# 0.727

anova(fs[[1]], fs[[3]])$Pr[2]
# 0.785

anova(fs[[1]], fs[[4]])$Pr[2]
# 0.613

anova(fs[[1]], fs[[5]])$Pr[2]
# 1

save.image("fragrate_vs_diversity_collected.Rdata")
