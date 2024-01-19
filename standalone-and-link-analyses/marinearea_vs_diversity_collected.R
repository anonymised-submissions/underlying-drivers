load("marinearea_vs_diversity_collected.Rdata")
# 24 runs

library(layeranalyzer)


#####################################
# Marine area vs origination:
#####################################

# 1= null model, 2=caus,m1->o1, 3=caus,m1->o2, 
# 4=corr,m1-o1, 5=corr, m1-o2

# Hypothesis tests:

anova(mo[[1]], mo[[2]])$Pr[2]
# 0.0128

anova(mo[[1]], mo[[3]])$Pr[2]
# 0.0140

anova(mo[[1]], mo[[4]])$Pr[2]
# 0.00104

anova(mo[[1]], mo[[5]])$Pr[2]
# 0.00189

# Bonferroni-limit: 0.05/42=0.00119
# That test of corr,m1-o1 actually snuck under that!



# parameter estimates:
for(m in 1:5)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(mo[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}

# corr_log10_surface_area,1_origts,1 : 0.959960





#####################################
# Marine area vs extinction:
#####################################


# 1= null model, 2=caus,t1->e1, 3=caus,t1->e2, 
# 4=corr,t1-e1, 5=corr, t1-e2

# Hypothesis tests:

anova(me[[1]], me[[2]])$Pr[2]
# 0.327

anova(me[[1]], me[[3]])$Pr[2]
# 0.339

anova(me[[1]], me[[4]])$Pr[2]
# 0.020

anova(me[[1]], me[[5]])$Pr[2]
# 0.0015


# parameter estimates:
for(m in 1:5)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(me[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}

# corr_log10_surface_area,1_extts,2 : 0.892


#####################################
# Marine area vs sampling:
#####################################

# 1= null model, 2=caus,t1->s1, 3=caus,t1->s2, 
# 4=corr,t1-s1, 5=corr, t1-s2

# Hypothesis tests:

anova(ms[[1]], ms[[2]])$Pr[2]
# 0.227

anova(ms[[1]], ms[[3]])$Pr[2]
# 0.357

anova(ms[[1]], ms[[4]])$Pr[2]
# 0.444

anova(ms[[1]], ms[[5]])$Pr[2]
# 0.180


# parameter estimates:
for(m in 1:5)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(ms[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}


