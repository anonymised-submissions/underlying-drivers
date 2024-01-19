

load("sealevel_vs_diversity_collected.Rdata")
# 25 runs

library(layeranalyzer)



#####################################
# Sealevel vs origination:
#####################################

# 1=null model, 2=caus,s1->o1, 3=caus,s2->o1, 
# 4=caus, s1->o2, 5=caus, s2->o2, 6=corr,s2-o1, 
# 7=corr, s2-o2



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

# parameter estimates:
for(m in 1:7)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(so[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}




#####################################
# Sealevel vs extinction:
#####################################

# 1=null model, 2=caus,s1->e1, 3=caus,s2->e1, 
# 4=caus, s1->e2, 5=caus, s2->e2, 6=corr,s2-e1, 
# 7=corr, s2-e2


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


# parameter estimates:
for(m in 1:7)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(se[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}



#####################################
# Sealevel (s) vs sampling (S):
#####################################

# 1=null model, 2=caus,s1->S1, 3=caus,s2->S1, 
# 4=caus, s1->S2, 5=caus, s2->S2, 6=corr,s2-S1, 
# 7=corr, s2-S2



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



# parameter estimates:
for(m in 1:7)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(ss[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}



