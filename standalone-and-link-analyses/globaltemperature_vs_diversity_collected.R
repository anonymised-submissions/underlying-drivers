load("globaltemperature_vs_diversity_collected.Rdata")
# 32 runs


#####################################
# Global temperature vs origination:
#####################################

# 1= null model, 2=caus,t1->o1, 3=caus,t1->o2, 
# 4=corr,t1-o1, 5=corr, t1-o2

# Hypothesis tests:

anova(to[[1]], to[[2]])$Pr[2]
# 0.0684

anova(to[[1]], to[[3]])$Pr[2]
# 0.0687

anova(to[[1]], to[[4]])$Pr[2]
# 0.0297

anova(to[[1]], to[[5]])$Pr[2]
# 0.0282


# parameter estimates:
for(m in 1:5)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(to[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}




#####################################
# Global temperature vs extinction:
#####################################


# 1= null model, 2=caus,t1->e1, 3=caus,t1->e2, 
# 4=corr,t1-e1, 5=corr, t1-e2

# Hypothesis tests:

anova(te[[1]], te[[2]])$Pr[2]
#  1

anova(te[[1]], te[[3]])$Pr[2]
# 1

anova(te[[1]], te[[4]])$Pr[2]
# 0.104

anova(te[[1]], te[[5]])$Pr[2]
# 0.208

# parameter estimates:
for(m in 1:5)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(te[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}


#####################################
# Global temperature vs sampling:
#####################################

# 1= null model, 2=caus,t1->s1, 3=caus,t1->s2, 
# 4=corr,t1-s1, 5=corr, t1-s2
# Hypothesis tests:

anova(ts[[1]], ts[[2]])$Pr[2]
# 0.258

anova(ts[[1]], ts[[3]])$Pr[2]
# 0.305

anova(ts[[1]], ts[[4]])$Pr[2]
# 1

anova(ts[[1]], ts[[5]])$Pr[2]
# 0.118

# parameter estimates:
for(m in 1:5)
{
  print.srcref(sprintf("Model %d:",m))
  print(summary.layered(ts[[m]]))

  if(m<7)
  {
    print.srcref("") 
    print.srcref("###############################")
    print.srcref("") 
  }
}



