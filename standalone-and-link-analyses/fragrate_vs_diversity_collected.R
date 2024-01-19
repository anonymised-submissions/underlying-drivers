library(layeranalyzer)
load("fragrate_vs_diversity_collected.Rdata")



#####################################
# Fragmentation rate vs origination:
#####################################

# 1= null model, 2=1caus, 3=2caus, 4=1corr, 5=2corr 

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

# Hypothesis tesfs:

anova(fs[[1]], fs[[2]])$Pr[2]
# 0.727

anova(fs[[1]], fs[[3]])$Pr[2]
# 0.785

anova(fs[[1]], fs[[4]])$Pr[2]
# 0.613

anova(fs[[1]], fs[[5]])$Pr[2]
# 1
