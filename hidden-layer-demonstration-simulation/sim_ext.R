source("numeric_stoch_diff.R")

t1=0.0174
sigma1=1
sigma2=0.0299
init1=-4.336193
init2=-3.000669

lv.der=function(x)
  {
    ret=rep(0,2)
    ret[1]=-(x[1]-x[2])/t1
    ret[2]=0
    return(ret)
  }

lv.diff=function(x)
  {
    g.ret=matrix(rep(0,2*2),nrow=2)
    g.ret[1,1]=sigma1
    g.ret[2,2]=sigma2
    return(g.ret)
  }


x0.lv=c(init1,init2)

traj.lv=explicit.stoch(x0.lv, lv.der, lv.diff,530,0.005)
traj.lv$t=traj.lv$t-530

plot(traj.lv$t, traj.lv$x[,1], type="l",lwd=1,col="grey")
lines(traj.lv$t, traj.lv$x[,2], type="l",lwd=3,col="black")

plot(traj.lv$t, traj.lv$x[,2], type="l",lwd=2,col="black")

meas=sort(sample(1:length(traj.lv$t),90))
t.meas=traj.lv$t[meas]
v.meas=traj.lv$x[meas,1]+0.1*rnorm(90)

plot(t.meas,v.meas)
lines(traj.lv$t, traj.lv$x[,2], type="l",lwd=3,col="black")


library(layeranalyzer)
simext=layer.data.series(time.points=t.meas, 
                      value.points=v.meas,
                      std.dev=rep(0.1,90),
                      name="simext")

pr=layer.prior(mu=c(-5,5), dt=c(0.1,10), sigma=c(0.0001,1),
    init = c(-5,5) , lin = c(-1,1), beta = c(-1,1), 
    obs = c(0.01,0.5))

simstruct=layer.series.structure(simext,
  numlayers=2,no.pull=TRUE,init.time=-1000,prior=pr)

res=layer.analyzer(simstruct, 
  num.MCMC = 1000, spacing = 10, burnin = 10000,
  do.model.likelihood = FALSE,
  #silent.mode=FALSE,  
  #do.maximum.likelihood = TRUE, 
  maximum.likelihood.numstart = 600, 
  smoothing.specs = list(do.smoothing = TRUE, 
        smoothing.time.diff = 1, smoothing.start = -1000, 
        smoothing.end = 0, num.smooth.per.mcmc = 10,
        do.return.smoothing.samples = FALSE) )


png("reconstruct_sim_2layer.png",width=1800,height=1400)
par(cex=5)
plot(t.meas,v.meas, xlab="Time",ylab="Value")
lines(traj.lv$t, traj.lv$x[,1], type="l",lwd=1,col="grey")
lines(traj.lv$t, traj.lv$x[,2], type="l",lwd=4,col="black")
lines(res$process.time.points, res$process.mean[2,],
  col="blue",lwd=4)
lines(res$process.time.points, res$process.lower95[2,],
  col="green",lwd=2)
lines(res$process.time.points, res$process.upper95[2,],
  col="green",lwd=2)
points(t.meas,v.meas,lwd=3)
dev.off()


save(traj.lv, file="simulated_time_series.Rdata")
save(t.meas, file="sim_measurement_time_points.Rdata")
save(v.meas, file="sim_measurements.Rdata")
save(res, file="layer_analyzer_result.Rdata")






