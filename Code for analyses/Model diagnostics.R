
#Running model diagnostics from 'Clavel_model_shifts_jumps.R' and 'Discrete_trait_models.R'
library(geiger)
library(coda)
?coda.options

setwd("D:/Checkpoint - ICVM/Asymmetry project")
chain <- load.rjmcmc("relaxedBM.testmcmc-rjmcmcREMOVED.log")
chain2 <- load.rjmcmc("jump-relaxedBM.relaxedBM.testmcmc-jumprjmcmcREMOVED.log")
#Extra things from Ryan Felice to look at trace of the chain 
mcmc(chain$log) #which column do you want to look at?

str(chain$log)

#plot the trace of the chain - we want 'SE' 
plot(mcmc(chain$log[,8]))#for caterpillar plot
#plot the trace of the chain 
plot(mcmc(chain$log[,]))#for caterpillar plot

#Looking at ESS 
effectiveSize(mcmc(chain$log)) 
effectiveSize(mcmc(chain2$log))

#Gelman diagnostics 
#Look for the convergence of the chains 
mychains=mcmc.list(mcmc(chain$log[,8]), mcmc(chain2$log[,8])) 
gelman.diag(mychains) #what does this number mean? See text below. We get 1.02 and 1.09

#Ellen - The gelman.diag gives you the scale reduction factors for each parameter. 
#Approximate convergence is diagnosed when the upper limit is close to 1. 
#A factor of 1 means that between variance and within chain variance are equal. 
#Larger values mean that there is still a notable difference between chains. 
#A values of 1.1 or below is generally accepted. 

#Look at the convergence of the chains - one chain is red, one is black
gelman.plot((mychains))
#Iterations of the two chains (should overlay one another for good convergence)
plot(mychains)
© 2020 GitHub, Inc.
