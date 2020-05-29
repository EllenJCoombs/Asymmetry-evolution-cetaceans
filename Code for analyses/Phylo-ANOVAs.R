


library(ape)
library(nlme)

#Phylo ANOVAs 
#We ran phylogenetically corrected ANOVAs on each of the different scenarios using the R package nlme (v.3.1-137)
#and function ‘gls’ to test for correlations between the level of asymmetry seen in the skull and the potential scenarios (or regimes)
#‘gls’ allows for a more flexible model with better power. 
#Simulations were run using a ‘Pagel’s Lambda’ (λ) correlation structure (corPagel) in the ape package which is derived from a 
#Brownian Motion model by multiplying the covariances by λ. 


#Do this for all of the different 'regimes' 
fit<-gls(sum.radii ~ as.factor(age), correlation=corBrownian(phy=tree), data = data) # or better use corPagel
anova(fit)

#Use corPagel for Pagal's lambda 
