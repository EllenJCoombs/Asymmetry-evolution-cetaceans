
#Phylo ANOVAs 
#We ran phylogenetically corrected ANOVAs on each of the different scenarios using the R package nlme (v.3.1-137)
#and function ‘gls’ to test for correlations between the level of asymmetry seen in the skull and the potential scenarios (or regimes)
#‘gls’ allows for a more flexible model with better power. 
#Simulations were run using a ‘Pagel’s Lambda’ (λ) correlation structure (corPagel) in the ape package which is derived from a 
#Brownian Motion model by multiplying the covariances by λ. 

library(ape)
library(nlme)

#Do this for all of the different 'regimes' 
fit<-gls(sum.radii ~ as.factor(regime_split), correlation=corPagel(phy=tree), data = data) # can use corBrownian
anova(fit)

#Use corPagel for Pagal's lambda 

#Adjusting for false positive results with Benjamini-Hochberg 
#Post hoc tests - to check that this is not a false signal based on the large number
#of models we have run

p <- c(X.XX, X.XX, X.XX, X.XX, X.XX) # original p-values from ANOVAs 
p.adjust(p,method="BH") #ouputed Benjamin-Hochberg corrected results 
