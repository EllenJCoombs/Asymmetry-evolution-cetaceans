

CODE CURRENTLY UNDER REVIEW - PLEASE DO NOT USE UNTIL CITATION HAS BEEN MADE 

#Running evolutionary models on traits (asymmetry)
#The first model focuses on 'shifts'
#The second model focuses on 'jumps'

#These models were run on:
#1. the full landmark radii dataset
#2. the landmark radii dataset with the rostrum landamrks removed 
#3. the phylogeny that includes only species that appear in a character matrix 

library(ape)
library(coda) #mcmc
library(geiger)
library(phytools)
library(mvMORPH) #contmap 
library(evobiR) # to reorder data

setwd("D:/Checkpoint - ICVM/Asymmetry project")
#Read in subtree 
tree <- read.nexus('tree10.nexus')
#read in sum radii for specimens or your trait data 
data <- read.csv('xx.csv', row.names = 1) #turning row into species names 

#Make sure trait data and phylogeny are in the same order and double-check
Reordered <- ReorderData(tree_rearranged, data, taxa.names="row names")
data <- Reordered 

########################################
#                                      #
#   REARRANGE TREE - if required       #
#                                      #
########################################

#Need to convert tree from polytomy to binary - remove polytomies 
subtree <- multi2di(tree, random = TRUE)
subtree_unroot <- unroot(subtree)

subtree_unroot1 <- multi2di(subtree_unroot)
data <- data[subtree_unroot1$tip.label, ] #this takes the label from the dataset
data <- as.vector(data$sum.radii) #pulls out a vector (not matrix) using the $ sign - just the radii here 
names(data) = subtree_unroot1$tip.label

# Set parameters 
prop = 1.5        # tuning parameter for the proposal window (it's just a starting value it will be improved in the model fit I think)
iterations = 10e6  # number of iterations for the chain (5 million)
sampling = 10000   # sampling frequency of the parameters (i.e. when you store the parameter along the chain)
model = "rbm"     # relaxed brownian motion (the model used in Eastman et al. 2011. More or less similar to BayesTrait)
filename = paste("testmcmc-rjmcmcREARRANGED.log", sep="", collapse="")

# Half-Cauchy distribution (see also Gelman 2006) - I used it for the prior density of the rate scalar instead of the default exponential distribution
dhalfCauchy <- function(x, scale=1, log=F){
 if(any(x<0)) return(-Inf)
 density <- 2 /(pi * (1+(x/scale)^2))
 if(log==TRUE) density<-log(density/scale)
 return(density/scale)
}

# Define priors for the hyperparameter and measurment error (SE)
ratePrior <- function(x) dhalfCauchy(x, 25, log=TRUE)
sePrior <- function(x) dhalfCauchy(x, 25, log=TRUE)

# run the rjmcmc
rjmcmc.bm(tree_rearranged, data, prop.width=prop, ngen=iterations, samp=sampling, filebase=filename,
         simple.start=TRUE, type=model, dlnRATE=ratePrior, dlnSE=sePrior) # note: here I assume you're estimating a nuisance parameter as well

# retrieve the chain
chain_rearranged <- load.rjmcmc("relaxedBM.testmcmc-rjmcmcREARRANGED.log")

result=plot(chain_rearranged, par="shift", type= "fan", legend = F, cex = 0.6)

#Extra things from Ryan Felice to look at trace of the chain 
mcmc(chain110$log) #which column do you want to look at? 
#plot the trace of the chain 
plot(mcmc(chain$log[,8]))#for caterpillar plot 

#we can see then for diagnostic plots and retrieving parameters...etc
#if you want to use the default priors just remove dlnRATE=ratePrior, dlnSE=sePrior from the rjmcmc function

#shifts are cladewise 
plotShifts(reorder(subtree_unroot1,"cladewise") , chain=result$median.rates, color=c("blue","light blue", "green", "yellow", "orange", "dark orange", "red"))

# plot the trace of the chain
plot(mcmc(chain$log[,1:3]))

#Running the tree with the traits on the branches - by colour 
scaled_tree=reorder(subtree_unroot1,"cladewise")
scaled_tree$edge.length = scaled_tree$edge.length * result$median.rates
reconstructed_radii = fastAnc(scaled_tree, x = trait)

traits_plot <- contMap(tree_rearranged, trait, user=reconstructed_radii, cex=0.05)

#scaled_tree=reorder(subtree_unroot1,"cladewise")
#scaled_tree$edge.length = scaled_tree$edge.length * result$median.rates
#reconstructed_radii = fastAnc(scaled_tree, x = trait)
#contMap(subtree_unroot1, trait, user=reconstructed_radii, cex=0.6)

#invert colors so that red is fast 
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}
plot(setMap(traits_plot, invert=TRUE), outline=FALSE, lwd=4, fsize = 0.5, type = "fan") #label.offset doesn't work on an unrooted tree


##############################
#                            #
#   SECOND MODEL WITH JUMPS  #
#                            #
##############################

data <- read.csv('xx.csv', row.names = 1) #turning row into species names 

#Need to convert tree from polytomy to binary 
subtree <- multi2di(tree, random = TRUE)
subtree_unroot <- unroot(subtree)

subtree_unroot1 <- multi2di(subtree_unroot)
data <- data[subtree_unroot1$tip.label, ] #this takes the label from the dataset
data <- as.vector(data$sum.radii) #pulls out a vector (not matrix) using the $ sign - just the radii here 
names(data) = subtree_unroot1$tip.label

# Some parameters
prop = 1.5        # tuning parameter for the proposal window (it's just a starting value it will be improved in the model fit I think)
iterations = 10e6  # number of iterations for the chain
sampling = 10000 
model = "jump-rbm"     # relaxed brownian motion (the model used in Eastman et al. 2011. More or less similar to BayesTrait)
filename2 = paste("relaxedBM.testmcmc-jumprjmcmcREARRANGED.log", sep="", collapse="")

# run the rjmcmc
rjmcmc.bm(tree_rearranged, data, prop.width=prop, ngen=iterations, samp=sampling, filebase=filename2,
          simple.start=TRUE, type=model) # note: here I assume you're estimating a nuisance parameter as well

# retrieve the chain
chain_rearranged2 <- load.rjmcmc("jump-relaxedBM.relaxedBM.testmcmc-jumprjmcmcREARRANGED.log")

result=plot(chain210, par="shift", legend = F, cex = 0.6)
result=plot(chain_rearranged2, par="jumps", type = "fan", legend = F, cex = 0.6) #jumps 

#we can see then for diagnostic plots and retrieving parameters...etc
#if you want to use the default priors just remove dlnRATE=ratePrior, dlnSE=sePrior from the rjmcmc function

###############################################
#                                             #
#  Be sure to have run 'Clavel_plotShifts.R'  #
#                                             #
###############################################

#Looking at jump probability
plotShifts(reorder(subtree_unroot1,"cladewise") , chain=result$jump.probability, color=c("blue","green","orange","red"))

# plot the trace of the chain
plot(mcmc(chain210$log[,1:3]))
plot(mcmc(chain210$jumps[,1:3]))

plot(tree);
axisPhylo()
mcmc(chain210$log) #which column do you want to look at? We want the root 
#plot the trace of the chain 
plot(mcmc(chain210$log[,8]))
