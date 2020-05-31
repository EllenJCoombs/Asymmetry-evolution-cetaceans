

CODE CURRENTLY UNDER REVIEW - PLEASE DO NOT USE UNTIL CITATION HAS BEEN MADE 


#Code modified from Julien Clavel 
#mvMORPH 
#See: Clavel et al 2015 - mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data

#Models testing whether changes in cetacean cranial asymmetry are associated with other discrete traits

###########################################                                   
#      Ancestral state reconstruction     #
#          + PCMs                         #
###########################################

#If not already loaded 
library(ape)
library(coda) #mcmc
library(geiger)
library(phytools)
library(mvMORPH)


# Function to transform ancestral states reconstruction to a SIMMAP like tree (to be added to mvMORPH)
paintAllTree <- function(tree, ancestral, tips){

 names_grps <- colnames(ancestral$lik.anc)
 statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
 combined = as.character(c(tips, statesNodes))
 treebis=tree
 for(i in sort(tree$edge[,2])){
   treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
 }

 return(treebis)
}


# Make a reconstruction using ape 
#(note that you can do directly the reconstruction using make.simmap, but because it is stochastic 
#you probably need to perform the model fitting on several reconstructed ancestral histories and average over them)
# Just a worked example
library(mvMORPH) # for plot tree + models 
#tree=pbtree(n=50)
#my_categories = rep(c("mysticetes","odontocetes"), each=25)
#trait <- rTraitCont(tree)

# the categories
#sta = as.character(my_categories)
#names(sta) = tree$tip.label

tree=subtree_unroot1
#range(subtree_unroot1$edge.length) - have a look at the rainge/tree type
data <- read.csv('sum_radius_sub_TEST.csv', row.names = 1)
data <- data[subtree_unroot1$tip.label, ]
sta = as.character(data$suborder)
names(sta) = tree$tip.label
trait <- data$sum.radii #add trait data 
names(trait) = tree$tip.label

#Add a very small length to the 0 branch length 
tree$edge.length[which(tree$edge.length == 0)] = 1e-5


## set zero-length branches to be 1/1000000 total tree length
#dst<-subtree_unroot1
#dst$edge.length[dst$edge.length==0]<-max(nodeHeights(subtree_unroot1))*1e-6
#fit3<-ace(sta,dst,type="discrete",model="ER")
#fit3

#Hypothesis1 
# ancestral state reconstruction (you can change the model to the complex ARD or the simple ER)
fitSYMTEST<-ace(sta, tree, model="ER", type="discrete")

# plot the tree states - spiral tree
plotTree(tree,type="fan",fsize=0.3,ftype="i",lwd=1)
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitSYMTEST$lik.anc, cex=0.2)

# transform it to simmap
tree2TEST = paintAllTree(tree, fitSYMTEST, sta)
cols = rainbow(ncol(tree2TEST$mapped.edge)); names(cols) = colnames(tree2TEST$mapped.edge)
plot(tree2TEST, col=cols)
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)


######hypothesis2 WITH REGIME #############################
sta2 = as.character(data$regime)
names(sta2) = tree$tip.label

# ancestral state reconstruction (you can change the model to the complex ARD or the simple ER)
fitSYM2<-ace(sta2, tree, model="ER", type="discrete")

# plot the tree states - spiral tree
plotTree(tree,type="fan",fsize=0.6,ftype="i",lwd=1) #can do type = 'fan' #lwd is line width 
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitSYM2$lik.anc, cex=0.2)

# transform it to simmap
tree3TEST = paintAllTree(tree, fitSYM2, sta2)
cols = rainbow(ncol(tree3TEST$mapped.edge)); names(cols) = colnames(tree3TEST$mapped.edge)
plot(tree3TEST, col=cols)
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)


##############################################################
#hypothesis3 - higher asymmetry in the odonts under different regimes 
sta3 = as.character(data$regime_split)
names(sta3) = tree$tip.label

# ancestral state reconstruction (you can change the model to the complex ARD or the simple ER)
fitSYM3<-ace(sta3, tree, model="ER", type="discrete")

# plot the tree states - spiral tree
plotTree(tree,fsize=0.6,ftype="i",lwd=1) #can do type = 'fan' #lwd is line width 
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM3$lik.anc, cex=0.2)

# transform it to simmap
tree4TEST = paintAllTree(tree, fitSYM3, sta3)
cols = rainbow(ncol(tree4TEST$mapped.edge)); names(cols) = colnames(tree4TEST$mapped.edge)
plot(tree4TEST, col=cols)
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)

#############################################################
#hypothesis4 - the presence and absence of echolocation 

sta4 = as.character(data$PA_echo)
names(sta4) = tree$tip.label

#ancestral state reconstruction (you can change the model to the complex ARD or the simple ER)
fitSYM4<-ace(sta4, tree, model="ER", type="discrete")

# plot the tree states - spiral tree
plotTree(tree,type="fan",fsize=0.6,ftype="i",lwd=1) #can do type = 'fan' #lwd is line width 
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM4$lik.anc, cex=0.2)

# transform it to simmap
tree5TEST = paintAllTree(tree, fitSYM4, sta4)
cols = rainbow(ncol(tree5TEST$mapped.edge)); names(cols) = colnames(tree5TEST$mapped.edge)
plot(tree5TEST, col=cols)
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)

##################################################################
######hypothesis5 WITH BAND (ECHOLOCATION)
sta5 = as.character(data$frequency)
names(sta5) = tree$tip.label

# ancestral state reconstruction (you can change the model to the complex ARD or the simple ER)
fitSYM5<-ace(sta5, tree, model="ER", type="discrete")

# plot the tree states - spiral tree
plotTree(tree,type="fan",fsize=0.6,ftype="i",lwd=1) #can do type = 'fan' #lwd is line width 
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM5$lik.anc, cex=0.2)

# transform it to simmap
tree6TEST = paintAllTree(tree, fitSYM5, sta5)
cols = rainbow(ncol(tree6TEST$mapped.edge)); names(cols) = colnames(tree6TEST$mapped.edge)
plot(tree6TEST, col=cols)
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)


######################################################################
#hypothesis4 - echolocation with better data e.g. hearing 
#sta4 = as.character(data$band_echo)
#names(sta4) = tree$tip.label

# ancestral state reconstruction (you can change the model to the complex ARD or the simple ER)
#fitSYM4<-ace(sta4, tree, model="ER", type="discrete")

# plot the tree states - spiral tree
#plotTree(tree,type="fan",fsize=0.6,ftype="i",lwd=1) #can do type = 'fan' #lwd is line width 
#nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitSYM4$lik.anc, cex=0.2)

# transform it to simmap
#tree5TEST = paintAllTree(tree, fitSYM4, sta4)
#cols = rainbow(ncol(tree5TEST$mapped.edge)); names(cols) = colnames(tree5TEST$mapped.edge)
#plot(tree5TEST, col=cols)
#add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)

###############################################################
# Model fit (here with mvMORPH but you can try OUwie)
fit1 <- mvBM(tree2TEST, trait, model="BM1")
fit2 <- mvBM(tree2TEST, trait, model="BMM")
fit3 <- mvBM(tree2TEST, trait, model="BMM", param=list(trend=TRUE))
fit4 <- mvBM(tree2TEST, trait, model="BMM", param=list(smean=FALSE))
fit5 <- mvBM(tree2TEST, trait, model="BM1", param=list(smean=FALSE))
fit6 <- mvOU(tree2TEST, trait, model="OU1")
fit7 <- mvOU(tree2TEST, trait, model="OUM")

#Regime model fits
fit8 <- mvBM(tree3TEST, trait, model="BMM")
fit9 <- mvBM(tree3TEST, trait, model="BMM", param=list(trend=TRUE))
fit10 <- mvBM(tree3TEST, trait, model="BMM", param=list(smean=FALSE))
fit11 <- mvOU(tree3TEST, trait, model="OUM")

#Regime split models 
fit12 <- mvBM(tree4TEST, trait, model="BMM")
fit13 <- mvBM(tree4TEST, trait, model="BMM", param=list(trend=TRUE))
fit14 <- mvBM(tree4TEST, trait, model="BMM", param=list(smean=FALSE))
fit15 <- mvOU(tree4TEST, trait, model="OUM")

#PA echo models 
fit16 <- mvBM(tree5TEST, trait, model="BMM")
fit17 <- mvBM(tree5TEST, trait, model="BMM", param=list(trend=TRUE))
fit18 <- mvBM(tree5TEST, trait, model="BMM", param=list(smean=FALSE))
fit19 <- mvOU(tree5TEST, trait, model="OUM")

#Frequency plots 
fit20 <- mvBM(tree6TEST, trait, model="BMM")
fit21 <- mvBM(tree6TEST, trait, model="BMM", param=list(trend=TRUE))
fit22 <- mvBM(tree6TEST, trait, model="BMM", param=list(smean=FALSE))
fit23 <- mvOU(tree6TEST, trait, model="OUM")

#PA echo plot 
#fit24 <- mvBM(tree7, trait, model="BMM")
#fit25 <- mvBM(tree7, trait, model="BMM", param=list(trend=TRUE))
#fit26 <- mvBM(tree7, trait, model="BMM", param=list(smean=FALSE))
#fit27 <- mvOU(tree7, trait, model="OUM")

# ? if you want rate dependence to climatic changes see fit_t_env in RPANDA, for the optimum ask Andy :slightly_smiling_face:
#Compare the fit of all of the above: 

results <- AIC(fit1);AIC(fit2);AIC(fit3);AIC(fit4);AIC(fit5);AIC(fit6);AIC(fit7);AIC(fit8);AIC(fit9);AIC(fit10);AIC(fit11);AIC(fit12);AIC(fit13);AIC(fit14);AIC(fit15);AIC(fit16);AIC(fit17);AIC(fit18);AIC(fit19);AIC(fit20);AIC(fit21);AIC(fit22);AIC(fit23)

#Extra models on Tree 3 (regime model)
#AIC(fit8); AIC(fit9); AIC(fit10); AIC(fit11)

#Extra models on Tree 4 (frequency bands model)
#AIC(fit12); AIC(fit13); AIC(fit14); AIC(fit15)

# Compute the Akaike weights:
results <- list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13,fit14,fit15,fit16,fit17,fit18,fit19,fit20,fit21,fit22,fit23)
weights <- aicw(results)
weights


#example make.simmap
#make simulations of different discrete traits 
tree3SIMMAP = make.simmap(tree, sta2, model="SYM", nsim = 10) #sta is the regime or model/variable
#par(mfrow=c(1,2)) #rows, columns
plot(tree3SIMMAP, type = "fan", col=cols, fsize = 0.5); #adding fan type 
axisPhylo()
multfit <- lapply(tree3SIMMAP, function(x) mvBM(x, trait, model="BMM"))


tree4SIMMAP = make.simmap(tree, sta3, model="SYM", nsim = 10) #sta is the regime or model/variable
#par(mfrow=c(1,2)) #rows, columns
plot(tree4SIMMAP, type = "fan", col=cols, fsize = 0.5);
axisPhylo()
multfit <- lapply(tree4SIMMAP, function(x) mvBM(x, trait, model="BMM"))


sigma_estimates <- sapply(multfit, function(x) x$sigma)
boxplot(t(sigma_estimates), col=cols)
AICweights = aicw(multfit)$aicweights
sigma_estimates*AICweights

#plot the overview of all the simulations 
pd <- summary(tree4SIMMAP)
pd


plot(pd,fsize=0.6,ftype="i")

