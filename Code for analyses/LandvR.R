
#Code adapted from Guillerme, T and Weisbecker V (2019). 
#landvR: Tools for measuring landmark position variation. Zenodo. doi:10.5281/zenodo.2620785

#Looking at asymmetry in odontocetes, mysticetes, archaeocetes, and terrestrial artiodactyls
#Run for all specimens (172) 
#Average specimen also 

#Code to run to pull out radii for assessing asymmetry in cetaceans

library(Rvcg)
library(rgl)
library(Morpho)
library(rgl)
library(geomorph)
library(paleomorph)

#Reading in landmarks only - these are .pts files
###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 

setwd("X:xxxxxxx/xxxxxx") 

#Mesh + spheres 
ply=ply2mesh(file="Delphinus delphis AMNH 75332 AB.ply")
shade3d(ply, col=bone1)
spheres3d(ptsarray[,,2],radius = 4,col=datcol2) 

#shows the colour spheres 
text3d(dataAB[,,2],text=1:dim(dataAB)[[1]]) #add text numbers 

#lollipop plot 
#Add label = TRUE if you want the numbers 
plotRefToTarget(ptsarray[,,1],MirroredAC[,,1],method="vector")


#=========== Combine datasets (AB and AC) to look at morphospace 


all_array <- abind(ptsarray, MirroredAC, along = 3)

Y.gpa3=gpagen(all_array) #Remove non-shape aspects 
data3=Y.gpa3$coords #Subset out the coords 
size=Y.gpa3$Csize

#plot morphospace
PCA=plotTangentSpace(all_array, axis1=1, axis2=2, label=dimnames(data)[[3]])



#========================================#
#      1. READ IN THE MANUAL LMS         #
#========================================#


###=========== LOADING DATA SET 1: WHOLE LANDMARKED SKULL 
#Read in landmarks manually placed on the whole skull 

ntaxa<-172 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarray<-array(dim=c(123,3,172)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Need the .plys for this 
#[3] stays the same 
dimnames(ptsarray)[3]<-list(
  substr(dir("./ply",pattern=".ply"),1,(nchar(dir("./ply",pattern=".ply"))-4)))
arraylm<-ptsarray #this is your array


##### MISSING LANDMARKS #########
#do this first and then rearrange landmarks 
#Landmarks need to be 'reshuffled' because we added 4 extra LMs on the nasal to be reflect that one midline landmark would 
#not work in the asymmetrical odontocetes 

arraylm[which(arraylm==9999)] <- NA
arraylm <- estimate.missing(arraylm,method="TPS")


##### MOVING NEW LANDMARKS INTO THE CORRECT POSITION ##############

#there are 123 landmarks with the last 4 added last (for asymmetry) 
#these last 4 landmarks had to be slotted in to the correct positions 
#120 ---> 67 
#121 ---> 70 
#122 ---> 85 
#123 ---> 85 

#Slot in the landmarks (see above)
#You can use abind to bind datasets eg. the odonts and mysts 
arranged=abind::abind(arraylm[c(1:66),,], arraylm[120,,],arraylm[c(67,68),,],
                      arraylm[121,,], arraylm[c(69:82),,], arraylm[c(122,123),,],
                      arraylm[c(83:119),,],
                      along= 1)

#show a print out with the new landmark numbers to check correct positioning 
#the last number [,,x] is the species number 
#check 
text3d(arranged[,,2], text=1:dim(arranged)[1])

#careful not to have missing data 
text3d(arranged[,,2], text=1:dim(arranged)[1])

#This would be to change the names according to the species (as above with .ply)
#dimnames(arranged)[3]=species

#############################
#                           #
#  PROCRUSTES THE DATA      #
#                           #
#############################

arranged=gpagen(arranged) #Remove non-shape aspects 
arranged_data=arranged$coords #Subset out the coords 
size=arranged$Csize #centroid size
#PlotTangentSpace if you want to see a quick morphospace of the skulls 

PCA <- plotTangentSpace(arranged$coords, legend=TRUE) 

#arranged_data is the Procrusted coords dataset of manually placed landmarks


############################################################
#                                                          #
#     2. Load second dataset - half LM skull to mirror     # 
#                                                          #                            
############################################################

setwd("X:xxxxxxx/xxxxxx") 
#NB the plys have been transformed to ASCII (from binary) so that landmarks can be visualised on the mesh 
#See code 'Binary_ASCII_ply.R" if needed for this step

ntaxa<-172 ## number of specimens (extant only) - NB can also put this in the code below (x,y,z) 
#data set .pts from Checkpoint

ptslist<-dir(pattern='.pts',recursive=T)
ptsarrayAC<-array(dim=c(66,3,34)) #dim=c(number of landmarks and semilandmarks, number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarrayAC[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}

#Set to the whole project directory for this - the code wants to look for 'ply' folder 
dimnames(ptsarrayAC)[3]<-list(
  substr(dir("./ply",pattern=".ply"),1,(nchar(dir("./ply",pattern=".ply"))-4)))###donne nom de scan a ptsarray et -4 pour retirer derniere lettre de nom
arraylmAC<-ptsarrayAC
arraylmAC


##### CHANGING MISSING LANDMARKS BEFORE MIRRORING #########
#Otherwise you get weird numbers that aren't 9999 mirroring
#estimate.missing from Geomorph is used for estimating missing landmarks 

arraylmAC[which(arraylmAC==9999)] <- NA
arraylmAC <- estimate.missing(arraylmAC,method="TPS")


#MIRROR THESE LANDMARKS over the central line plane of the skull 

########## SYMMETRISATION TO IMPROVE THE SHAPE ANALYSES #########################
#Ellen's code 

midline<-as.integer(c(38,40,48,49,51,54,55,56,61)) # LM that are on the midline + parasphenoid curve points + NO patch point
#got length(midline)= 9 points on the midline

left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

lengmatrice=dim(arraylmAC)[1]*2-length(midline)#-length(nasalfrontal) #should be the length with the both sides, 1 is the column and 2 
#just means that we are duplicating the data to be on both sides of the skull 

Matrice=array(NA,dim = c(lengmatrice,3,172)) #3 is the dimensions (x, y, z), 2 is specimen number 
Matrice[1:dim(arraylmAC)[1],,]=arraylmAC

#left.lm <- c(1:37,39,41:47,50,52,53,57:60,62:66)
#left.lm <- c(2,3,5:18,21:37,39,41:47,50,52,53,57:60,62:66)
#exclude midline points. Last number = last number of newpts 

#Check left.lm and midline [left.lm,,x] = species number
spheres3d(arraylmAC[left.lm,,10],radius=4) #left LMs
spheres3d(arraylmAC[midline,,10],radius=5,col='red') #midline

right.lm <- c(67:123) #left.lm +1:lenmatrice

bilat.landmarks <- cbind(left.lm, right.lm) #one column is the rows of the right side LM and the other column the rows of the left side

MirroredAC=mirrorfill(A=Matrice,  l1=midline, l2=bilat.landmarks) # the NA rows are now filled so that the numbers are the same on both
#sides of the skull 
MirroredAC
#deformGrid3d(MirroredAC[67:123,,2], Matrice[,,2], ngrid=0) #This shows you the new mirroed landmarks 

#These visualisations are done before Procrusted data 
#This shows the original landmarks
spheres3d(MirroredAC[c(1:66),,1],col=2,radius=4)
#check dimensions

#############################
#                           #
#  PROCRUSTES THE DATA      #
#                           #
#############################

MirroredAC=gpagen(MirroredAC) #Remove non-shape aspects 
Mirrored_data=Mirrored_data$coords #Subset out the coords 

#Checking if landmark ordering is correct 
#Comparing landmarks on both sides oft he skull 

col=rainbow(length(1:dim(Mirrored_data)[1]))
shapes3d(Y.gpa$coords[,,1], color=col)

#==========

if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

############# LANDVR ###################
#Check the values
MirroredAC_data[1,,2] == arranged_data[1,,2]
# FALSE FALSE FALSE
## This is due to rounding, in fact they have the same 9 digits - round them 
round(MirroredAC_data[1,,2], digits = 9) == round(arranged_data[1,,2], digits = 9)
# TRUE TRUE TRUE

#I’ve updated landvR to version 0.3 where the coordinates.difference function now have a tolerance optional argument.
#You can use the following to get the 0 difference results:

differences_between_lms <- coordinates.difference(coordinates = MirroredAC_data[,,65],
                                                  reference = arranged_data[,,65],
                                                  type = "spherical",
                                                  rounding = 9)

#Remove errornous missing landmarks (these should be zero because they are static)
differences_between_lms[[1]][1:66, 1:3] <- c(0.000000, 0.000000, 0.000000)


#Ellen's own colour function 
colfunc <- colorRampPalette(c("red", "yellow", "white"))
colfunc(10)
plot(rep(1,10),col=colfunc(10),pch=19,cex=3)

get.col.spectrum <- landvR::procrustes.var.plot(arranged_data[,,65], MirroredAC_data[,,65], col.val = differences_between_lms[[1]][,1], col = colfunc)


test=differences_between_lms[[1]][,1] #this is a test for specimen 1 to look at the differences between lms 
test

##### LOOKING AT AN AVERAGE SPECIMEN ######
N=123 #number of landmarks 
specs=166 #number of specimens 
all_combined=array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)

i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = Mirrored_data[,,i],
                                            reference = arranged_data[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}


#55, 56, 57, 59, 60 are all missing data and should be zero 
all_combined[1:66, 1:3] <- c(0.000000, 0.000000, 0.000000)
#write.csv(all_combined, file = 'all_combined.csv')

radii=all_combined[,1,] #looking at the second column (usually x,y,z) but here it is the radii, aziumuth, and polar 
radii_mean=apply(radii, c(1), mean) #c(1) looking at the first column which is the radii 
#test=all_combined[[1]][,,1] #this is a test for sepcimen 1 to look at the differences between lms 
#test

radii=all_combined[,1,] #second coloumn of whole dataset with just the radii [,1,]


#Looking at the average radii compared to specimen 21 (or an average specimen)
get.col.spectrum <- landvR::procrustes.var.plot(arranged_data[,,9], Mirrored_data[,,9], col.val = radii_mean, col = colfunc)
#datcol2<-c(rep("black",66),get.col.spectrum)
#open3d()

#############################
#                           #
#     MORPHOSPACE           #
#                           #
#############################

#Load species data .csv
#This should be a .csb with the specimen name, family, age, and any other data you want to plot by
Specimen_data_ALL=read.csv("D:\Checkpoint - ICVM\Species_data_asymmetry", header=TRUE, sep=",")
Specimen_data_ALL <- read.csv('Species_data_asymmetry_ALL.csv')

#Morphospace
PCA=plotTangentSpace(arranged_data, axis1=1, axis2=2, label = Specimen_data_ALL$species)

gp <- plotTangentSpace(arranged_data, 
                       groups = as.factor(paste(Specimen_data_ALL$suborder))) 

#PCA now has PC scores and variations 
library(ggplot2) #plot PCA
library(viridis)
library(ggfortify) # For prcomp 
#Manual plotting of data 
#To change the aethetics 
myplot = plotTangentSpace(arranged_data)
attributes(myplot) # shows extractable parts
PC.scores = myplot$pc.scores
df <- PC.scores
autoplot(prcomp(df)) #a rough and ready plot 

#Plotting suborder 
a <- autoplot(prcomp(df), data = Specimen_data_ALL, colour = 'suborder', size = 4)
a = a + scale_colour_manual(values=c("#D55E00", "#0072B2", "#73BFB8", "black")) 
a = a + theme_light()
a = a + labs(x="PC1 (38.5%)",y="PC2 (23.7%)", axes=FALSE, cex.lab = 10)
a = a + theme(axis.title.x = element_text(size = rel(1.15))) 
a = a + theme(axis.title.y = element_text(size = rel(1.15))) 
a = a +labs(color="suborder") ## change legend title
a


#Now read in the whole dataset 
#With lands 1-66 removed (we just want to look at the radii difference for the mirrored landmarks i.e. 67-123)
#This is for a PCA of the radii

library(factoextra)
#for prcomp (PCA)

#Compute PCA 
#Load landmark data 
landmarks <- read.csv('radii_X_PCA.csv')

# IF ERRORS #
#landmarks[,6] <- as.numeric(as.character(landmarks[,6]))
#omit.na(landmarks)
#any(is.na(landmarks)) #Want this to be FALSE 
#landmarks[is.na(landmarks)] <- 0.015850205 #Or whatever the value is 

myPr <- prcomp(landmarks[1:174,9:65], scale = TRUE) #pull out the columns of data which are nummerical only

#If the above isn't working, it's likely there is a numeric/factor problem
landmarks <- as.numeric(landmarks$X67)

plot(myPr, type = 'b') #boxplot to look at variation 
biplot(myPr, scale = 0) #explore the data 

#extract PC scores 
str(myPr)
myPr$x

summary(myPr) #look at PC variations 

Landradii <- cbind(landmarks, myPr$x[,1:2]) #pull out PC1 and PC2

library(ggplot2) #plot PCA
library(viridis)
library(ggfortify)
library(RColorBrewer)

#plot
a <- ggplot(Landradii, aes(PC1, PC2, col = suborder, fill = suborder, labels = TRUE)) + 
  geom_point(shape = 21) + 
  xlab('PC 1 (XX.X%)') + 
  ylab('XX.X%)') + 
  geom_point(aes(size = sum.radii)) +
  theme_classic() +
  #geom_text(aes(label= ID), hjust=1, vjust=2) - with numbers of specimens 
  #geom_text_repel(aes(label = test2$ID)) - with numbers of specimens repelled 
  
a 
#a + scale_x_reverse() #if you want to reverse axis 
#a + scale_y_reverse() #if you want to reverse axis
 
a + scale_colour_viridis_d() + 
  scale_x_reverse()#viridis

#variations of the morphospace 
a + scale_colour_viridis_d(option = "inferno") #plasma
#geom_text repel to repel text 

a + scale_colour_viridis_d(option = "plasma")
a + scale_colour_brewer(palette="YlOrRd")
a

library(vegan)
library(ggplot2)
library(ggConvexHull)
install.packages('ggConvexHull')

#Adding a convex hull 
ggplot(Landradii,aes(x = PC1, y = PC2,col = frequency)) +
  geom_convexhull(alpha = 0.3,aes(fill = frequency)) + 
  xlab('PC 1 (XX.X%)') + 
  ylab('PC2 (XX.X%)') +
  geom_point() +
  theme_minimal() + 
  scale_x_reverse()
  

