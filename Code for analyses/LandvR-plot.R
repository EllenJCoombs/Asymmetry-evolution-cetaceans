

#if(!require(devtools)) install.packages("devtools")
#library(devtools)
#install_github("TGuillerme/landvR")
library(landvR)

#Guillerme T, Weisbecker V. landvR: tools for measuring landmark position
#variation. 2019. Zenodo. . https://doi.org/10.5281/zenodo.2620785

setwd("C:/Users/carlb/Box Sync/Carla Bardua/Ellen_whales")

load("./arranged_data.R")
load("./MirroredAC_data.R")

##first, find the coord differences for all specimens
N=123 #number of landmarks 
specs=172 #number of specimens 
all_combined=array(dim=c(N,3,specs)) #3 is the columns of data we need (radii, azimuth, polar)
i=1
for (i in 1:specs)
{
  all_differences <- coordinates.difference(coordinates = MirroredAC_data[,,i],
                                            reference = arranged_data[,,i],
                                            type = "spherical",
                                            rounding = 9)
  
  all_combined[,,i]=all_differences[[1]]
  
  i=i+1
}
#1-60 are all fixed landmarks make sure they are zero 
all_combined[1:66, 1:3, 1:172] <- c(0.000000, 0.000000, 0.000000)

species_data <- read.csv("./Species_data_asymmetry_REMOVED_ALL.csv")
odontocetes=which(species_data$suborder=="odontocete")
archaeocete=which(species_data$suborder=="archaeocete")
mysticete=which(species_data$suborder=="mysticete")

##now for each group, find the mean of the arranged data, mirrored data and coord diffs
myst_arranged_data=arranged_data[,,mysticete]
myst_arranged_data_mean=apply(myst_arranged_data, c(1,2), mean)
myst_MirroredAC_data=MirroredAC_data[,,mysticete]
myst_MirroredAC_data_mean=apply(myst_MirroredAC_data, c(1,2), mean)
myst_all_combined=all_combined[,,mysticete]
myst_all_combined_mean=apply(myst_all_combined, c(1,2), mean)

odont_arranged_data=arranged_data[,,odontocetes]
odont_arranged_data_mean=apply(odont_arranged_data, c(1,2), mean)
odont_MirroredAC_data=MirroredAC_data[,,odontocetes]
odont_MirroredAC_data_mean=apply(odont_MirroredAC_data, c(1,2), mean)
odont_all_combined=all_combined[,,odontocetes]
odont_all_combined_mean=apply(odont_all_combined, c(1,2), mean)

arch_arranged_data=arranged_data[,,archaeocete]
arch_arranged_data_mean=apply(arch_arranged_data, c(1,2), mean)
arch_MirroredAC_data=MirroredAC_data[,,archaeocete]
arch_MirroredAC_data_mean=apply(arch_MirroredAC_data, c(1,2), mean)
arch_all_combined=all_combined[,,archaeocete]
arch_all_combined_mean=apply(arch_all_combined, c(1,2), mean)

all_means=abind::abind(arch_all_combined_mean, odont_all_combined_mean, myst_all_combined_mean, 
                       along=3)
dim(all_means) ##should be 123,3,3 as 3 means

radius = 1 # i.e 1 = radius, 2 = azimuth, 3 = polar
min(all_means[,radius,]) 
max(all_means[,radius,]) 
##use these to manually input range of b, round up a little for max value to make sure it's covered by the range
ii <- cut(all_means[,radius,], b = 0:0.0125, breaks = 100,
          include.lowest = TRUE)

length(ii) ## e.g. length of 3 groups would be 3*123=369

##assign colours
n_cols=100 ## 100 as 100 breaks- can change this
colors <-  colorRampPalette(c("white","yellow", "red"))(n_cols)[ii]

arch_cols=colors[c(1:123)]
odont_cols=colors[c(124:246)]
myst_cols=colors[c(247:369)]

get.col.spectrum <- landvR::procrustes.var.plot(arch_arranged_data_mean, 
                                                arch_MirroredAC_data_mean, 
                                                col.val = arch_all_combined_mean[,radius], 
                                                col = arch_cols)

get.col.spectrum <- landvR::procrustes.var.plot(odont_arranged_data_mean, 
                                                odont_MirroredAC_data_mean, 
                                                col.val = odont_all_combined_mean[,radius], 
                                                col = odont_cols)

get.col.spectrum <- landvR::procrustes.var.plot(myst_arranged_data_mean, 
                                                myst_MirroredAC_data_mean, 
                                                col.val = myst_all_combined_mean[,radius], 
                                                col = myst_cols)
