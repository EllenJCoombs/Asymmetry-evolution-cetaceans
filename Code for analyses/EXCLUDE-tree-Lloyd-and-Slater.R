
#Running all analyses on Lloyd and Slater's 'EXCLUDE' tree - please use Lloyd and Slater reference if using these phylogenies: 
#Lloyd GT and Slater GJ (in prep). A total-group phylogenetic metatree of
#Cetacea and the importance of fossil data in diversification analyses

treeEXCLUDE <- read.tree('Exclude_map.tre')
#read in sum radii for specimens 
treeREMOVED <- read.nexus('treeREMOVED.nexus')

#tree EXCLUDE has loads of tips that aren't in my data - but I can't use these as we have no data on them 
#So...look at what's in REMOVED and drop the extra tips that aren't in EXCLUDE 
comparePhylo(treeEXCLUDE, treeREMOVED)

#Drop the tips that aren't in treeEXCLUDE
treeEXCLUDE_trim <- drop.tip(tree, c("Orycternocetus_crocodilinus", "Mesoplodon_traversii", "Berardius_minimus", 
"Sotalia_guianensis", "Globicephala_baereckii", "Kentriodon_fuchsii", "Zarhinocetus_donnamatsonae", "Schizodelphis_barnesi", 
"Schizodelphis_squalodontoides", "Patriocetid_new_genus_ChM_PV4753", "Patriocetus_sp_MB_Ma._42882", "Balaenoptera_floridana", 
"Balaenoptera_colcloughi", "Pakicetus_attocki", "Neophocaena_asiaeorientalis", "Orcaella_heinsohni", "Cephalorhynchus_hectori_maui", 
"Mesoplodon_hotaula", "Sousa_plumbea", "Sousa_teuszii", "Sousa_sahulensis")) ## drop species in tree but NOT from data

write.nexus(treeEXCLUDE_trim, file = 'treeEXCLUDE_trim.nexus')
plot(treeEXCLUDE_trim, cex = 0.6)
