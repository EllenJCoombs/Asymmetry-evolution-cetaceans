# Wonky whales: The evolution of cranial asymmetry in cetaceans
:whale2: :chart_with_upwards_trend:


__Authors:__
[Ellen J. Coombs](mailto:ellen.coombs.14@ucl.ac.uk),
[Julien Clavel](https://github.com/JClavel),
[Travis Park](https://github.com/travispark),
Morgan Churchill,
Anjali Goswami 

__To cite the paper__: NA
> 
Available at: 
If using any of this code or data please cite the paper above and the repo
__To cite this repo__: NA
> 
Github link + Zenodo
<insert Zenodo badge> 

## Data :bar_chart: :chart_with_downwards_trend: 

The data provided here are: 
1. the raw radii data for each specimen (n = 172); split by suborder, with details on ranking 
2. specimen data; ID, species, Family, geological age, echolocation frequency categories
3. Sum radii data for landmarks 67-123
4. Data with rostrum removed 

## Analysis :chart_with_upwards_trend: :whale2:
In this repository you will find raw data (.csv files), analyses (details on analyses), code for analyses (code supplied as .R files)
 :file_folder:
* **Raw data**
- CSV name and description 
- CSV name and description 
- CSV name and description 
 :file_folder:
* **Analyses**
- Details on analyses run e.g. links to papers, packages 
 :file_folder:
* **Code for analyses**
`Binary-ASCII-ply.R`
`LandvR.R`
`Build-phylogeny.R`
`Clavel-plotShifts.R`
`Clavel-models-shifts-jump.R`
`Discrete-trait-models.R`
`Model-diagnostics.R` 
`Phylo-MANOVAs`

## License :page_with_curl:
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/EllenJCoombs/Asymmetry-evolution-cetaceans/blob/master/LICENSE) file for details

## Session Info :clipboard:
For reproducibility purposes, here is the output of `devtools::session_info()` used to perform the analyses in the publication. 

```{r}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
[3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] phangorn_2.5.5          subplex_1.5-4           gtools_3.5.0            spam_2.2-2             
 [5] tidyselect_0.2.5        purrr_0.3.3             lattice_0.20-35         phytools_0.6-99        
 [9] vctrs_0.2.3             colorspace_1.3-2        expm_0.999-4            viridisLite_0.3.0      
[13] yaml_2.2.1              rlang_0.4.4             pillar_1.4.3            glue_1.3.1             
[17] ggfortify_0.4.5         lifecycle_0.1.0         stringr_1.4.0           dotCall64_1.0-0        
[21] munsell_0.5.0           combinat_0.0-8          gtable_0.2.0            coda_0.19-1            
[25] parallel_3.5.0          Rcpp_1.0.1              corpcor_1.6.9           scales_1.0.0           
[29] plotrix_3.7-6           clusterGeneration_1.3.4 scatterplot3d_0.3-41    mvMORPH_1.1.0          
[33] gridExtra_2.3           fastmatch_1.1-0         mnormt_1.5-5            ggplot2_3.2.1          
[37] stringi_1.1.7           ggrepel_0.8.1           animation_2.6           dplyr_0.8.4            
[41] numDeriv_2016.8-1.1     grid_3.5.0              quadprog_1.5-7          tools_3.5.0            
[45] magrittr_1.5            maps_3.3.0              lazyeval_0.2.1          tibble_2.1.3           
[49] tidyr_1.0.2             crayon_1.3.4            ape_5.3                 pkgconfig_2.0.3        
[53] MASS_7.3-49             Matrix_1.2-14           ggConvexHull_0.1.0      viridis_0.5.1          
[57] assertthat_0.2.0        rstudioapi_0.11         R6_2.2.2                glassoFast_1.0         
[61] igraph_1.2.4.1          nlme_3.1-137            compiler_3.5.0  
```
