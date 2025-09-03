# Arctic-cod-gillraker-linked-diet

Scripts repository used to analyse Arctic cod gill rakers and their link to diet and habitat preference. This repository contains the relevant scripts needed to analyse the data files for the manuscript titled:

## Assessing the link between ontogenetic shifts in gill raker and diet composition of Arctic cod (Boreogadus saida) in the Canadian Beaufort Sea-Amundsen Gulf

Authored by: Marie Launay, Billy Plaitis, Andrew R. Majewski, Andrea Niemi, James D. Reist, María Quintela, Torild Johansen, and Denis Roy

#
### Scripts:

To see the main scripts for this research, go to the 'Scripts' folder. Script description:

* Map.R: This script uses ggOceanMap, written by Mikko Vihtakari from the IMR, Tromsø, Norway, to map out the sampled location in the Beaufort Sea (Vihtakari M (2024). _ggOceanMaps: Plot Data on Oceanographic Maps using 'ggplot2'_. R package version 2.2.0, <https://CRAN.R-project.org/package=ggOceanMaps>). The basemap is uploaded, and the coordinates are loaded through a .csv file that contains both the lat/long in decimal degrees, but also the lat/long of where the supporting features should be drawn on the map. See the citation for more information on its use.

* ACPower.R: This script is used to estimate the sample size required to see a given effect size. It then plots the results using ggplot to visualise the relationship between sample size and estimated effect size at a given significance level (here set to alpha = 0.05).

* AC_age_gr_dietsize.R: Script that reads in the Arctic cod size, Gillraker and diet item sizes as measured during dissection under the stereomicroscope in the lab. Arctic cod size is converted to age using the von Bertalanffy growth model with the relationship and parameters derived by Forster et al. 2021 (as per 
Malizia et al. 2023). Gill raker (GR) data is used to calculate GR density, and diet item sizes are used to estimate the mean diet item sizes for each individual.

* AC_diet_item.R: A script that reads in Arctic cod diet data compiled during diet analyses. Data are listed as the sizes of 18 different possible diet categories recovered from the stomachs of each individual. Each item was measured to get diet sizes, but was also counted. So entries are on one or more single rows per individual. The data thus have to be compiled by individual, and then assessed using the points method as used in Genner et al. 2003 and in Roy et al. 2007.

* AC_SIA.R: This script reads in Script that reads in stable isotope data determined from individual Arctic cod collected from the Beaufort Sea and assesses them for differences among the estimated age classes (determined from their length and calculated using the von Bertalanffy growth equation). The data are plotted and tested for adherence to normality and, when they do used in linear mixed effect models testing differences among age classes but also using location as a random factor. Results of among age class differences (and subsequent location differences) are then tested in pairwise assessments. Plotted data have 95% ellipses plotted, and these are used to calculate the ellipse eccentricity. Eccentricity gives an estimate of diet specificity or constraints relative to the isotopes used to describe it. The last part of the script uses SIBER (Stable Isotopes Bayesian Ellipses in R) to redo some of the basics – but importantly uses a Bayesian framework to assess diet overlap among age classes. The last bit takes some time to run – even on a relatively fast machine.

## Note:
Scripts written here are to help with reproducibility and to be used with the data avaiable from the Borealis data repository at https://doi.org/10.5683/SP3/0RKRGH 



