# Worldwide vulnerability of local pollinator abundance and crop pollination to land-use and climate change

This repository contains all the scripts used for the analysis carried out in the below paper, to be submitted:

> Millard _et al_., (in prep). Worldwide vulnerability of local pollinator abundance and crop pollination to land-use and climate change.

Note that this repository is written relative to an R project file (.Rproj). Those wanting to reproduce this analysis should download the whole repo, and then open via the .Rproj file.

To run the code in this repo you will need to download a series of datasets:
1. The PREDICTS pollinator subset from either here (https://github.com/Joemillard/Global_effects_of_land-use_intensity_on_local_pollinator-biodiversity/tree/main/outputs) or here (https://figshare.com/articles/dataset/Global_effects_of_land-use_intensity_on_local_pollinator_biodiversity/12815669/2)
2. Klein et al (2007) pollination dependence ratios from this GitHub repo in the folder (data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/)
3. Estimates of global crop production from here (http://www.earthstat.org/harvested-area-yield-175-crops/)
4. Historical estimates of climate change from here, same as in the previous paper but here just tmp (https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/cruts.1905011326.v4.03/)
5. ISIMIP anomaly projections of future climate change (folder is ISIMIPAnomalies.tar.gz) from onedrive shared by Tim here (https://liveuclac-my.sharepoint.com/:f:/r/personal/ucbttne_ucl_ac_uk/Documents/LandUseClimate/FutureClimate?csf=1&web=1&e=A1ihT0)

The required scripts are as below. All 'models' scripts introduce the models built for the interaction of climate change and land use. All 'map' scripts represent either the current or future geographic distribution of climate change, crop production, or proportional production risk. All 'projection' scripts project future change over time in either proportional production risk or production risk. 

```R/00_functions.R```<br>
```R/01_map_climate_data.R``` -- Figure 2<br> 
```R/02a_models_insect-pollinator_insect-non-pollinator.R``` -- Figure 1<br>
```R/02b_models_vertebrate-pollinator_vertebrate-non-pollinator.R``` -- Figure 1<br>
```R/02c_models_insect-pollinator_jack-knife.R``` -- Figure S5<br>
```R/03_map_pollinator_dependent_production.R``` -- Figure 2<br>
```R/04_map_proportional_production_risk.R``` -- Figure 4<br>
```R/05a_projections_production_risk.R``` -- Figure 3<br>
```R/05b_projections_production_risk_validation.R``` -- Figure S2<br>
```R/05c_projections_production_risk_data_quality.R``` -- Figure S4<br>
```R/06_projections_proportional_production_risk_country.R``` -- Figure 4<br>
```R/07_projections_proportional_production_risk_crop_cells.R``` -- Figure 5<br>
```R/08_projections_proportional_production_risk_crop_total.R``` -- Figure S3<br>