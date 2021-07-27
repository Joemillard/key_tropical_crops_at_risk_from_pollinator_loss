# Worldwide vulnerability of local pollinator abundance and crop pollination to land-use and climate change

This repository contains all the scripts used for the analysis carried out in the below paper, to be submitted:

> **Millard _et al_., (in prep). Worldwide vulnerability of local pollinator abundance and crop pollination to land-use and climate change.

Note that this repository is written relative to an R project file (.Rproj). Those wanting to reproduce this analysis should download the whole repo, and then open via the .Rproj file.

To run the code in this repo you will need to download a series of datasets:
1. The PREDICTS pollinator subset from either here (https://github.com/Joemillard/Global_effects_of_land-use_intensity_on_local_pollinator-biodiversity/tree/main/outputs) or here (https://figshare.com/articles/dataset/Global_effects_of_land-use_intensity_on_local_pollinator_biodiversity/12815669/2)
2. Pollination dependence ratios from here
3. Estimates of global crop production from here
4. Historical estimates of climate change from here
5. Projections of future climate change from here

The required script are as follows:

```R/00_functions.R```<br>
```R/01_map_climate_data.R```<br>
```R/02a_models_insect-pollinator_insect-non-pollinator.R```<br>
```R/02b_models_vertebrate-pollinator_vertebrate-non-pollinator.R```<br>
```R/02c_models_insect-pollinator_jack-knife.R```<br>
```R/03_map_pollinator_dependent_production.R```<br>
```R/04_map_proportional_production_risk.R```<br>
```R/05a_projections_production_risk.R```<br>
```R/05b_projections_production_risk_validation.R```<br>
```R/05c_projections_production_risk_data_quality.R```<br>
```R/06_projections_proportional_production_risk_country.R```<br>
```R/07_projections_proportional_production_risk_crop_cells.R```<br>
```R/08_projections_proportional_production_risk_crop_total.R```<br>