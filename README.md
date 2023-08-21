# Key tropical crops at risk from pollinator loss due to climate change and land use

### Introduction

This repository contains all the scripts used for the analysis carried out in the below paper:

> Millard _et al_., (in review). Key tropical crops at risk from pollinator loss due to climate change and land use

Note that this repository is written relative to an R project file (.Rproj). Those wanting to reproduce this analysis should download the whole repo, and then open via the .Rproj file.

To run the code in this repo you will need to download a series of datasets:
1. The PREDICTS pollinator subset from here (https://doi.org/10.5281/zenodo.7385950)
2. Klein et al (2007) pollination dependence ratios from this GitHub repo in the folder (https://github.com/Joemillard/Worldwide-vulnerability-of-local-pollinator-abundance-and-crop-pollination-to-land-use-and-climate/tree/master/data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops)
3. Estimates of global crop production from Monfreda et al (2008) (see here to download http://www.earthstat.org/harvested-area-yield-175-crops/)
4. Historical estimates of climate change from Harris et al (2020) (see here to download https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/cruts.1905011326.v4.03/)
5. ISIMIP anomaly projections of future climate change from Warszawski et al (2014)
6. Esimates of global crop prices from the Food and Agriculture Organization of the United Nations (FAO) (see https://www.fao.org/faostat/en/#data/PP to download)
7. Virtual pollination trade flow data from Silva et al (2021) (see here to download https://github.com/virtual-pollination-trade/virtual-biotic-pollination-flow)
8. Estimates of global population in the years 2015-2019, taken from Our World in Data (see https://ourworldindata.org/population-sources)

------------

### Structure

The required scripts are as below, with figure numbers from the revised submission on 21/08/2022. All 'models' scripts introduce the models built for the interaction of climate change and land use. All 'map' scripts represent either the current or future geographic distribution of climate change, crop production, or proportional production risk. All 'projection' scripts project future change over time in either proportional production risk or production risk. 

- ```R/00_functions.R```                                                      
- ```R/01_map_climate_data.R```                                               
- ```R/01a_map_predicts_sites.R```                                             
- ```R/01b_predicts_descriptive.R```                                           
- ```R/02a_models_insect-pollinator_insect-non-pollinator.R```            
- ```R/02c_models_insect-pollinator_jack-knife.R```                         
- ```R/02d_models_insect-pollinator_family_jack-knife.R```                    
- ```R/02e_models_insect-pollinator_active_season_validation_baseline_site.R```
- ```R/02f_models_insect_pollinator_active_season_validation_PREDICTS_site.R```
- ```R/02g_models_insect-pollinator_insect-non-pollinator_mean_annual_temp.R```
- ```R/02h_models_insect-pollinator_insect-non-pollinator_sampling method.R```
- ```R/03_map_pollinator_dependent_production.R```                             
- ```R/04_map_proportional_production_risk.R```                            
- ```R/04a_map_proportional_production_risk_standard_dev.R```                
- ```R/04b_map_proportional_production_risk_richness.R```                   
- ```R/04c_map_proportional_production_risk_bees.R```                         
- ```R/05a_projections_production_risk.R```                                   
- ```R/05b_projections_production_risk_validation.R```                        
- ```R/05c_projections_production_risk_data_quality.R```                      
- ```R/05d_projections_production_risk_abundance-service_validation.R```      
- ```R/05e_projections_production_risk_abundance-service_main-figure.R```    
- ```R/06_projections_proportional_production_risk_country.R```               
- ```R/06a_projections_total_production_risk_country.R```                      
- ```R/07_projections_proportional_production_risk_crop_cells.R```            
- ```R/08_projections_proportional_production_risk_crop_total.R```           
- ```R/09_trade_flow_formatting.R```                                       
- ```R/10_projections_import_production_risk.R```

------------

### References
Food and Agriculture Organization of the United Nations (FAO). Producer prices & crops and livestock products. (2022).

Harris, I., Osborn, T. J., Jones, P. & Lister, D. 2020. Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. Sci. Data 7, 109 

Klein, A.M., Vaissière, B.E., Cane, J.H., Steffan-Dewenter, I., Cunningham, S.A., Kremen, C. and Tscharntke, T., 2007. Importance of pollinators in changing landscapes for world crops. Proceedings of the royal society B: biological sciences, 274(1608), pp.303-313.

Mathieu, E., Rodés-Guirao, L., Roser, M., What are the sources for Our World in Data’s population estimates? Our World in Data (2023) (available at https://ourworldindata.org/population-sources).

Monfreda, C., Ramankutty, N. and Foley, J.A., 2008. Farming the planet: 2. Geographic distribution of crop areas, yields, physiological types, and net primary production in the year 2000. Global biogeochemical cycles, 22(1).

Silva, F. D. S., Carvalheiro, L. G. , Aguirre-Gutiérrez, J., Lucotte, M., Guidoni-Martins, K., Mertens, F. Virtual pollination trade uncovers global dependence on biodiversity of developing countries. Science Advances. 7, eabe6636 (2021).

Warszawski, L. et al. 2014. The inter-sectoral impact model intercomparison project (ISI-MIP): project framework. Proc. Natl Acad. Sci. USA 111, 3228–3232

------------

### Acknowledgements

Thanks to Richard Gregory, Robin Freeman, and colleagues on the GLITRS project for useful discussions and comments. Thanks to Andy Purvis and Melinda Mills, who both provided support and guidance. Thanks also to Opeyemi Adedoja, Sabrina Gavini, Esther Kioko, Michael Kuhlmann, Zong-Xin Ren, and Manu Saunders, who provided input on our set of likely pollinating species, and to all contributors of the PREDICTS database.

------------

### License

MIT License, Copyright (c) 2023 Joe Millard

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.