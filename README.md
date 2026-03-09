# Urban Stability Project – Data and Code Repository

#1. Overview

This repository contains the data and code required to reproduce the analyses presented in:

“Urbanization weakens biodiversity–stability relationships in butterfly communities.”
Colom, P., García-Callejas, D., Schmucki, R., Stefanescu, C., van Swaay, C. & Melero, Y.

The study investigates how urbanization affects the relationship between diversity and ecological stability in butterfly communities across three European metropolitan regions:

- Barcelona (Catalonia, Spain)
- London (United Kingdom)
- Randstad (The Netherlands)

Using standardized butterfly monitoring data, we test how urban landscape context influences diversity, species asynchrony, population stability, and overall community stability.

-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

#2. Project Structure

Urban_Stability
│
├── code/
│   ├── 1a_mapping_and_extracting_land_data_lnd.R
│   ├── 1b_mapping_and_extracting_land_data_rnd.R
│   ├── 1c_mapping_and_extracting_land_data_bcn.R
│   ├── 2_annual_abundance_indices.R
│   ├── 3_diversity_stability_metrics.R
│   └── 4_stability_models.R
│
├── input/
│   ├── city_boundaries/
│   ├── data/
│   ├── ghsl/
│   └── land_cover/
│
├── output/
│   ├── data/
│   ├── figures/
│   └── results/
│
├── Urban_Stability.Rproj
├── .gitignore
└── README.md

-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

#3. Scripts Description and Order of Execution

  #1a_mapping_and_extracting_land_data_lnd.R
  #1b_mapping_and_extracting_land_data_rnd.R
  #1c_mapping_and_extracting_land_data_bcn.R

Purpose
Prepare spatial datasets for each metropolitan region (London, Randstad, Barcelona).

Main tasks
Filter well-sampled butterfly monitoring transects
Classify transects as urban (inside city boundaries) or rural (outside)
Extract landscape variables around each site
Compute built-up surface metrics
Compute land-cover diversity
Produce spatial maps used in the manuscript

Input
input/data/ebms_transect_coord.csv
input/data/ebms_visit.csv
input/ghsl/*.tif
input/land_cover/*.tif
input/city_boundaries/*.gpkg

Output
output/data/built_up_bcn.csv
output/data/built_up_lnd.csv
output/data/built_up_rnd.csv
output/data/land_diversity_bcn.csv
output/data/land_diversity_lnd.csv
output/data/land_diversity_rnd.csv

Maps produced
output/figures/urb_map_bcn.png
output/figures/urb_map_lnd.png
output/figures/urb_map_rnd.png

output/figures/landcover_map_bcn.png
output/figures/landcover_map_lnd.png
output/figures/landcover_map_rnd.png

-------------------------------------------------------------------------------------------------------------

  #2_annual_abundance_indices.R

Purpose
Compute annual butterfly abundance indices for each monitoring transect.

Input
input/data/ebms_count.csv
input/data/ebms_visit.csv
input/data/coord_length_transects_CBMS.csv
input/data/ubms_visits.csv

Output
output/data/annual_abundance_bcn.csv
output/data/annual_abundance_lnd.csv
output/data/annual_abundance_rnd.csv

Each row corresponds to a site × species × year combination.

-------------------------------------------------------------------------------------------------------------

  #3_diversity_stability_metrics.R

Purpose
Calculate biodiversity and stability metrics for each site.

Metrics calculated
Species richness
Shannon diversity
Phylogenetic diversity
Functional diversity
Community stability
Population stability
Species asynchrony

Input
output/data/annual_abundance_*.csv
input/data/species_trait_table.csv
input/data/phylogeny.nwk

Output
output/data/diversity_stability_bcn.csv
output/data/diversity_stability_lnd.csv
output/data/diversity_stability_rnd.csv

These datasets contain the core ecological metrics used in the models.

-------------------------------------------------------------------------------------------------------------

  #4_stability_models.R

Purpose
Test the effect of biodiversity and urbanization on community stability.

Analyses performed
Linear models relating diversity and stability
Generalized linear models testing the effect of urban context
Structural Equation Models (SEM) evaluating causal pathways linking:
Urbanization
Biodiversity
Species asynchrony
Population stability
Community stability

Input
output/data/diversity_stability_bcn.csv
output/data/diversity_stability_lnd.csv
output/data/diversity_stability_rnd.csv

Output
output/results/GLM_diversity_stability_context.csv
output/results/SEM_paths_all_regions_buffers.csv

Figures produced
output/figures/fig_diversity_stability_context.png
output/figures/fig_diversity_stability_context_region.png

-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

#4. Relevant session info

The analyses were conducted using the following R version and package versions:

R version 4.5.1 (2025-06-13 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Package versions
                    package  version
sf                       sf   1.0.21
terra                 terra   1.8.54
dplyr                 dplyr    1.1.4
ggplot2             ggplot2    3.5.2
osmdata             osmdata    0.2.5
viridis             viridis    0.6.5
ggspatial         ggspatial   1.1.10
rbms                   rbms    1.2.0
lubridate         lubridate    1.9.4
here                   here    1.0.2
extrafont         extrafont     0.19
exactextractr exactextractr   0.10.0
tidyr                 tidyr    1.3.1
purrr                 purrr    1.2.0
vegan                 vegan    2.7.1
FD                       FD 1.0.12.3
picante             picante    1.8.2
ape                     ape    5.8.1
stringr             stringr    1.5.1
tibble               tibble    3.3.0
lavaan               lavaan   0.6.19
lavaanPlot       lavaanPlot    0.8.1
lme4                   lme4   1.1.37
lmerTest           lmerTest    3.1.3
ggeffects         ggeffects    2.3.0
patchwork         patchwork    1.3.1
broom                 broom    1.0.8

-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

#5. Data Description

This section describes the datasets included in the repository that are necessary to reproduce the statistical analyses presented in the manuscript.
Raw butterfly monitoring data and annual abundance indices are not included due to data-sharing restrictions (see section 6 - Reproducibility Notes)

  #built_up_bcn.csv
  #built_up_lnd.csv
  #built_up_rnd.csv

Description
Built-up surface surrounding each butterfly monitoring transect extracted from the Global Human Settlement Layer (GHSL).

Built-up area was calculated within circular buffers of different radii around each transect centroid to quantify the degree of urbanization at multiple spatial scales.

Source
Global Human Settlement Layer (GHSL).

Structure
Each row corresponds to a monitoring transect.

| Variable  | Description                                                           |
| --------- | --------------------------------------------------------------------- |
| SITE_ID   | Monitoring transect identifier                                        |
| CONTEXT   | Spatial context relative to the city boundary ("inside" or "outside") |
| built1000 | Built-up surface within a 1000 m radius buffer around the transect    |
| built2000 | Built-up surface within a 2000 m radius buffer                        |
| built5000 | Built-up surface within a 5000 m radius buffer                        |


Built-up surface values correspond to the total area classified as built-up within each buffer based on GHSL built-up layers.

-------------------------------------------------------------------------------------------------------------



-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

#6. Reproducibility Notes

BMS data and annual abundance indices are not publicly available in this repository
Access requires a signed data-sharing agreement with the European Butterfly Monitoring Scheme (eBMS)                              
Data requests can be submitted through: https://butterfly-monitoring.net/     