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

-------------------------------------------------------------------------------------------------------------

  #GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_.tif
  #GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_.tif

Description
Global raster layers of built-up surface used to quantify urbanization around butterfly monitoring transects.
These rasters were used to calculate the amount of built-up area surrounding each transect within circular buffers of different radii.

Two spatial resolutions were used during processing:
3 arc-seconds (~100 m) resolution rasters for detailed spatial representation
30 arc-seconds (~1 km) resolution rasters used for efficient extraction during buffer calculations

Source
Built-up surface layers from the Global Human Settlement Layer (GHSL).

Dataset
GHSL Built-Up Surface 2020 (GHS-BUILT-S), Release 2023A.

Structure
Each file is a GeoTIFF raster tile in geographic coordinates (WGS84, EPSG:4326) representing the proportion of built-up surface within each grid cell.

Use in the analysis
These rasters were used to calculate built-up surface within 1000 m, 2000 m, and 5000 m buffers around butterfly monitoring transects. The resulting values are provided in the following datasets:
built_up_bcn.csv
built_up_lnd.csv
built_up_rnd.csv

The original GHSL rasters are publicly available from the European Commission Joint Research Centre.

-------------------------------------------------------------------------------------------------------------

  #ESA_WorldCover_10m_2021_*.tif

Description
Global land-cover raster tiles used to quantify landscape composition and habitat heterogeneity surrounding butterfly monitoring transects.
These layers correspond to the ESA WorldCover 2021 global land-cover product at 10 m spatial resolution.

Source
European Space Agency (ESA) WorldCover dataset.

Dataset
ESA WorldCover 10 m 2021 (Version 2.0).

Structure
Each GeoTIFF raster contains categorical land-cover classes representing different habitat types (e.g. forest, cropland, grassland, built-up areas, water bodies).

Use in the analysis
These rasters were used to compute landscape diversity around each butterfly monitoring transect by extracting land-cover composition within circular buffers of different radii.

-------------------------------------------------------------------------------------------------------------

  #lc_raster_low_bcn.tif
  #lc_raster_low_lnd.tif
  #lc_raster_low_rnd.tif

Description
Region-specific land-cover rasters derived from the ESA WorldCover dataset.

These rasters correspond to spatial subsets of the global ESA WorldCover dataset covering the three metropolitan study regions:

Barcelona metropolitan region (lc_raster_low_bcn)

London metropolitan region (lc_raster_low_lnd)

Randstad region (Netherlands) (lc_raster_low_rnd)

Use in the analysis
These rasters were used to calculate landscape diversity (Shannon diversity of land-cover classes) within 1000 m, 2000 m, and 5000 m buffers surrounding each monitoring transect. The resulting values are included in:

land_diversity_bcn.csv

land_diversity_lnd.csv

land_diversity_rnd.csv
-------------------------------------------------------------------------------------------------------------

  #phylogeny.nwk

Description
Phylogenetic tree of European butterfly species used to calculate phylogenetic diversity metrics.

Source
The phylogeny corresponds to the European butterfly phylogenetic hypothesis published by Wiemers et al. (2020).

Structure
Newick-format phylogenetic tree containing the butterfly species included in the analyses.

Use in the analysis
This phylogeny was used to calculate the Mean Phylogenetic Distance (MPD) metric for each butterfly community.

Reference
Wiemers, M., Chazot, N., Wheat, C. W., Schweiger, O., & Wahlberg, N. (2020). A complete time-calibrated multi-gene phylogeny of the European butterflies. ZooKeys, 938, 97.

-------------------------------------------------------------------------------------------------------------

  #species_trait_table.csv

Description
Species trait dataset used to calculate functional diversity metrics for butterfly communities.

Source
Trait data were derived from the comprehensive European butterfly trait database published by Middleton-Welling et al. 2020

Structure
Each row corresponds to a butterfly species and columns represent ecological and morphological traits used in the functional diversity analyses.

Use in the analysis
These traits were used to calculate Functional Dispersion (FDis) for each butterfly community.

Reference
Middleton-Welling, J., Dapporto, L., García-Barros, E., Wiemers, M., Nowicki, P., Plazio, E., et al. (2020). A new comprehensive trait database of European and Maghreb butterflies, Papilionoidea. Scientific Data, 7, 351.

-------------------------------------------------------------------------------------------------------------

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

  #land_diversity_bcn.csv
  #land_diversity_lnd.csv
  #land_diversity_rnd.csv

Description
Landscape diversity surrounding each butterfly monitoring transect derived from land-cover maps.

Landscape diversity was calculated within circular buffers of different radii around each transect centroid to quantify habitat heterogeneity at multiple spatial scales.

Source
ESA WorldCover land-cover dataset.

Structure
Each row corresponds to a monitoring transect.

| Variable    | Description                                                           |
| ----------- | --------------------------------------------------------------------- |
| SITE_ID     | Monitoring transect identifier                                        |
| CONTEXT     | Spatial context relative to the city boundary ("inside" or "outside") |
| landdiv1000 | Landscape diversity within a 1000 m radius buffer around the transect |
| landdiv2000 | Landscape diversity within a 2000 m radius buffer                     |
| landdiv5000 | Landscape diversity within a 5000 m radius buffer                     |


Landscape diversity values correspond to the Shannon diversity of land-cover classes within each buffer, calculated from ESA WorldCover land-cover maps. These values quantify habitat heterogeneity surrounding each monitoring transect.

-------------------------------------------------------------------------------------------------------------

  #diversity_stability_bcn.csv
  #diversity_stability_lnd.csv
  #diversity_stability_rnd.csv

Description
Community stability and biodiversity metrics calculated for butterfly communities at each monitoring transect. These variables were derived from annual butterfly abundance indices over the study period and quantify multiple components of community dynamics, including community stability, species asynchrony, population stability, and biodiversity.

These datasets provide the core ecological variables used in the structural equation models examining how biodiversity and landscape context influence the temporal stability of butterfly communities.

Source
Derived from standardized butterfly monitoring data from regional Butterfly Monitoring Schemes.

Structure
Each row corresponds to a monitoring transect.

| Variable                  | Description                                                                          |
| ------------------------- | ------------------------------------------------------------------------------------ |
| SITE_ID                   | Monitoring transect identifier                                                       |
| mean_community_abundance  | Mean total butterfly abundance across all years at the transect                      |
| sd_community_abundance    | Standard deviation of total community abundance across years                         |
| n_years                   | Number of years of monitoring data used to calculate community metrics               |
| community_stability       | Temporal stability of total community abundance (mean divided by standard deviation) |
| mean_population_abundance | Mean abundance of species populations across years                                   |
| mean_population_stability | Mean temporal stability of species populations                                       |
| wm_population_stability   | Abundance-weighted mean population stability                                         |
| species_richness          | Mean number of butterfly species recorded per year                                   |
| shannon_diversity         | Shannon diversity index of butterfly communities                                     |
| evenness_sw               | Community evenness derived from the Shannon index                                    |
| species_asynchrony        | Temporal asynchrony among species population dynamics                                |
| FDis                      | Functional dispersion of butterfly communities based on species traits               |
| MPD                       | Mean phylogenetic distance among species in the community                            |

These metrics summarize different mechanisms linking biodiversity to community stability, including population-level stability and compensatory dynamics among species.
Functional and phylogenetic diversity metrics were calculated from species trait data and phylogenetic relationships among European butterflies.

-------------------------------------------------------------------------------------------------------------

  #amsterdam_boundary.gpkg
  #barcelona_boundary.gpkg
  #london_boundary.gpkg
  #rotterdam_boundary.gpkg
  #thehague_boundary.gpkg
  #utrecht_boundary.gpkg

Description
Administrative boundaries of the cities used to define the spatial urban context of butterfly monitoring transects. These polygons were used to classify transects as located either inside or outside the urban core of each metropolitan region.
For the Randstad region, boundaries for the main cities (Amsterdam, Rotterdam, The Hague, and Utrecht) were combined to represent the metropolitan urban area.

Source
City administrative boundaries obtained from OpenStreetMap.

Structure
Each file contains a polygon representing the administrative boundary of a city.

Use in the analysis
Transect centroids were spatially intersected with these boundaries to classify sites into two spatial contexts:

inside: transects located within the city boundary
outside: transects located outside the city boundary

This classification was used to generate the CONTEXT variable included in the landscape datasets.
These boundaries were used exclusively to classify monitoring transects into urban and rural spatial contexts.


-------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------

#6. Reproducibility Notes

BMS data and annual abundance indices are not publicly available in this repository
Access requires a signed data-sharing agreement with the European Butterfly Monitoring Scheme (eBMS)                              
Data requests can be submitted through: https://butterfly-monitoring.net/     