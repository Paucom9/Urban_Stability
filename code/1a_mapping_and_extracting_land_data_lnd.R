# ============================================================================ #
#   1a_mapping_and_extracting_land_data_lnd.R                               #
#   Author: Pau Colom                                                          #
#                                                                              #
#   This script prepares spatial data for the London study system.             #
#   It identifies well-sampled EBMS butterfly monitoring transects,            #
#   classifies them according to their position relative to the                #
#   Greater London boundary, and extracts landscape variables                  #
#   describing urbanization and habitat composition around each site.          #
#                                                                              #
#   Input*:                                                                    #
#   - EBMS transect coordinates:                                               #
#         input/data/ebms_transect_coord.csv                                   #
#   - EBMS visit data:                                                         #
#         input/data/ebms_visit.csv                                            #
#   - GHSL built-up surface rasters (3 arcsec and 30 arcsec):                  #
#         input/ghsl/                                                          #
#   - ESA WorldCover land cover tiles                                          #
#         input/land_cover/                                                    #
#   - London boundary:                                                         #
#         input/city_boundaries/london_boundary.gpkg                           #
#                                                                              #
#   Output:                                                                    #
#   - Built-up surface dataset:                                                #
#         output/data/built_up_lnd.csv                                         #
#   - Land cover diversity dataset:                                            #
#         output/data/land_diversity_lnd.csv                                   #
#   - Maps:                                                                    #
#         output/figures/urb_map_lnd.png                                       #
#         output/figures/landcover_map_lnd.png                                 #
#                                                                              #
#   *Note: Raw EBMS data are not publicly available in this repository.        #
#   Access requires a signed data-sharing agreement with the                   #
#   European Butterfly Monitoring Scheme (eBMS).                               #
#   Data requests can be submitted through:                                    #
#   https://butterfly-monitoring.net/                                          #
# ============================================================================ #



#1. Libraries ####

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(osmdata)
library(viridis)
library(ggspatial)
library(extrafont)
library(here)
library(exactextractr)


#2. Define project folders #####

input_dir  <- here("input")
output_dir <- here("output")


#3. Load EBMS data ####

m_coord <- read.csv(
  file.path(input_dir, "data", "ebms_transect_coord.csv")
)

m_visit <- read.csv(
  file.path(input_dir, "data", "ebms_visit.csv")
)


#4. Filter well-sampled sites ####

good_sites <- m_visit %>%
  
  filter(year >= 2017, year <= 2023) %>%
  
  group_by(bms_id, transect_id, year) %>%
  summarise(n_visits = n(), .groups = "drop") %>%
  
  filter(n_visits >= 10) %>%
  
  group_by(bms_id, transect_id) %>%
  summarise(years_with_10 = n(), .groups = "drop") %>%
  
  filter(years_with_10 >= 6)


m_coord_clean <- m_coord %>%
  
  inner_join(good_sites, by = c("bms_id", "transect_id")) %>%
  
  filter(!is.na(transect_lon), !is.na(transect_lat))


#5. Convert to spatial points ####

sites_sf <- st_as_sf(
  m_coord_clean,
  coords = c("transect_lon", "transect_lat"),
  crs = 3035
)

sites_sf <- st_transform(sites_sf, 4326)


#6. Get London boundary ####
# (download once from OSM or load from input folder)

boundary_file <- file.path(input_dir, "city_boundaries", "london_boundary.gpkg")

if (!file.exists(boundary_file)) {
  
  london_boundary <- opq("London") %>%
    add_osm_feature(key = "name", value = "Greater London") %>%
    osmdata_sf() %>%
    .$osm_multipolygons %>%
    slice(1)
  
  st_write(london_boundary, boundary_file, quiet = TRUE)
  
} else {
  
  london_boundary <- st_read(boundary_file, quiet = TRUE)
  
}


#7. Create London buffer (60 km) ####

lnd_buffer <- london_boundary %>%
  st_transform(27700) %>%
  st_buffer(60000) %>%
  st_transform(4326)



#8. Classify sites (inside / outside) ####

inside_flag <- st_within(
  sites_sf,
  london_boundary,
  sparse = FALSE
)[,1]

sites_sf$context <- ifelse(
  inside_flag,
  "inside",
  "outside"
)

sites_sf$context <- factor(
  sites_sf$context,
  levels = c("inside", "outside")
)

sites_sf$in_buffer <- st_within(
  sites_sf,
  lnd_buffer,
  sparse = FALSE
)[,1]

sites_buffer <- sites_sf[sites_sf$in_buffer, ]



#9. Load GHSL built-up rasters ####
# (30ss for extracting data; 3ss for plotting)

ghsl_dir <- file.path(input_dir, "ghsl")

# --- 30 arcsec resolution (analysis) ---
ghsl_30_files <- c(
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_V1_0_R4_C18.tif",
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_V1_0_R4_C19.tif"
)

ghsl_extract_30ss <- ghsl_30_files |>
  lapply(\(f) rast(file.path(ghsl_dir, f))) |>
  do.call(mosaic, args = _)

# --- 3 arcsec resolution (map plotting) ---
ghsl_3_files <- c(
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R4_C18.tif",
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R4_C19.tif"
)

ghsl_extract_3ss <- ghsl_3_files |>
  lapply(\(f) rast(file.path(ghsl_dir, f))) |>
  do.call(mosaic, args = _)


#10. Extract built-up surface values within different buffers ####
# (using 30 arcsec resolution)

built_values <- terra::extract(
  ghsl_extract_30ss,
  vect(sites_buffer)
)[,2]

pts_bng <- st_transform(sites_buffer, 27700)

results_df <- sites_buffer %>%
  st_drop_geometry() %>%
  dplyr::select(transect_id, context)

buffers <- c(1000, 2000, 5000)

na_counts <- numeric(length(buffers)) 

for (i in seq_along(buffers)) {
  r <- buffers[i]
  
  buf_bng <- st_buffer(pts_bng, dist = r)
  buf_wgs <- st_transform(buf_bng, crs(ghsl_extract_30ss))
  
  vals <- terra::extract(
    ghsl_extract_30ss,
    vect(buf_wgs),
    fun = mean,
    na.rm = TRUE
  )[,2]
  
  na_counts[i] <- sum(is.na(vals))
  
  if (any(is.na(vals))) {
    warning(sprintf("Radius %dm: %d NA values detected", r, na_counts[i]))
  }
  
  results_df[[paste0("BUILT_UP_", r, "m")]] <- vals
  
  setTxtProgressBar(pb, i)
}

built_up_df <- results_df %>%
  rename(
    SITE_ID = transect_id,
    CONTEXT = context,
    built1000 = BUILT_UP_1000m,
    built2000 = BUILT_UP_2000m,
    built5000 = BUILT_UP_5000m,
  )

write.csv(
  built_up_df,
  file.path(output_dir, "data", "built_up_lnd.csv"),
  row.names = FALSE
)


#11. Prepare GHSL raster for plotting ####
# (using 3 arcsec resolution)

ghsl_crop <- crop(
  ghsl_extract_3ss,
  vect(lnd_buffer)
)

ghsl_mask <- mask(
  ghsl_crop,
  vect(lnd_buffer)
)

ghsl_df <- as.data.frame(
  ghsl_mask,
  xy = TRUE
)

names(ghsl_df)[3] <- "built"

ghsl_df$built_log <- log1p(ghsl_df$built)


#12. Load land cover tiles ####

lc_files <- list.files(
  file.path(input_dir, "land_cover"),
  full.names = TRUE
)

lc_tiles <- lc_files[
  grepl("N51W003|N48W003|N51E000|N48E00", lc_files)
]

lc_raster_lnd <- do.call(mosaic, lapply(lc_tiles, rast))

lc_crop_lnd <- terra::crop(lc_raster_lnd, vect(lnd_buffer))

lc_mask_lnd <- terra::mask(lc_crop_lnd, vect(lnd_buffer))

lc_raster_low_lnd <- terra::aggregate(
  lc_mask_lnd,
  fact = 25,
  fun = function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
)

lc_df_lnd <- as.data.frame(lc_raster_low_lnd, xy = TRUE)
names(lc_df_lnd)[3] <- "class"


water_df_lnd <- lc_df_lnd %>%
  dplyr::filter(class == 80)


#13. Map: built-up surface ####

map_urb_lnd <- ggplot() +
  
  geom_raster(
    data = ghsl_df,
    aes(x = x, y = y, fill = built_log),
    alpha = 0.5
  ) +
  
  geom_sf(data = london_boundary, fill = NA, color = "black", size = 1.5) +
  
  geom_sf(data = sites_buffer, aes(color = context), size = 1) +
  
  geom_raster(data = water_df_lnd, aes(x = x, y = y), fill = "white", alpha = 1) +
  
  geom_sf(data = lnd_buffer, fill = NA, color = "black", size = 1.5) +
  
  scale_fill_viridis(option = "D", name = "Built-up surface") +
  
  scale_color_manual(
    values = c(inside = "red", outside = "black"),
    guide = "none"
  ) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_minimal(base_family = "Garamond", base_size = 18) +
  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  
  annotation_scale(location = "bl")

map_urb_lnd


ggsave(
  file.path(output_dir, "figures", "urb_map_lnd.png"),
  map_urb_lnd,
  width = 5.5,
  height = 5,
  dpi = 1200,
  bg = "white"
)


#14. Land cover map ####

lc_colors <- c(
  "10" = "#1a9850",
  "20" = "#66c2a5",
  "30" = "#cccc66",
  "40" = "#ffcc33",
  "50" = "#4D4D4D",
  "60" = "#ffff99",
  "80" = "white",
  "90" = "#2c7bb6"
)

labels_lc <- c(
  "10" = "Woodland",
  "20" = "Shrubland",
  "30" = "Grassland",
  "40" = "Cropland",
  "50" = "Built-up",
  "60" = "Sparse vegetation",
  "80" = "Water",
  "90" = "Wetland"
)

map_lc_lnd <- ggplot() +
  
  geom_raster(data = lc_df_lnd, aes(x = x, y = y, fill = factor(class)), alpha = 0.85) +
  
  
  geom_sf(data = lnd_buffer, fill = NA, color = "black", size = 1.2) +
  
  geom_sf(
    data = london_boundary,
    fill = "#e31a1c",
    alpha = 0.4,
    color = "#e31a1c",
    size = 2
  ) +
  
  geom_sf(
    data = subset(sites_buffer, context == "outside"),
    shape = 21,
    fill = "#0072B2",
    color = "black",
    size = 1.5
  ) +
  
  geom_sf(
    data = subset(sites_buffer, context == "inside"),
    shape = 21,
    fill = "#D55E00",
    color = "black",
    size = 1.5
  ) +
  
  scale_fill_manual(
    values = lc_colors,
    labels = labels_lc,
    name = "Land cover"
  ) +
  
  coord_sf(crs = 4326) +
  
  theme_minimal(base_family = "Garamond", base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.title = element_blank(),
  ) +
  
  annotation_scale(location = "bl")

map_lc_lnd


ggsave(
  file.path(output_dir, "figures", "lnd_landcover_map.png"),
  map_lc_lnd,
  width = 5.5,
  height = 5,
  dpi = 1200,
  bg = "white"
)

#15. Extract land cover diversity ####

pts_wgs <- st_transform(pts_bng, crs(lc_raster_lnd))

lc_raster_lnd

calc_habdiv <- function(raster, geom, valid_classes = c(10,20,30,40,60)) {
  cover_tab <- exact_extract(raster, geom, function(values, coverage_fraction) {
    tapply(coverage_fraction, values, sum, na.rm = TRUE)
  })
  
  # convertir a vector amb noms
  if (is.matrix(cover_tab)) {
    props <- cover_tab[,1]
    names(props) <- rownames(cover_tab)
  } else {
    props <- cover_tab
  }
  
  # filtrar i normalitzar
  names_num <- as.numeric(names(props))
  props <- props[names_num %in% valid_classes]
  props <- props / sum(props)
  
  if (length(props) == 0 || all(is.na(props))) return(0)
  
  vegan::diversity(props, index = "shannon")
}

for (r in buffers) {
  cat("\n=== Buffer:", r, "m ===\n")
  
  buf_bng <- st_buffer(pts_bng, dist = r)
  buf_wgs <- st_transform(buf_bng, crs(lc_raster_lnd))
  
  pb <- txtProgressBar(min = 0, max = nrow(buf_wgs), style = 3)
  
  habdiv_vals <- sapply(1:nrow(buf_wgs), function(i) {
    val <- calc_habdiv(lc_raster_lnd, buf_wgs[i, ])
    setTxtProgressBar(pb, i)
    cat("\nSite", i, "→ landdiv:", val, "\n")
    return(val)
  })
  
  close(pb)
  
  sites_buffer[[paste0("landdiv_", r, "m")]] <- habdiv_vals
}

landdiv_df <- sites_buffer %>%
  st_drop_geometry() %>%
  dplyr::select(transect_id, context, starts_with("landdiv_")) 


write.csv(
  landdiv_df,
  file.path(output_dir, "data", "land_diversity_lnd.csv"),
  row.names = FALSE
)

