# ============================================================================ #
#   1b_mapping_and_extracting_land_data_rnd.R                             #
#   Author: Pau Colom                                                          #
#                                                                              #
#   This script prepares spatial data for the Randstad study system.           #
#   It identifies well-sampled EBMS butterfly monitoring transects,            #
#   classifies them according to their position relative to the main           #
#   Randstad metropolitan cities (Amsterdam, Rotterdam, The Hague,             #
#   and Utrecht), and extracts landscape variables describing                  #
#   urbanization and habitat composition around each site.                     #
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
#   - City boundaries:                                                         #
#         input/city_boundaries/amsterdam_boundary.gpkg                        #
#         input/city_boundaries/rotterdam_boundary.gpkg                        #
#         input/city_boundaries/thehague_boundary.gpkg                         #
#         input/city_boundaries/utrecht_boundary.gpkg                          #
#                                                                              #
#   Output:                                                                    #
#   - Built-up surface dataset:                                                #
#         output/data/rnd_built_up.csv                                    #
#   - Land cover diversity dataset:                                            #
#         output/data/rnd_land_diversity.csv                              #
#   - Maps:                                                                    #
#         output/figures/rnd_urb_map.png                                  #
#         output/figures/rnd_landcover_map.png                            #
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


#2. Define project folders ####

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


#6. Get Randstad city boundaries ####
# (download once from OSM or load from input folder)

# 6. Get Randstad city boundaries ####
# (download once from OSM or load from input folder)

amsterdam_file <- file.path(input_dir, "city_boundaries", "amsterdam_boundary.gpkg")
rotterdam_file <- file.path(input_dir, "city_boundaries", "rotterdam_boundary.gpkg")
thehague_file  <- file.path(input_dir, "city_boundaries", "thehague_boundary.gpkg")
utrecht_file   <- file.path(input_dir, "city_boundaries", "utrecht_boundary.gpkg")

# --- Amsterdam ---
if (!file.exists(amsterdam_file)) {
  
  amsterdam_boundary <- opq("Amsterdam") %>%
    add_osm_feature(key = "name", value = "Amsterdam") %>%
    osmdata_sf() %>%
    .$osm_multipolygons %>%
    slice(1)
  
  st_write(amsterdam_boundary, amsterdam_file, quiet = TRUE)
  
} else {
  
  amsterdam_boundary <- st_read(amsterdam_file, quiet = TRUE)
}

# --- Rotterdam ---
if (!file.exists(rotterdam_file)) {
  
  rotterdam_boundary <- opq("Rotterdam") %>%
    add_osm_feature(key = "name", value = "Rotterdam") %>%
    osmdata_sf() %>%
    .$osm_multipolygons %>%
    slice(1)
  
  st_write(rotterdam_boundary, rotterdam_file, quiet = TRUE)
  
} else {
  
  rotterdam_boundary <- st_read(rotterdam_file, quiet = TRUE)
}

# --- The Hague ---
if (!file.exists(thehague_file)) {
  
  thehague_boundary <- opq("Den Haag") %>%
    add_osm_feature(key = "name", value = "Den Haag") %>%
    osmdata_sf() %>%
    .$osm_multipolygons %>%
    slice(1)
  
  st_write(thehague_boundary, thehague_file, quiet = TRUE)
  
} else {
  
  thehague_boundary <- st_read(thehague_file, quiet = TRUE)
}

# --- Utrecht ---
if (!file.exists(utrecht_file)) {
  
  utrecht_boundary <- opq("Utrecht") %>%
    add_osm_feature(key = "boundary", value = "administrative") %>%
    add_osm_feature(key = "admin_level", value = "8") %>%
    add_osm_feature(key = "name", value = "Utrecht") %>%
    osmdata_sf() %>%
    .$osm_multipolygons %>%
    slice(1)
  
  st_write(utrecht_boundary, utrecht_file, quiet = TRUE)
  
} else {
  
  utrecht_boundary <- st_read(utrecht_file, quiet = TRUE)
}

#7. Create Randstad buffer (60 km) ####

cities_geom <- st_sfc(
  st_geometry(st_transform(amsterdam_boundary, 28992))[[1]],
  st_geometry(st_transform(rotterdam_boundary, 28992))[[1]],
  st_geometry(st_transform(thehague_boundary, 28992))[[1]],
  st_geometry(st_transform(utrecht_boundary, 28992))[[1]],
  crs = 28992
)

cities_union <- st_union(st_combine(st_make_valid(cities_geom)))

rnd_buffer <- st_buffer(cities_union, dist = 60000) %>%
  st_transform(4326)

#8. Classify sites (inside / outside) ####

all_boundaries <- st_sfc(
  st_geometry(amsterdam_boundary)[[1]],
  st_geometry(rotterdam_boundary)[[1]],
  st_geometry(thehague_boundary)[[1]],
  st_geometry(utrecht_boundary)[[1]],
  crs = 4326
) %>% st_sf()

inside_flag <- st_within(
  sites_sf,
  all_boundaries,
  sparse = FALSE
)

sites_sf$context <- ifelse(
  apply(inside_flag, 1, any),
  "inside",
  "outside"
)

sites_sf$context <- factor(
  sites_sf$context,
  levels = c("inside", "outside")
)

sites_sf$in_buffer <- st_within(
  sites_sf,
  rnd_buffer,
  sparse = FALSE
)[,1]

sites_buffer <- sites_sf[sites_sf$in_buffer, ]

#9. Load GHSL built-up rasters ####
# (30ss for extracting data; 3ss for plotting)

ghsl_dir <- file.path(input_dir, "ghsl")

ghsl_3_files <- c(
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R4_C18.tif",
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R4_C19.tif"
)

ghsl_extract_3ss <- ghsl_3_files |>
  lapply(\(f) rast(file.path(ghsl_dir, f))) |>
  do.call(mosaic, args = _)


lc_files <- list.files(
  file.path(input_dir, "land_cover"),
  full.names = TRUE
)

#10. Load land cover tiles ####

lc_tiles <- lc_files[
  grepl("N51E000|N51E003|N51E006", lc_files)
]

lc_raster_rnd <- do.call(mosaic, lapply(lc_tiles, rast))

lc_crop_rnd <- terra::crop(lc_raster_rnd, vect(rnd_buffer))

lc_mask_rnd <- terra::mask(lc_crop_rnd, vect(rnd_buffer))

lc_raster_low_rnd <- terra::aggregate(
  lc_mask_rnd,
  fact = 25,
  fun = function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
)

lc_df_rnd <- as.data.frame(lc_raster_low_rnd, xy = TRUE)

names(lc_df_rnd)[3] <- "class"

water_df_rnd <- lc_df_rnd %>%
  dplyr::filter(class == 80)


#11. Prepare GHSL raster for plotting ####

ghsl_crop <- crop(
  ghsl_extract_3ss,
  vect(rnd_buffer)
)

ghsl_mask <- mask(
  ghsl_crop,
  vect(rnd_buffer)
)

ghsl_df <- as.data.frame(
  ghsl_mask,
  xy = TRUE
)

names(ghsl_df)[3] <- "built"

ghsl_df$built_log <- log1p(ghsl_df$built)


#12. Map built-up surface ####

map_urb_rnd <- ggplot() +
  
  geom_raster(
    data = ghsl_df,
    aes(x = x, y = y, fill = built_log),
    alpha = 0.5
  ) +
  
  geom_sf(data = all_boundaries, fill = NA, color = "black", size = 1.5) +
  
  geom_sf(data = sites_buffer, aes(color = context), size = 1) +
  
  geom_raster(
    data = water_df_rnd,
    aes(x = x, y = y),
    fill = "white",
    alpha = 1
  ) +
  
  geom_sf(data = rnd_buffer, fill = NA, color = "black", size = 1.5) +
  
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


ggsave(
  file.path(output_dir, "figures", "rnd_urb_map.png"),
  map_urb_rnd,
  width = 5.5,
  height = 5,
  dpi = 1200,
  bg = "white"
)


#13. Map land-cover ####

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

map_lc_rnd <- ggplot() +
  
  geom_raster(
    data = lc_df_rnd,
    aes(x = x, y = y, fill = factor(class)),
    alpha = 0.85
  ) +
  
  geom_sf(data = rnd_buffer, fill = NA, color = "black", size = 1.2) +
  
  geom_sf(
    data = all_boundaries,
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
    axis.title = element_blank()
  ) +
  
  annotation_scale(location = "bl")


ggsave(
  file.path(output_dir, "figures", "rnd_landcover_map.png"),
  map_lc_rnd,
  width = 5.5,
  height = 5,
  dpi = 1200,
  bg = "white"
)


#14. Extract land-cover diversity ####

pts_bng <- st_transform(sites_buffer, crs(lc_raster_rnd))


calc_habdiv <- function(raster, geom, valid_classes = c(10,20,30,40,60)) {
  
  cover_tab <- exact_extract(raster, geom, function(values, coverage_fraction) {
    tapply(coverage_fraction, values, sum, na.rm = TRUE)
  })
  
  if (is.matrix(cover_tab)) {
    props <- cover_tab[,1]
    names(props) <- rownames(cover_tab)
  } else {
    props <- cover_tab
  }
  
  names_num <- as.numeric(names(props))
  props <- props[names_num %in% valid_classes]
  props <- props / sum(props)
  
  if (length(props) == 0 || all(is.na(props))) return(0)
  
  vegan::diversity(props, index = "shannon")
}


buffers <- c(1000, 2000, 5000)


for (r in buffers) {
  
  buf_bng <- st_buffer(pts_bng, dist = r)
  buf_wgs <- st_transform(buf_bng, crs(lc_raster_rnd))
  
  pb <- txtProgressBar(min = 0, max = nrow(buf_wgs), style = 3)
  
  habdiv_vals <- sapply(1:nrow(buf_wgs), function(i) {
    
    val <- calc_habdiv(lc_raster_rnd, buf_wgs[i, ])
    
    setTxtProgressBar(pb, i)
    
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
  file.path(output_dir, "data", "rnd_land_diversity.csv"),
  row.names = FALSE
)

