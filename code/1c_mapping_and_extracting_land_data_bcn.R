# ============================================================================ #
#   1c_mapping_and_extracting_land_data_bcn.R                                  #
#   Author: Pau Colom                                                          #
#                                                                              #
#   This script prepares spatial data for the Barcelona study system.          #
#   It merges CBMS and uBMS transects, identifies well-sampled sites,          #
#   classifies them according to their position relative to the                #
#   Barcelona boundary, and extracts landscape variables describing            #
#   urbanization and habitat composition around each site.                     #
#                                                                              #
#   Input*:                                                                    #
#   - Monitoring datasets:                                                     #
#         input/data/s_index_18_24_ubms_cbms.csv                               #
#         input/data/coord_length_transects_CBMS.csv                           #
#         input/data/ubms_visits.csv                                           #
#         input/data/ubms_sites.csv                                            #
#                                                                              #
#   - GHSL built-up surface rasters                                            #
#         input/ghsl/                                                          #
#                                                                              #
#   - ESA WorldCover land cover tiles                                          #
#         input/land_cover/                                                    #
#                                                                              #
#   - Barcelona boundary                                                       #
#         input/city_boundaries/barcelona_boundary.gpkg                        #
#                                                                              #
#   Output:                                                                    #
#   - Built-up surface dataset                                                 #
#         output/data/built_up_bcn.csv                                         #
#                                                                              #
#   - Land cover diversity dataset                                             #
#         output/data/land_diversity_bcn.csv                                   #
#                                                                              #
#   - Maps                                                                     #
#         output/figures/urb_map_bcn.png                                       #
#         output/figures/landcover_map_bcn.png                                 #
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
library(stringr)
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

#3. Load monitoring datasets ####

df <- read.csv(
  file.path(input_dir,"data","s_index_18_24_ubms_cbms.csv"),
  sep=";"
)

cbms_sites <- read.csv(
  file.path(input_dir,"data","coord_length_transects_CBMS.csv"),
  sep=";"
)

ubms_visits <- read.csv(
  file.path(input_dir,"data","ubms_visits.csv")
)


#4. Basic preprocessing ####

df$SITE_ID <- as.factor(df$SITE_ID)
cbms_sites$SITE_ID <- as.factor(cbms_sites$SITE_ID)

cbms_sites$SCHEME <- "CBMS"

df_full <- left_join(df, cbms_sites, by=c("SITE_ID","SCHEME"))

df_full <- df_full %>%
  mutate(SINDEX = ifelse(SINDEX < 0, NA, SINDEX)) %>%
  dplyr::select(
    SCHEME,SITE_ID,SPECIES,YEAR,SINDEX,transect_length
  )


#5. Expand species–site–year combinations ####

site_years <- df_full %>% distinct(SITE_ID,YEAR)

species_sites <- df_full %>%
  group_by(SITE_ID,SPECIES) %>%
  summarise(ever_present = any(SINDEX>0,na.rm=TRUE),.groups="drop") %>%
  filter(ever_present) %>%
  dplyr::select(-ever_present)

expanded <- species_sites %>% left_join(site_years,by="SITE_ID")

df_expanded <- expanded %>%
  left_join(df_full,by=c("SITE_ID","SPECIES","YEAR")) %>%
  mutate(SINDEX = ifelse(is.na(SINDEX),0,SINDEX))


#6. Filter uBMS sites ####

df_ubms <- df_expanded %>% filter(SCHEME=="uBMS")

ubms_visits <- ubms_visits %>%
  filter(!grepl("_A$",transect_id)) %>%
  mutate(transect_id=gsub("_T","",transect_id))

visit_counts <- ubms_visits %>%
  group_by(transect_id,year) %>%
  summarise(n_visits=n(),.groups="drop") %>%
  rename(SITE_ID=transect_id)

good_sites <- visit_counts %>%
  filter(n_visits>=10) %>%
  group_by(SITE_ID) %>%
  summarise(years_10visits=n()) %>%
  filter(years_10visits>=5)

df_ubms_filtered <- df_ubms %>%
  filter(SITE_ID %in% good_sites$SITE_ID)


#7. Filter CBMS sites ####

df_cbms <- df_expanded_clean %>% filter(SCHEME=="CBMS")

df_cbms <- df_cbms %>%
  mutate(SINDEX_std = SINDEX * 1000 / transect_length)

balearic_sites <- c(60,61,101,145,146,155,170,171,173)

df_cbms <- df_cbms %>%
  filter(!SITE_ID %in% balearic_sites)

valid_sites <- c(
  8,21,26,29,33,34,40,58,68,69,75,76,80,83,88,89,95,106,107,108,
  114,118,120,121,128,129,134,135,136,147,148,149,150,152,157,158,
  159,161,174,175,178,179,184,191,192,196,213,223,226,230,233,234,
  237,238,239
)

df_cbms <- df_cbms %>% filter(SITE_ID %in% valid_sites)

df_cbms <- df_cbms %>%
  group_by(SITE_ID) %>%
  filter(n_distinct(YEAR)>=6) %>%
  ungroup()

df_cbms <- df_cbms %>%
  mutate(SITE_ID=paste0("ES-CTBMS.",SITE_ID))


#8. Load transect coordinates ####

m_coord_ubms <- read.csv(
  file.path(input_dir,"data","ubms_sites.csv"),
  sep=";"
)

m_coord_ubms_clean <- na.exclude(m_coord_ubms)

m_coord_cbms <- read.csv(
  file.path(input_dir,"data","coord_length_transects_CBMS.csv"),
  sep=";"
)

m_coord_cbms <- m_coord_cbms %>%
  mutate(transect_id=paste0("ES-CTBMS.",SITE_ID))

m_coord_ubms$transect_id <- substr(
  m_coord_ubms$transect_id,
  1,
  nchar(m_coord_ubms$transect_id)-2
)


#9. Convert coordinates to sf ####

ubms_sf <- st_as_sf(
  m_coord_ubms_clean,
  coords=c("transect_longitude","transect_latitude"),
  crs=4326
)

ubms_sf_3035 <- st_transform(ubms_sf,3035)

coords_3035 <- st_coordinates(ubms_sf_3035)

ubms_sf_3035$transect_lon <- coords_3035[,1]
ubms_sf_3035$transect_lat <- coords_3035[,2]


#10. Merge CBMS and uBMS transects ####

cbms_df <- m_coord_cbms %>%
  mutate(source="CBMS") %>%
  select(transect_id,transect_length,transect_lon,transect_lat,source)

ubms_df <- ubms_sf_3035 %>%
  st_drop_geometry() %>%
  mutate(source="uBMS") %>%
  select(transect_id,transect_length,transect_lon,transect_lat,source)

merged_transects <- bind_rows(cbms_df,ubms_df)

m_coord_sf <- st_as_sf(
  merged_transects,
  coords=c("transect_lon","transect_lat"),
  crs=3035
)

m_coord_sf_wgs <- st_transform(m_coord_sf,4326)


#11. Get Barcelona boundary ####

boundary_file <- file.path(
  input_dir,
  "city_boundaries",
  "barcelona_boundary.gpkg"
)

if (!file.exists(boundary_file)) {
  
  bcn_boundary <- opq("Barcelona") %>%
    add_osm_feature(key="name",value="Barcelona") %>%
    osmdata_sf() %>%
    .$osm_multipolygons %>%
    slice(1)
  
  st_write(barcelona_boundary,boundary_file,quiet=TRUE)
  
} else {
  
  barcelona_boundary <- st_read(boundary_file,quiet=TRUE)
}


#12. Create Barcelona buffer (60 km) ####

bcn_buffer <- barcelona_boundary %>%
  st_transform(25831) %>%
  st_buffer(60000) %>%
  st_transform(4326)


#13. Classify sites (inside / outside) ####

inside_flag <- st_within(
  m_coord_sf_wgs,
  barcelona_boundary,
  sparse=FALSE
)[,1]

m_coord_sf_wgs$context <- ifelse(inside_flag,"inside","outside")

m_coord_sf_wgs$context <- factor(
  m_coord_sf_wgs$context,
  levels=c("inside","outside")
)

m_coord_sf_wgs$in_buffer <- st_within(
  m_coord_sf_wgs,
  bcn_buffer,
  sparse=FALSE
)[,1]

sites_buffer <- m_coord_sf_wgs[m_coord_sf_wgs$in_buffer,]


#14. Load GHSL rasters ####

ghsl_dir <- file.path(input_dir,"ghsl")

ghsl_30_files <- c(
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_V1_0_R5_C19.tif"
)

ghsl_extract_30ss <- rast(
  file.path(ghsl_dir,
            "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_V1_0_R5_C19.tif")
)

ghsl_3_files <- c(
  "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R5_C19.tif"
)

ghsl_extract_3ss <- rast(
  file.path(ghsl_dir,
            "GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R5_C19.tif")
)


#15. Extract built-up surface ####

pts_etrs <- st_transform(sites_buffer,25831)

results_df <- sites_buffer %>%
  st_drop_geometry() %>%
  select(transect_id,context)

buffers <- c(1000,2000,5000)

for(r in buffers){
  
  buf_etrs <- st_buffer(pts_etrs,dist=r)
  buf_wgs  <- st_transform(buf_etrs,crs(ghsl_extract_30ss))
  
  vals <- terra::extract(
    ghsl_extract_30ss,
    vect(buf_wgs),
    fun=mean,
    na.rm=TRUE
  )[,2]
  
  results_df[[paste0("BUILT_UP_",r,"m")]] <- vals
}

built_up_df <- results_df %>%
  rename(
    SITE_ID = transect_id,
    CONTEXT = context,
    built1000 = BUILT_UP_1000m,
    built2000 = BUILT_UP_2000m,
    built5000 = BUILT_UP_5000m
  )

write.csv(
  built_up_df,
  file.path(output_dir,"data","built_up_bcn.csv"),
  row.names=FALSE
)


#16. Prepare raster for plotting ####

ghsl_crop <- crop(ghsl_extract_3ss,vect(bcn_buffer))
ghsl_mask <- mask(ghsl_crop,vect(bcn_buffer))

ghsl_df <- as.data.frame(ghsl_mask,xy=TRUE)

names(ghsl_df)[3] <- "built"
ghsl_df$built_log <- log1p(ghsl_df$built)


#17. Load land cover tiles ####

lc_files <- list.files(
  file.path(input_dir,"land_cover"),
  full.names=TRUE
)

lc_tiles <- lc_files[
  grepl("N39E000",lc_files)
]

lc_raster_bcn <- mosaic(
  sprc(lapply(lc_tiles,rast))
)

lc_crop_bcn <- crop(lc_raster_bcn,vect(bcn_buffer))
lc_mask_bcn <- mask(lc_crop_bcn,vect(bcn_buffer))

lc_raster_low_bcn <- aggregate(
  lc_mask_bcn,
  fact=25,
  fun=function(x){
    ux<-unique(x)
    ux[which.max(tabulate(match(x,ux)))]
  }
)

lc_df_bcn <- as.data.frame(lc_raster_low_bcn,xy=TRUE)

names(lc_df_bcn)[3] <- "class"

water_df_bcn <- lc_df_bcn %>% filter(class==80)


#18. Map: built-up surface ####

map_urb_bcn <- ggplot() +
  
  geom_raster(data=ghsl_df,aes(x=x,y=y,fill=built_log),alpha=0.5) +
  geom_sf(data=barcelona_boundary,fill=NA,color="black",size=1.5) +
  geom_sf(data=sites_buffer,aes(color=context),size=1) +
  geom_raster(data=water_df_bcn,aes(x=x,y=y),fill="white",alpha=1) +
  geom_sf(data=bcn_buffer,fill=NA,size=1.5) +
  
  scale_fill_viridis(option="D",name="Built-up surface") +
  scale_color_manual(values=c(inside="red",outside="black"),guide="none") +
  
  theme_minimal(base_family="Garamond",base_size=18) +
  theme(
    panel.grid=element_blank(),
    panel.border=element_rect(colour="black",fill=NA),
    axis.title=element_blank(),
    legend.position="none"
  ) +
  
  annotation_scale(location="bl")

ggsave(
  file.path(output_dir,"figures","urb_map_bcn.png"),
  map_urb_bcn,
  width=5.5,
  height=5,
  dpi=1200,
  bg="white"
)


#19. Land cover map ####

lc_colors <- c(
  "10"="#1a9850",
  "20"="#66c2a5",
  "30"="#cccc66",
  "40"="#ffcc33",
  "50"="#4D4D4D",
  "60"="#ffff99",
  "80"="white",
  "90"="#2c7bb6"
)

labels_lc <- c(
  "10"="Woodland",
  "20"="Shrubland",
  "30"="Grassland",
  "40"="Cropland",
  "50"="Built-up",
  "60"="Sparse vegetation",
  "80"="Water",
  "90"="Wetland"
)

map_lc_bcn <- ggplot() +
  
  geom_raster(data=lc_df_bcn,aes(x=x,y=y,fill=factor(class)),alpha=0.85) +
  geom_sf(data=bcn_buffer,fill=NA,color="black",size=1.2) +
  geom_sf(data=barcelona_boundary,fill="#e31a1c",alpha=0.4,color="#e31a1c",size=2) +
  
  geom_sf(data=subset(sites_buffer,context=="outside"),
          shape=21,fill="#0072B2",color="black",size=1.5) +
  
  geom_sf(data=subset(sites_buffer,context=="inside"),
          shape=21,fill="#D55E00",color="black",size=1.5) +
  
  scale_fill_manual(values=lc_colors,labels=labels_lc,name="Land cover") +
  
  coord_sf(crs=4326) +
  
  theme_minimal(base_family="Garamond",base_size=12) +
  theme(
    panel.grid=element_blank(),
    panel.border=element_rect(colour="black",fill=NA),
    axis.title=element_blank()
  ) +
  
  annotation_scale(location="bl")

ggsave(
  file.path(output_dir,"figures","landcover_map_bcn.png"),
  map_lc_bcn,
  width=5.5,
  height=5,
  dpi=1200,
  bg="white"
)


#20. Extract land cover diversity ####

pts_etrs <- st_transform(sites_buffer,crs(lc_raster_bcn))

calc_habdiv <- function(raster,geom,valid_classes=c(10,20,30,40,60)){
  
  cover_tab <- exact_extract(raster,geom,function(values,coverage_fraction){
    tapply(coverage_fraction,values,sum,na.rm=TRUE)
  })
  
  if(is.matrix(cover_tab)){
    props<-cover_tab[,1]
    names(props)<-rownames(cover_tab)
  }else{
    props<-cover_tab
  }
  
  names_num <- as.numeric(names(props))
  props <- props[names_num %in% valid_classes]
  props <- props/sum(props)
  
  if(length(props)==0 || all(is.na(props))) return(0)
  
  vegan::diversity(props,index="shannon")
}

for(r in buffers){
  
  buf_etrs <- st_buffer(pts_etrs,dist=r)
  buf_wgs  <- st_transform(buf_etrs,crs(lc_raster_bcn))
  
  habdiv_vals <- sapply(1:nrow(buf_wgs),function(i){
    calc_habdiv(lc_raster_bcn,buf_wgs[i,])
  })
  
  sites_buffer[[paste0("landdiv_",r,"m")]] <- habdiv_vals
}

landdiv_df <- sites_buffer %>%
  st_drop_geometry() %>%
  select(transect_id,context,starts_with("landdiv_"))

write.csv(
  landdiv_df,
  file.path(output_dir,"data","land_diversity_bcn.csv"),
  row.names=FALSE
)
