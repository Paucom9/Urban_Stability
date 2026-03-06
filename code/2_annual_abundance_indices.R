# ============================================================================ #
#   2_annual_abundance_indices.R                                               #
#   Author: Pau Colom                                                          #
#                                                                              #
#   This script calculates annual butterfly population indices (SINDEX)        #
#   from raw eBMS monitoring data using the rbms workflow.                     #
#   Indices are calculated separately for the three study systems:             #
#   Barcelona, London, and Randstad.                                           #
#                                                                              #
#   Input*:                                                                    #
#   - EBMS monitoring datasets:                                                #
#         input/data/ebms_count.csv                                            #
#         input/data/ebms_visit.csv                                            #
#                                                                              #
#   - Lists of valid transects per region (derived from land data):            #
#         output/data/bcn_built_up.csv                                         #
#         output/data/lnd_built_up.csv                                         #
#         output/data/rnd_built_up.csv                                         #
#                                                                              #
#   Output:                                                                    #
#   - Annual abundance indices (SINDEX)                                        #
#         output/data/annual_abundance_bcn.csv                                 #
#         output/data/annual_abundance_lnd.csv                                 #
#         output/data/annual_abundance_rnd.csv                                 #
#                                                                              #
#   *Note: Raw EBMS data are not publicly available in this repository.        #
#   Access requires a signed data-sharing agreement with the                   #
#   European Butterfly Monitoring Scheme (eBMS).                               #
#   Data requests can be submitted through:                                    #
#   https://butterfly-monitoring.net/                                          #
# ============================================================================ #



#1. Libraries ####

library(rbms)
library(dplyr)
library(lubridate)
library(here)


#2. Define project folders ####

input_dir  <- here("input")
output_dir <- here("output")


#3. Load EBMS monitoring data ####

m_count_raw <- read.csv(
  file.path(input_dir, "data", "ebms_count.csv")
)

m_visit_raw <- read.csv(
  file.path(input_dir, "data", "ebms_visit.csv")
)


#4. Harmonise monitoring datasets ####

m_count_all <- m_count_raw %>%
  transmute(
    SITE_ID  = as.character(transect_id),
    DATE     = as.character(visit_date),
    SPECIES  = as.character(species_name),
    COUNT    = as.integer(count),
    BMS_ID   = bms_id
  )

m_visit_all <- m_visit_raw %>%
  transmute(
    SITE_ID = as.character(transect_id),
    DATE    = as.character(visit_date),
    BMS_ID  = bms_id
  )


#5. Function to calculate annual indices (SINDEX) ####

calc_sindex <- function(m_count, m_visit, sites,
                        init_year, last_year = 2023,
                        min_fc = 0.10) {
  
  m_count <- m_count %>% filter(SITE_ID %in% sites)
  m_visit <- m_visit %>% filter(SITE_ID %in% sites)
  
  species_list <- sort(unique(m_count$SPECIES))
  
  ts_date <- rbms::ts_dwmy_table(
    InitYear = init_year,
    LastYear = last_year,
    WeekDay1 = "monday"
  )
  
  ts_season <- rbms::ts_monit_season(
    ts_date,
    StartMonth = 4, EndMonth = 9, StartDay = 1,
    CompltSeason = TRUE,
    Anchor = TRUE, AnchorLength = 2, AnchorLag = 2,
    TimeUnit = "d"
  )
  
  ts_visit <- rbms::ts_monit_site(ts_season, m_visit)
  
  out <- vector("list", length(species_list))
  names(out) <- species_list
  
  for (sp in species_list) {
    
    message("Processing: ", sp)
    
    ts_count <- rbms::ts_monit_count_site(ts_visit, m_count, sp)
    
    if (nrow(ts_count) == 0) next
    
    fc <- rbms::flight_curve(
      ts_count,
      NbrSample = 300,
      MinVisit = 10,
      MinOccur = 3,
      MinNbrSite = 1,
      MaxTrial = 4,
      GamFamily = "nb",
      SpeedGam = FALSE,
      CompltSeason = TRUE,
      TimeUnit = "d"
    )
    
    impt <- rbms::impute_count(ts_count, fc, TimeUnit = "d")
    
    sidx <- rbms::site_index(impt, MinFC = min_fc)
    sidx$SPECIES <- sp
    
    out[[sp]] <- sidx
  }
  
  sindex_all <- bind_rows(out)
  
  never_detected <- sindex_all %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(all_zero = all(SINDEX == 0), .groups = "drop") %>%
    filter(all_zero)
  
  anti_join(sindex_all, never_detected,
            by = c("SITE_ID", "SPECIES"))
}


#6. Define regions ####

regions <- list(
  
  barcelona = list(
    bms = c("ES-CTBMS", "ES-uBMS"),
    init_year = 2018,
    sites = read.csv(
      file.path(output_dir, "data", "bcn_built_up.csv")
    )$SITE_ID,
    out = file.path(
      output_dir, "data", "annual_abundance_bcn.csv"
    )
  ),
  
  london = list(
    bms = "UKBMS",
    init_year = 2017,
    sites = read.csv(
      file.path(output_dir, "data", "lnd_built_up.csv")
    )$SITE_ID,
    out = file.path(
      output_dir, "data", "annual_abundance_lnd.csv"
    )
  ),
  
  randstad = list(
    bms = "NLBMS",
    init_year = 2017,
    sites = read.csv(
      file.path(output_dir, "data", "rnd_built_up.csv")
    )$SITE_ID,
    out = file.path(
      output_dir, "data", "annual_abundance_rnd.csv"
    )
  )
)


#7. Run index calculation for each region ####

for (r in names(regions)) {
  
  cfg <- regions[[r]]
  
  message("=== Running region: ", r, " ===")
  
  m_count_r <- m_count_all %>%
    filter(BMS_ID %in% cfg$bms)
  
  m_visit_r <- m_visit_all %>%
    filter(BMS_ID %in% cfg$bms)
  
  res <- calc_sindex(
    m_count = m_count_r,
    m_visit = m_visit_r,
    sites = cfg$sites,
    init_year = cfg$init_year
  )
  
  write.csv(res, cfg$out, row.names = FALSE)
}
