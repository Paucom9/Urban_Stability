# Clean environment
rm(list = ls())

# Charge packages
library(dplyr)
library(tidyr)
library(vegan)

# Load data
df<- read.csv("C:/Users/Usuario/Desktop/Urban Stability Project/Urban_Stability/Data/s_index_18_24_ubms_cbms.csv", sep = ";")
cbms_sites<- read.csv("C:/Users/Usuario/Desktop/Urban Stability Project/Urban_Stability/Data/coord_length_transects_CBMS.csv", sep = ";")


# Merge data

df$SITE_ID <- as.factor(df$SITE_ID)
cbms_sites$SITE_ID <- as.factor(cbms_sites$SITE_ID)
cbms_sites$SCHEME <- "CBMS"

df_m<- left_join(df, cbms_sites, by = c("SITE_ID", "SCHEME"))


#Standarize SINDEX of CBMS sites

df_m$SINDEX[df_m$SCHEME == "CBMS"] <- 
  (df_m$SINDEX[df_m$SCHEME == "CBMS"] / df_m$transect_length[df_m$SCHEME == "CBMS"]) * 1000


# Select variables
df_final <- df_m[, c("SCHEME", "SITE_ID", "SPECIES", "YEAR", "SINDEX")]

# Remove negative SINDEX values
df_final <- df_final %>% mutate(SINDEX = ifelse(SINDEX < 0, NA, SINDEX))


#### Calculate variables ####

calculate_variables <- function(df) {
  # Preprocess: set negative SINDEX values to NA
  df <- df %>% mutate(SINDEX = ifelse(SINDEX < 0, NA, SINDEX))
  
  ## 1. Population stability: per (SITE_ID, SPECIES)
  pop_stats <- df %>%
    group_by(SITE_ID, SPECIES) %>%
    summarize(
      mean_abundance = mean(SINDEX, na.rm = TRUE),
      sd_abundance   = sd(SINDEX, na.rm = TRUE),
      cv             = sd_abundance / mean_abundance,
      stability      = ifelse(cv == 0, NA, 1 / cv),
      .groups = "drop"
    )
  
  ## 2. Community stability: aggregate by SITE_ID and YEAR
  community_year <- df %>%
    group_by(SITE_ID, YEAR) %>%
    summarize(total_abundance = sum(SINDEX, na.rm = TRUE), .groups = "drop")
  
  comm_stats <- community_year %>%
    group_by(SITE_ID) %>%
    summarize(
      mean_community_abundance = mean(total_abundance, na.rm = TRUE),
      sd_community_abundance   = sd(total_abundance, na.rm = TRUE),
      cv_community_abundance   = sd_community_abundance / mean_community_abundance,
      stability  = ifelse(cv_community_abundance == 0, NA, 1 / cv_community_abundance),
      .groups = "drop"
    )
  
  ## 3. Species synchrony (community-level metric):
  # For each population (species and site) we already have sd_abundance
  pop_sd <- pop_stats %>% select(SITE_ID, SPECIES, sd_abundance)
  
  # Sum of species' SDs per site
  sum_pop_sd <- pop_sd %>%
    group_by(SITE_ID) %>%
    summarize(sum_pop_sd = sum(sd_abundance, na.rm = TRUE), .groups = "drop")
  
  # Variance of the community total abundance for each site
  comm_var <- community_year %>%
    group_by(SITE_ID) %>%
    summarize(var_comm = var(total_abundance, na.rm = TRUE), .groups = "drop")
  
  species_sync <- left_join(sum_pop_sd, comm_var, by = "SITE_ID") %>%
    mutate(species_synchrony = var_comm / (sum_pop_sd^2)) %>%
    select(SITE_ID, species_synchrony)
  
  ## 4. Other community-level metrics:
  ### (a) Species richness: count unique species per site
  species_richness <- df %>%
    group_by(SITE_ID) %>%
    summarize(species_richness = n_distinct(SPECIES), .groups = "drop")
  
  ### (b) Shannon diversity using vegan
  # Aggregate total abundance per SITE_ID and SPECIES
  species_abundance_site <- df %>%
    group_by(SITE_ID, SPECIES) %>%
    summarize(total_abundance = sum(SINDEX, na.rm = TRUE), .groups = "drop")
  
  # Pivot to create a community matrix with SITE_ID as rows and species as columns
  species_diversity <- species_abundance_site %>%
    pivot_wider(names_from = SPECIES, values_from = total_abundance, 
                values_fill = list(total_abundance = 0)) %>%
    mutate(across(-SITE_ID, ~ifelse(. < 0, 0, .))) %>%  # ensure non-negative
    mutate(shannon_diversity = diversity(as.matrix(select(., -SITE_ID)), index = "shannon")) %>%
    select(SITE_ID, shannon_diversity)
  
  ### (c) Mean population abundance: average of species' mean abundances per site
  mean_pop_abundance <- pop_stats %>%
    group_by(SITE_ID) %>%
    summarize(mean_population_abundance = mean(mean_abundance, na.rm = TRUE), .groups = "drop")
  
  ### (d) Mean population stability: average of species' stability values per site
  mean_pop_stability <- pop_stats %>%
    group_by(SITE_ID) %>%
    summarize(mean_population_stability = mean(stability, na.rm = TRUE), .groups = "drop")
  
  # Combine all "other metrics" into one data frame
  other_metrics <- species_richness %>%
    left_join(species_diversity, by = "SITE_ID") %>%
    left_join(mean_pop_abundance, by = "SITE_ID") %>%
    left_join(mean_pop_stability, by = "SITE_ID") %>%
    left_join(species_sync, by = "SITE_ID")
  
  # Now combine community stability (comm_stats) with the other metrics
  combined_comm_stats <- left_join(comm_stats, other_metrics, by = "SITE_ID")
  
  return(list(
    population_stats  = pop_stats,
    community_stats   = combined_comm_stats
  ))
}

# Calculate stability metrics
results <- calculate_variables(df_final)

# Inspect the outputs
head(results$population_stats)   # Population stats (species*site)
head(results$community_stats)    # Combined community-level metrics per SITE_ID
