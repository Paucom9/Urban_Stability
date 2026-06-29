# ============================================================================ #
#   3_diversity_stability_metrics.R                                            #
#   Author: Pau Colom                                                          #
#                                                                              #
#   This script calculates butterfly diversity and community stability         #
#   metrics for each monitoring site. Metrics include population stability,    #
#   community stability, species asynchrony, taxonomic diversity,              #
#   functional diversity (FDis), and phylogenetic diversity (MPD).             #
#                                                                              #
#   Input*:                                                                    #
#   - Annual abundance indices                                                 #
#         output/data/annual_abundance_bcn.csv                                 #
#         output/data/annual_abundance_lnd.csv                                 #
#         output/data/annual_abundance_rnd.csv                                 #
#                                                                              #
#   - Species traits                                                           #
#   (extracted from Middleton-Welling et al. 2020; Scientific Data, 7(1), 351) #  
#         input/data/trait_table.csv                                           #
#                                                                              #
#   - Butterfly phylogeny                                                      #
#   (extracted from Wiemers et al. 2020; ZooKeys 938, 97)                      #                                                     
#         input/data/EUROPEAN_BUTTERFLIES_FULLMCC_DROPTIPED.nwk                #
#                                                                              #
#   Output:                                                                    #
#   - Diversity and stability metrics                                          #
#         output/data/diversity_stability_bcn.csv                              #
#         output/data/diversity_stability_lnd.csv                              #
#         output/data/diversity_stability_rnd.csv                              #
#                                                                              #
#   *Note: Annual abundance indices are not publicly available                #
#   in this repository.                                                        #
#   Access requires a signed data-sharing agreement with the                   #
#   European Butterfly Monitoring Scheme (eBMS).                               #
#   Data requests can be submitted through:                                    #
#   https://butterfly-monitoring.net/                                          #
# ============================================================================ #



#1. Libraries ####

library(dplyr)
library(tidyr)
library(purrr)
library(vegan)
library(FD)
library(picante)
library(ape)
library(stringr)
library(tibble)
library(here)


#2. Define project folders ####

input_dir  <- here("input", "data")
output_dir <- here("output", "data")


#3. Load traits and phylogeny ####

traits_raw <- read.csv(
  file.path(input_dir, "species_trait_table.csv"),
  sep = ";"
)

traits_raw$OvS <- ifelse(
  traits_raw$OvS %in% c("E","L","P","A"),
  traits_raw$OvS,
  "Mixed"
)

traits_raw$Voltinism3 <- with(
  traits_raw,
  ifelse(
    Vol_max <= 1.5, "Univoltine",
    ifelse(Vol_max == 2, "Bivoltine",
           ifelse(Vol_max >= 3, "Multivoltine", NA))
  )
)

traits <- traits_raw %>%
  mutate(
    SPECIES = str_replace_all(Taxon, "_", " "),
    Voltinism3 = factor(
      Voltinism3,
      levels = c("Univoltine","Bivoltine","Multivoltine")
    )
  ) %>%
  select(SPECIES, Wi = WIn, OvS, Voltinism3, TrC)


tree <- read.tree(
  file.path(
    input_dir,
    "phylogeny.nwk"
  )
)

tree$tip.label <- gsub("_"," ",tree$tip.label)


#4. Diversity and stability calculation function ####

calculate_diversity_stability <- function(df, tree, traits, min_species = 2){
  
  # ---- Helper functions ----
  
  calc_evar <- function(x){
    x <- x[x > 0]
    if(length(x) < 2) return(NA_real_)
    v <- var(log(x))
    1 - (2/pi) * atan(v)
  }
  
  clean_species <- function(x) {
    x %>%
      as.character() %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_squish()
  }
  
  safe_neg_log <- function(x) {
    ifelse(is.na(x) | x <= 0, NA_real_, -log(x))
  }
  
  
  # ---- Standardise names and remove negative values ----
  
  df <- df %>%
    mutate(
      SPECIES = clean_species(SPECIES),
      SINDEX_std = ifelse(SINDEX_std < 0, NA, SINDEX_std)
    )
  
  traits <- traits %>%
    mutate(
      SPECIES = clean_species(SPECIES),
      Wi = as.numeric(Wi),
      OvS = as.factor(OvS),
      Voltinism3 = as.factor(Voltinism3),
      TrC = as.factor(TrC)
    )
  
  tree$tip.label <- clean_species(tree$tip.label)
  
  
  # ---- Collapse abundances per site-year-species ----
  
  df_sum <- df %>%
    group_by(SITE_ID, YEAR, SPECIES) %>%
    summarise(
      SINDEX_std = sum(SINDEX_std, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ========================================================================== #
  # Community stability decomposition
  # ========================================================================== #
  
  stability_components <- purrr::map_dfr(unique(df_sum$SITE_ID), function(site){
    
    site_dat <- df_sum %>%
      filter(SITE_ID == site)
    
    comm_wide <- site_dat %>%
      dplyr::select(YEAR, SPECIES, SINDEX_std) %>%
      pivot_wider(
        names_from = SPECIES,
        values_from = SINDEX_std,
        values_fill = 0
      ) %>%
      arrange(YEAR)
    
    mat <- comm_wide %>%
      dplyr::select(-YEAR) %>%
      as.matrix()
    
    storage.mode(mat) <- "numeric"
    
    # Keep species with positive mean abundance
    mu_i <- colMeans(mat, na.rm = TRUE)
    mat <- mat[, mu_i > 0, drop = FALSE]
    
    if(nrow(mat) < 3 || ncol(mat) < min_species) return(NULL)
    
    # Species-level statistics
    mu_i    <- colMeans(mat, na.rm = TRUE)
    sigma_i <- apply(mat, 2, sd, na.rm = TRUE)
    CV_i    <- sigma_i / mu_i
    
    valid_sp <- !is.na(mu_i) & !is.na(sigma_i) & mu_i > 0
    
    mu_i    <- mu_i[valid_sp]
    sigma_i <- sigma_i[valid_sp]
    CV_i    <- CV_i[valid_sp]
    mat     <- mat[, valid_sp, drop = FALSE]
    
    if(length(mu_i) < min_species) return(NULL)
    
    # Community-level statistics
    total_abundance <- rowSums(mat, na.rm = TRUE)
    
    mu_com    <- mean(total_abundance, na.rm = TRUE)
    sigma_com <- sd(total_abundance, na.rm = TRUE)
    CV_com    <- sigma_com / mu_com
    
    if(
      is.na(mu_com) || is.na(sigma_com) ||
      mu_com <= 0 || sigma_com <= 0 ||
      sum(sigma_i, na.rm = TRUE) <= 0
    ){
      return(NULL)
    }
    
    # Relative abundance weights
    w_i <- mu_i / sum(mu_i, na.rm = TRUE)
    
    # Weighted mean population variability
    CV_w <- sum(w_i * CV_i, na.rm = TRUE)
    
    # Loreau & de Mazancourt synchrony index
    phi <- (sigma_com^2) / (sum(sigma_i, na.rm = TRUE)^2)
    
    if(is.na(CV_w) || CV_w <= 0 || is.na(phi) || phi <= 0){
      return(NULL)
    }
    
    # Classical log-transformed decomposition
    community_stability     <- safe_neg_log(CV_com)
    wm_population_stability <- safe_neg_log(CV_w)
    species_asynchrony      <- safe_neg_log(phi)
    
    tibble(
      SITE_ID = site,
      
      mean_community_abundance = mu_com,
      sd_community_abundance   = sigma_com,
      n_years = nrow(mat),
      
      community_CV = CV_com,
      community_stability = community_stability,
      
      mean_population_abundance = mean(mu_i, na.rm = TRUE),
      mean_population_CV = mean(CV_i, na.rm = TRUE),
      mean_population_stability = safe_neg_log(mean(CV_i, na.rm = TRUE)),
      
      wm_population_CV = CV_w,
      wm_population_stability = wm_population_stability,
      
      synchrony_phi = phi,
      species_asynchrony = species_asynchrony,
      
      stability_identity_error =
        community_stability -
        (wm_population_stability + 0.5 * species_asynchrony)
    )
  })
  
  
  # ---- Species richness ----
  
  richness <- df_sum %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(
      total_abundance = sum(SINDEX_std, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(total_abundance > 0) %>%
    group_by(SITE_ID) %>%
    summarise(
      species_richness = n(),
      .groups = "drop"
    )
  
  
  # ---- Shannon diversity ----
  
  shannon <- df_sum %>%
    group_by(SITE_ID, YEAR) %>%
    summarise(
      abund = list(tapply(SINDEX_std, SPECIES, sum, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      shannon_diversity = map_dbl(
        abund,
        ~ if(sum(.x, na.rm = TRUE) > 0) {
          vegan::diversity(.x, index = "shannon")
        } else {
          NA_real_
        }
      ),
      evenness_sw = map_dbl(
        abund,
        ~ calc_evar(.x)
      )
    ) %>%
    group_by(SITE_ID) %>%
    summarise(
      shannon_diversity = mean(shannon_diversity, na.rm = TRUE),
      evenness_sw = mean(evenness_sw, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ---- Community matrix for FDis and MPD ----
  
  comm <- df_sum %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(
      abund = sum(SINDEX_std, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = SPECIES,
      values_from = abund,
      values_fill = 0
    ) %>%
    column_to_rownames("SITE_ID")
  
  
  # ---- Functional diversity (FDis) ----
  
  traits_mat <- traits %>%
    dplyr::select(SPECIES, Wi, OvS, Voltinism3, TrC) %>%
    distinct(SPECIES, .keep_all = TRUE) %>%
    filter(
      !is.na(Wi),
      !is.na(OvS),
      !is.na(Voltinism3),
      !is.na(TrC)
    ) %>%
    column_to_rownames("SPECIES")
  
  sp <- intersect(colnames(comm), rownames(traits_mat))
  
  comm_fdis <- comm[, sp, drop = FALSE]
  traits_fdis <- traits_mat[sp, , drop = FALSE]
  
  fdis <- FD::dbFD(
    traits_fdis,
    comm_fdis,
    calc.FRic = FALSE,
    calc.CWM = FALSE,
    messages = FALSE
  )$FDis
  
  fdis_df <- tibble(
    SITE_ID = rownames(comm_fdis),
    FDis = as.numeric(fdis)
  )
  
  
  # ---- Phylogenetic diversity (MPD) ----
  
  sp_phy <- intersect(tree$tip.label, colnames(comm))
  
  tree_p <- ape::drop.tip(
    tree,
    setdiff(tree$tip.label, sp_phy)
  )
  
  comm_p <- comm[, tree_p$tip.label, drop = FALSE]
  
  phylo_df <- tibble(
    SITE_ID = rownames(comm_p),
    MPD = picante::mpd(
      comm_p,
      cophenetic(tree_p),
      abundance.weighted = TRUE
    )
  )
  
  
  # ---- Merge all metrics ----
  
  results <- stability_components %>%
    left_join(richness, by = "SITE_ID") %>%
    left_join(shannon, by = "SITE_ID") %>%
    left_join(fdis_df, by = "SITE_ID") %>%
    left_join(phylo_df, by = "SITE_ID") %>%
    mutate(
      mean_community_abundance  = log1p(mean_community_abundance),
      mean_population_abundance = log1p(mean_population_abundance)
    )
  
  return(results)
}


#5. Define study regions ####

regions <- list(
  
  bcn = list(
    abund = file.path(output_dir,"annual_abundance_bcn.csv"),
    out   = file.path(output_dir,"diversity_stability_bcn.csv")
  ),
  
  london = list(
    abund = file.path(output_dir,"annual_abundance_lnd.csv"),
    out   = file.path(output_dir,"diversity_stability_lnd.csv")
  ),
  
  randstad = list(
    abund = file.path(output_dir,"annual_abundance_rnd.csv"),
    out   = file.path(output_dir,"diversity_stability_rnd.csv")
  )
  
)


#6. Run calculations ####

for (r in names(regions)) {
  
  message("=== Processing region: ", r, " ===")
  
  cfg <- regions[[r]]
  
  # ---- Load abundance ----
  
  df <- read.csv(cfg$abund, sep = ifelse(r == "bcn", ";", ",")) %>%
    rename(YEAR = any_of(c("YEAR","M_YEAR"))) %>%
    mutate(
      SINDEX_std = ifelse(SINDEX < 0, NA, SINDEX)
    )
  
  # ---- Expand SITE × SPECIES × YEAR matrix ----
  
  site_years <- df %>%
    distinct(SITE_ID, YEAR)
  
  species_sites <- df %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(
      ever_present = any(SINDEX_std > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(ever_present) %>%
    select(-ever_present)
  
  expanded <- species_sites %>%
    left_join(site_years, by = "SITE_ID", relationship = "many-to-many")
  
  df <- expanded %>%
    left_join(df, by = c("SITE_ID","SPECIES","YEAR")) %>%
    mutate(
      SINDEX_std = ifelse(is.na(SINDEX_std), 0, SINDEX_std)
    )
  
  # ---- Add traits ----
  
  df <- df %>%
    left_join(traits, by = "SPECIES")
  
  # ---- Calculate stability metrics ----
  
  res <- calculate_diversity_stability(df, tree, traits)
  
  # ---- Save output ----
  
  write.csv(res, cfg$out, row.names = FALSE)
}

