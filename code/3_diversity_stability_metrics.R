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

calculate_diversity_stability <- function(df, tree, min_species = 2){
  
  # ---- Helper function (Smith & Wilson Evar) ----
  
  calc_evar <- function(x){
    x <- x[x > 0]
    if(length(x) < 2) return(NA)
    v <- var(log(x))
    1 - (2/pi) * atan(v)
  }
  
  
  # ---- Remove negative values ----
  
  df <- df %>%
    mutate(SINDEX_std = ifelse(SINDEX_std < 0, NA, SINDEX_std))
  
  
  # ---- Collapse abundances per site-year-species ----
  
  df_sum <- df %>%
    group_by(SITE_ID, YEAR, SPECIES) %>%
    summarise(
      SINDEX_std = sum(SINDEX_std, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ---- Population stability ----
  
  pop_stats <- df_sum %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(
      mean_abundance = mean(SINDEX_std, na.rm = TRUE),
      sd_abundance   = sd(SINDEX_std, na.rm = TRUE),
      n_years        = sum(!is.na(SINDEX_std)),
      .groups = "drop"
    ) %>%
    filter(mean_abundance > 0, sd_abundance > 0, n_years >= 3) %>%
    mutate(
      cv = sd_abundance / mean_abundance,
      stability = 1 / cv
    )
  
  pop_site <- pop_stats %>%
    group_by(SITE_ID) %>%
    summarise(
      mean_population_abundance = mean(mean_abundance, na.rm = TRUE),
      mean_population_stability = mean(stability, na.rm = TRUE),
      wm_population_stability   = weighted.mean(stability, mean_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ---- Community stability ----
  
  comm_year <- df_sum %>%
    group_by(SITE_ID, YEAR) %>%
    summarise(
      total_abundance = sum(SINDEX_std, na.rm = TRUE),
      .groups = "drop"
    )
  
  comm_stats <- comm_year %>%
    group_by(SITE_ID) %>%
    summarise(
      mean_community_abundance = mean(total_abundance, na.rm = TRUE),
      sd_community_abundance   = sd(total_abundance, na.rm = TRUE),
      n_years = n(),
      .groups = "drop"
    ) %>%
    filter(
      mean_community_abundance > 0,
      sd_community_abundance > 0,
      n_years >= 3
    ) %>%
    mutate(
      community_stability =
        1 / (sd_community_abundance / mean_community_abundance)
    )
  
  
  # ---- Species asynchrony ----
  
  site_ids <- unique(df_sum$SITE_ID)
  
  asynchrony <- map_dfr(site_ids, function(site){
    
    site_dat <- df_sum %>% filter(SITE_ID == site)
    
    comm <- site_dat %>%
      pivot_wider(
        names_from = SPECIES,
        values_from = SINDEX_std,
        values_fill = 0
      ) %>%
      arrange(YEAR)
    
    mat <- comm %>%
      select(where(is.numeric)) %>%
      as.matrix()
    
    mu <- colMeans(mat)
    mat <- mat[, mu > 0, drop = FALSE]
    
    if(ncol(mat) < min_species) return(NULL)
    
    sigma_i   <- apply(mat, 2, sd)
    sigma_com <- sd(rowSums(mat))
    
    phi <- (sigma_com^2) / (sum(sigma_i)^2)
    
    tibble(
      SITE_ID = site,
      species_asynchrony = -log(phi)
    )
  })
  
  
  # ---- Species richness ----
  
  richness <- df %>%
    group_by(SITE_ID) %>%
    summarise(
      species_richness = n_distinct(SPECIES),
      .groups = "drop"
    )
  
  
  # ---- Shannon diversity ----
  
  shannon <- df_sum %>%
    group_by(SITE_ID, YEAR) %>%
    nest() %>%
    mutate(
      shannon_diversity = map_dbl(
        data,
        ~ vegan::diversity(
          tapply(.x$SINDEX_std, .x$SPECIES, sum),
          index = "shannon")
      ),
      evenness_sw = map_dbl(
        data,
        ~ calc_evar(
          tapply(.x$SINDEX_std, .x$SPECIES, sum)
        )
      )
    ) %>%
    group_by(SITE_ID) %>%
    summarise(
      shannon_diversity = mean(shannon_diversity, na.rm = TRUE),
      evenness_sw = mean(evenness_sw, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # ---- Functional diversity (FDis) ----
  
  comm <- df %>%
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
  
  traits_mat <- traits %>%
    distinct(SPECIES, .keep_all = TRUE) %>%
    column_to_rownames("SPECIES")
  
  sp <- intersect(colnames(comm), rownames(traits_mat))
  
  fdis <- dbFD(
    traits_mat[sp, ],
    comm[, sp],
    calc.FRic = FALSE
  )$FDis
  
  fdis_df <- tibble(
    SITE_ID = rownames(comm),
    FDis = fdis
  )
  
  
  # ---- Phylogenetic diversity (MPD) ----
  
  sp_phy <- intersect(tree$tip.label, colnames(comm))
  
  tree_p <- drop.tip(
    tree,
    setdiff(tree$tip.label, sp_phy)
  )
  
  comm_p <- comm[, sp_phy]
  
  phylo_df <- tibble(
    SITE_ID = rownames(comm_p),
    MPD = mpd(
      comm_p,
      cophenetic(tree_p),
      abundance.weighted = TRUE
    )
  )
  
  
  # ---- Merge all metrics ----
  
  results <- comm_stats %>%
    left_join(pop_site, by = "SITE_ID") %>%
    left_join(richness, by = "SITE_ID") %>%
    left_join(shannon, by = "SITE_ID") %>%
    left_join(asynchrony, by = "SITE_ID") %>%
    left_join(fdis_df, by = "SITE_ID") %>%
    left_join(phylo_df, by = "SITE_ID") %>%
    mutate(
      community_stability        = log1p(community_stability),
      mean_community_abundance   = log1p(mean_community_abundance),
      mean_population_abundance  = log1p(mean_population_abundance),
      mean_population_stability  = log1p(mean_population_stability),
      wm_population_stability    = log1p(wm_population_stability)
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
  
  res <- calculate_diversity_stability(df, tree)
  
  # ---- Save output ----
  
  write.csv(res, cfg$out, row.names = FALSE)
}

