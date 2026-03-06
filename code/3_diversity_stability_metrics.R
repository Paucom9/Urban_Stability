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
#         output/data/bcn_annual_abundance.csv                                 #
#         output/data/lnd_annual_abundance.csv                                 #
#         output/data/rnd_annual_abundance.csv                                 #
#                                                                              #
#   - Urbanization variables                                                   #
#         output/data/barcelona_built_up.csv                                   #
#         output/data/london_built_up.csv                                      #
#         output/data/randstad_built_up.csv                                    #
#                                                                              #
#   - Landscape diversity datasets                                             #
#         output/data/barcelona_land_diversity.csv                             #
#         output/data/london_land_diversity.csv                                #
#         output/data/randstad_land_diversity.csv                              #
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
#   * Note: Annual abundance indices are not publicly available                #
#   in this repository.                                                        #
#   Access requires a signed data-sharing agreement with the                   #
#   European Butterfly Monitoring Scheme (eBMS).                               #
#   Data requests can be submitted through:                                    #
#   https://butterfly-monitoring.net/                                          #
# ============================================================================ #


# 1. Libraries ####

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

input_dir  <- here("input")
output_dir <- here("output")


#3. Load traits and phylogeny ####

traits_raw <- read.csv(
  file.path(input_dir, "data", "species_trait_table.csv"),
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
    "data",
    "phylogeny.nwk"
  )
)

tree$tip.label <- gsub("_"," ",tree$tip.label)




#4. Stability calculation function ####


calculate_stability <- function(df, tree, min_species = 2) {
  
  calc_evar <- function(x) {
    x <- x[x > 0]
    if (length(x) < 2) return(NA)
    v <- var(log(x))
    1 - (2 / pi) * atan(v)
  }
  
  df <- df %>%
    mutate(ABUND = ifelse(ABUND < 0, NA, ABUND))
  
  df_sum <- df %>%
    group_by(SITE_ID, YEAR, SPECIES) %>%
    summarise(
      ABUND = sum(ABUND, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Population stability ----
  
  pop_stats <- df_sum %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(
      mu = mean(ABUND, na.rm = TRUE),
      sd = sd(ABUND, na.rm = TRUE),
      n  = sum(!is.na(ABUND)),
      .groups = "drop"
    ) %>%
    filter(mu > 0, sd > 0, n >= 3) %>%
    mutate(stability = mu / sd)
  
  pop_site <- pop_stats %>%
    group_by(SITE_ID) %>%
    summarise(
      mean_population_abundance = mean(mu, na.rm = TRUE),
      mean_population_stability = mean(stability, na.rm = TRUE),
      wm_population_stability   = weighted.mean(stability, mu, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Community stability ----
  
  comm_year <- df_sum %>%
    group_by(SITE_ID, YEAR) %>%
    summarise(
      total = sum(ABUND, na.rm = TRUE),
      .groups = "drop"
    )
  
  comm_stats <- comm_year %>%
    group_by(SITE_ID) %>%
    summarise(
      mean_community_abundance = mean(total, na.rm = TRUE),
      sd_community_abundance   = sd(total, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(
      mean_community_abundance > 0,
      sd_community_abundance > 0,
      n >= 3
    ) %>%
    mutate(
      community_stability =
        mean_community_abundance /
        sd_community_abundance
    )
  
  # ---- Species asynchrony ----
  
  asynchrony <- map_dfr(
    unique(df_sum$SITE_ID),
    function(site) {
      
      mat <- df_sum %>%
        filter(SITE_ID == site) %>%
        pivot_wider(
          names_from = SPECIES,
          values_from = ABUND,
          values_fill = 0
        ) %>%
        arrange(YEAR) %>%
        select(where(is.numeric)) %>%
        as.matrix()
      
      mu <- colMeans(mat)
      mat <- mat[, mu > 0, drop = FALSE]
      
      if (ncol(mat) < min_species) return(NULL)
      
      phi <- sd(rowSums(mat))^2 /
        sum(apply(mat, 2, sd))^2
      
      tibble(
        SITE_ID = site,
        species_asynchrony = -log(phi)
      )
    }
  )
  
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
      shannon_diversity =
        map_dbl(
          data,
          ~ diversity(tapply(.x$ABUND, .x$SPECIES, sum))
        ),
      evenness_sw =
        map_dbl(
          data,
          ~ calc_evar(tapply(.x$ABUND, .x$SPECIES, sum))
        )
    ) %>%
    group_by(SITE_ID) %>%
    summarise(
      shannon_diversity = mean(shannon_diversity, na.rm = TRUE),
      evenness_sw = mean(evenness_sw, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Functional diversity ----
  
  comm <- df %>%
    group_by(SITE_ID, SPECIES) %>%
    summarise(
      abund = sum(ABUND, na.rm = TRUE),
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
  
  # ---- Phylogenetic diversity ----
  
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
  
  comm_stats %>%
    left_join(pop_site, by = "SITE_ID") %>%
    left_join(richness, by = "SITE_ID") %>%
    left_join(shannon, by = "SITE_ID") %>%
    left_join(asynchrony, by = "SITE_ID") %>%
    left_join(fdis_df, by = "SITE_ID") %>%
    left_join(phylo_df, by = "SITE_ID")
}



#5. Define study regions ####

regions <- list(
  
  bcn = list(
    sindex = file.path(output_dir,"data","sindex_bcn_18_23.csv"),
    built  = file.path(output_dir,"data","barcelona_built_up.csv"),
    habdiv = file.path(output_dir,"data","barcelona_land_diversity.csv"),
    out    = file.path(output_dir,"data","diversity_stability_bcn.csv")
  ),
  
  london = list(
    sindex = file.path(output_dir,"data","sindex_london_17_23.csv"),
    built  = file.path(output_dir,"data","london_built_up.csv"),
    habdiv = file.path(output_dir,"data","london_land_diversity.csv"),
    out    = file.path(output_dir,"data","diversity_stability_london.csv")
  ),
  
  randstad = list(
    sindex = file.path(output_dir,"data","sindex_randstad_17_23.csv"),
    built  = file.path(output_dir,"data","randstad_built_up.csv"),
    habdiv = file.path(output_dir,"data","randstad_land_diversity.csv"),
    out    = file.path(output_dir,"data","diversity_stability_randstad.csv")
  )
  
)


#6. Run calculations ####

for (r in names(regions)) {
  
  message("=== Processing region: ", r, " ===")
  
  cfg <- regions[[r]]
  
  df <- read.csv(cfg$sindex) %>%
    rename(YEAR = M_YEAR) %>%
    left_join(read.csv(cfg$built), by = "SITE_ID") %>%
    left_join(read.csv(cfg$habdiv), by = c("SITE_ID","CONTEXT")) %>%
    left_join(traits, by = "SPECIES")
  
  if (!"SINDEX_std" %in% names(df)) {
    df <- df %>% mutate(SINDEX_std = SINDEX)
  }
  
  df <- df %>%
    mutate(
      ABUND = SINDEX_std,
      urban_context = as.numeric(CONTEXT == "inside")
    )
  
  res <- calculate_stability(df, tree)
  
  res_std <- res %>%
    mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
  
  write.csv(res_std, cfg$out, row.names = FALSE)
}
