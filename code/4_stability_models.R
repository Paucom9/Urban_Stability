# ============================================================================ #
#   4_stability_models.R                                                       #
#   Author: Pau Colom                                                          #
#                                                                              #
#   This script analyses diversity–stability relationships across the          #
#   three study regions (Barcelona, London, Randstad).                         #
#   It fits structural equation models (SEMs) to evaluate the mechanisms       #
#   linking landscape structure, biodiversity, and community stability         #
#   across spatial buffers. It also fits GLMs to test interactions between     #
#   biodiversity metrics and urban context across cities.                      #
#                                                                              #
#   Input:                                                                     #
#   - Diversity and stability metrics                                          #
#         output/data/diversity_stability_bcn.csv                              #
#         output/data/diversity_stability_london.csv                           #
#         output/data/diversity_stability_randstad.csv                         #
#                                                                              #
#   - Landscape variables                                                      #
#         output/data/barcelona_built_up.csv                                   #
#         output/data/london_built_up.csv                                      #
#         output/data/randstad_built_up.csv                                    #
#                                                                              #
#         output/data/barcelona_land_diversity.csv                             #
#         output/data/london_land_diversity.csv                                #
#         output/data/randstad_land_diversity.csv                              #
#                                                                              #
#   Output:                                                                    #
#   - SEM path coefficients                                                    #
#         output/data/SEM_paths_all_regions_buffers.csv                        #
#                                                                              #
#   - GLM interaction results                                                  #
#         output/data/GLM_diversity_stability_interactions.csv                 #
# ============================================================================ #



#1. Libraries ####


library(dplyr)
library(tidyr)
library(purrr)
library(lavaan)
library(lavaanPlot)
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
library(patchwork)
library(broom)
library(here)


#2. Define project folders ####


input_dir  <- here("input")
output_dir <- here("output")


#3. Load diversity–stability datasets ####


ds <- list(
  
  BCN = read.csv(
    file.path(output_dir, "data", "diversity_stability_bcn.csv")
  ),
  
  LND = read.csv(
    file.path(output_dir, "data", "diversity_stability_london.csv")
  ),
  
  RND = read.csv(
    file.path(output_dir, "data", "diversity_stability_randstad.csv")
  )
)


#4. Load and prepare landscape variables ####


prep_landscape <- function(built, habdiv) {
  
  built %>%
    distinct(SITE_ID, .keep_all = TRUE) %>%
    mutate(
      urban_context = as.numeric(CONTEXT == "inside")
    ) %>%
    select(-CONTEXT) %>%
    left_join(
      habdiv %>%
        distinct(SITE_ID, .keep_all = TRUE) %>%
        select(-CONTEXT),
      by = "SITE_ID"
    )
}


landscape <- list(
  
  BCN = prep_landscape(
    read.csv(file.path(output_dir,"data","barcelona_built_up.csv")),
    read.csv(file.path(output_dir,"data","barcelona_land_diversity.csv"))
  ),
  
  LND = prep_landscape(
    read.csv(file.path(output_dir,"data","london_built_up.csv")),
    read.csv(file.path(output_dir,"data","london_land_diversity.csv"))
  ),
  
  RND = prep_landscape(
    read.csv(file.path(output_dir,"data","randstad_built_up.csv")),
    read.csv(file.path(output_dir,"data","randstad_land_diversity.csv"))
  )
)


#5. Merge diversity–stability and landscape data ####


results <- list(
  
  BCN = left_join(ds$BCN, landscape$BCN, by = "SITE_ID"),
  
  LND = left_join(ds$LND, landscape$LND, by = "SITE_ID"),
  
  RND = left_join(ds$RND, landscape$RND, by = "SITE_ID")
)


#6. SEM specification ####

sem_template <- function(buffer) {
  
  paste0('
  community_stability ~ species_asynchrony + wm_population_stability

  species_asynchrony ~ shannon_diversity + FDis + MPD + species_richness
  wm_population_stability ~ shannon_diversity + FDis + MPD + species_richness

  shannon_diversity ~ species_richness + built', buffer, ' + habdiv', buffer, '
  MPD ~ species_richness + built', buffer, ' + habdiv', buffer, '
  FDis ~ species_richness + built', buffer, ' + habdiv', buffer, '

  species_richness ~ built', buffer, ' + habdiv', buffer, ' + urban_context

  built', buffer, ' ~ urban_context
  habdiv', buffer, ' ~ urban_context

  species_asynchrony ~~ wm_population_stability
  built', buffer, ' ~~ habdiv', buffer, '
  MPD ~~ FDis
  MPD ~~ shannon_diversity
  FDis ~~ shannon_diversity

  species_asynchrony ~~ urban_context
  wm_population_stability ~~ urban_context
  community_stability ~~ urban_context
  ')
}


buffers <- c(1000, 2000, 5000)


#7. Run SEMs ####


sem_results <- expand.grid(
  region = names(results),
  buffer = buffers
) %>%
  mutate(
    model = map2(region, buffer, function(r, b) {
      
      df <- results[[r]]
      
      built_var  <- paste0("built", b)
      habdiv_var <- paste0("habdiv", b)
      
      required_vars <- c(
        "community_stability",
        "species_asynchrony",
        "wm_population_stability",
        "shannon_diversity",
        "FDis",
        "MPD",
        "species_richness",
        "urban_context",
        built_var,
        habdiv_var
      )
      
      if (!all(required_vars %in% names(df))) {
        message("Skipping SEM: ", r, " – buffer ", b)
        return(NULL)
      }
      
      sem(
        sem_template(b),
        data = df,
        missing = "fiml"
      )
      
    })
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  mutate(
    summary = map(
      model,
      ~ summary(.x,
                standardized = TRUE,
                fit.measures = TRUE)
    )
  )


#8. Extract SEM path coefficients ####


sem_paths <- sem_results %>%
  mutate(
    paths = map(
      model,
      ~ parameterEstimates(.x, standardized = TRUE) %>%
        filter(op == "~")
    )
  ) %>%
  select(region, buffer, paths) %>%
  unnest(paths)


write.csv(
  sem_paths,
  file.path(output_dir,"data","SEM_paths_all_regions_buffers.csv"),
  row.names = FALSE
)


#9. Prepare data for GLMs ####


all_cities <- bind_rows(
  mutate(results$BCN, city = "BCN"),
  mutate(results$LND, city = "LND"),
  mutate(results$RND, city = "RND")
)

all_cities$city <- factor(
  all_cities$city,
  levels = c("LND","RND","BCN")
)

all_cities$urban_context <- factor(
  all_cities$urban_context,
  levels = c(1,0),
  labels = c("Urban","Rural")
)


#10. GLMs ####

biodiv_vars <- c(
  "species_richness",
  "shannon_diversity",
  "FDis",
  "MPD"
)

glm_models <- map(
  biodiv_vars,
  ~ lm(
    as.formula(
      paste("community_stability ~", .x,
            "* urban_context + city")
    ),
    data = all_cities
  )
)

names(glm_models) <- biodiv_vars


glm_table <- map_df(
  names(glm_models),
  function(v) {
    
    broom::tidy(glm_models[[v]]) %>%
      mutate(variable = v)
    
  }
) %>%
  select(
    variable,
    term,
    estimate,
    std.error,
    statistic,
    p.value
  ) %>%
  mutate(
    estimate  = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2),
    p.value   = signif(p.value, 2)
  )


write.csv(
  glm_table,
  file.path(
    output_dir,
    "data",
    "GLM_diversity_stability_interactions.csv"
  ),
  row.names = FALSE
)

