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
#         output/data/diversity_stability_lnd.csv                              #
#         output/data/diversity_stability_rnd.csv                              #
#                                                                              #
#   - Landscape variables                                                      #
#         output/data/built_up_bcn.csv                                         #
#         output/data/built_up_lnd.csv                                         #
#         output/data/built_up_rnd.csv                                         #
#                                                                              #
#         output/data/land_diversity_bcn.csv                                   #
#         output/data/land_diversity_lnd.csv                                   #
#         output/data/land_diversity_rnd.csv                                   #
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

input_dir  <- here("input", "data")
output_dir <- here("output")


#3. Load diversity–stability datasets ####

ds <- list(
  
  BCN = read.csv(
    file.path(output_dir, "data", "diversity_stability_bcn.csv")
  ),
  
  LND = read.csv(
    file.path(output_dir, "data", "diversity_stability_lnd.csv")
  ),
  
  RND = read.csv(
    file.path(output_dir, "data", "diversity_stability_rnd.csv")
  )
)

# ---- Fix SITE_ID format ----

clean_bcn_id <- function(x) {
  
  x <- as.character(x)
  
  # remove duplicated prefixes
  x <- gsub("^ES-CTBMS\\.", "", x)
  
  # keep uBMS sites unchanged
  x <- ifelse(
    grepl("^ES_uBMS", x),
    x,
    paste0("ES-CTBMS.", x)
  )
  
  x
}

ds$BCN <- ds$BCN %>%
  mutate(SITE_ID = clean_bcn_id(SITE_ID))

ds$LND <- ds$LND %>%
  mutate(SITE_ID = as.character(SITE_ID))

ds$RND <- ds$RND %>%
  mutate(SITE_ID = as.character(SITE_ID))


#4. Load and prepare landscape variables ####

prep_landscape <- function(built, landdiv) {
  
  built %>%
    distinct(SITE_ID, .keep_all = TRUE) %>%
    mutate(
      urban_context = as.numeric(CONTEXT == "inside")
    ) %>%
    select(-CONTEXT) %>%
    left_join(
      landdiv %>%
        distinct(SITE_ID, .keep_all = TRUE) %>%
        select(-CONTEXT),
      by = "SITE_ID"
    )
}

landscape <- list(
  
  BCN = prep_landscape(
    read.csv(file.path(output_dir,"data","built_up_bcn.csv")),
    read.csv(file.path(output_dir,"data","land_diversity_bcn.csv"))
  ),
  
  LND = prep_landscape(
    read.csv(file.path(output_dir,"data","built_up_lnd.csv")),
    read.csv(file.path(output_dir,"data","land_diversity_lnd.csv"))
  ),
  
  RND = prep_landscape(
    read.csv(file.path(output_dir,"data","built_up_rnd.csv")),
    read.csv(file.path(output_dir,"data","land_diversity_rnd.csv"))
  )
)


#5. Merge diversity–stability and landscape data ####

results <- list(
  
  BCN = left_join(ds$BCN, landscape$BCN, by = "SITE_ID"),
  
  LND = left_join(ds$LND, landscape$LND, by = "SITE_ID"),
  
  RND = left_join(ds$RND, landscape$RND, by = "SITE_ID")
)
results <- map(results, \(df) df %>% filter(!is.na(urban_context)))

#6. Standardize SEM variables ####

sem_vars <- c(
  "community_stability",
  "species_asynchrony",
  "wm_population_stability",
  "shannon_diversity",
  "FDis",
  "MPD",
  "species_richness",
  "built1000","built2000","built5000",
  "landdiv1000","landdiv2000","landdiv5000"
)

results <- map(
  results,
  \(df) df %>%
    mutate(
      across(
        where(is.numeric) & !matches("urban_context"),
        ~ as.numeric(scale(.))
      )
    )
)

#7. SEM specification ####

sem_template <- function(buffer) {
  
  paste0('
  community_stability ~ species_asynchrony + wm_population_stability

  species_asynchrony ~ shannon_diversity + FDis + MPD + species_richness
  wm_population_stability ~ shannon_diversity + FDis + MPD + species_richness

  shannon_diversity ~ species_richness + built', buffer, ' + landdiv', buffer, '
  MPD ~ species_richness + built', buffer, ' + landdiv', buffer, '
  FDis ~ species_richness + built', buffer, ' + landdiv', buffer, '

  species_richness ~ built', buffer, ' + landdiv', buffer, ' + urban_context

  built', buffer, ' ~ urban_context
  landdiv', buffer, ' ~ urban_context

  species_asynchrony ~~ wm_population_stability
  built', buffer, ' ~~ landdiv', buffer, '
  MPD ~~ FDis
  MPD ~~ shannon_diversity
  FDis ~~ shannon_diversity

  species_asynchrony ~~ urban_context
  wm_population_stability ~~ urban_context
  community_stability ~~ urban_context
  ')
}

buffers <- c(1000, 2000, 5000)


#8. Run SEMs ####

sem_results <- expand.grid(
  region = names(results),
  buffer = buffers
) %>%
  mutate(
    model = map2(region, buffer, function(r, b) {
      
      df <- results[[r]]
      
      built_var  <- paste0("built", b)
      landdiv_var <- paste0("landdiv", b)
      
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
        landdiv_var
      )
      
      if (!all(required_vars %in% names(df))) {
        message("Skipping SEM: ", r, " – buffer ", b)
        return(NULL)
      }
      
      sem(
        sem_template(b),
        data = df,
        missing = "fiml",
        em.h1.iter.max = 10000
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


#9. Extract SEM path coefficients ####

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
  file.path(output_dir,"results","SEM_paths_all_regions_buffers.csv"),
  row.names = FALSE
)

# ---- Note on SEM figure visualization ----
# The SEM path diagrams presented in the manuscript were generated from the
# standardized path coefficients obtained here. For clarity and graphical
# consistency, the final diagrams shown in the figures were subsequently
# assembled and edited using an external image editor.


#10. Prepare data for GLMs ####

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


#11. GLMs ####

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
    "results",
    "GLM_diversity_stability_context.csv"
  ),
  row.names = FALSE
)


#12. Plot GLMM effects #### 

effects <- map(
  biodiv_vars,
  \(v) ggpredict(
    glm_models[[v]],
    terms = c(v, "urban_context")
  )
)

names(effects) <- biodiv_vars

make_plot <- function(eff, xvar, data, xlab){
  
  ggplot() +
    
    geom_point(
      data = data,
      aes_string(
        x = xvar,
        y = "community_stability",
        color = "urban_context"
      ),
      alpha = 0.25,
      size = 1
    ) +
    
    geom_line(
      data = eff,
      aes(
        x = x,
        y = predicted,
        color = group
      ),
      linewidth = 1
    ) +
    
    geom_ribbon(
      data = eff,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = group
      ),
      alpha = 0.25
    ) +
    
    scale_color_manual(
      values = c("Urban" = "#D55E00", "Rural" = "#0072B2"),
      name = "Urban context"
    ) +
    
    scale_fill_manual(
      values = c("Urban" = "#D55E00", "Rural" = "#0072B2"),
      guide = "none"
    ) +
    
    labs(
      x = xlab,
      y = "Community stability"
    ) +
    
    theme_classic(base_family = "Garamond", base_size = 14) +
    
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      axis.line = element_blank(),
      legend.position = "bottom"
    )
}


p1 <- make_plot(
  effects$species_richness,
  "species_richness",
  all_cities,
  "Species richness"
)

p2 <- make_plot(
  effects$shannon_diversity,
  "shannon_diversity",
  all_cities,
  "Species diversity"
) +
  ylab(NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p3 <- make_plot(
  effects$FDis,
  "FDis",
  all_cities,
  "Functional diversity"
)

p4 <- make_plot(
  effects$MPD,
  "MPD",
  all_cities,
  "Phylogenetic diversity"
) +
  ylab(NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

final_plot <- ((p1 | p2) / (p3 | p4)) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank()
  )

final_plot


ggsave(
  filename = file.path(
    output_dir,
    "figures",
    "fig_diversity_stability_context.png"
  ),
  plot = final_plot,
  width = 5.5,
  height = 5.5,
  dpi = 600,
  bg = "white"
)


#13. Plot Diversity–stability relationships by region #### 

# ---- Helper function ----

make_city_plot <- function(data, xvar, xlab){
  
  ggplot(
    data,
    aes(
      x = .data[[xvar]],
      y = community_stability,
      color = urban_context,
      fill  = urban_context
    )
  ) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.2) +
    
    scale_color_manual(
      values = c("Urban" = "#D55E00", "Rural" = "#0072B2"),
      name = "Urban context"
    ) +
    
    scale_fill_manual(
      values = c("Urban" = "#D55E00", "Rural" = "#0072B2"),
      guide = "none"
    ) +
    
    labs(
      x = xlab,
      y = "Community stability"
    ) +
    
    facet_wrap(~city, scales = "free_x") +
    
    theme_classic(base_family = "Garamond", base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(family = "Garamond", size = 12, face = "bold")
    )
}

# ---- Generate plots ----

r1 <- make_city_plot(all_cities, "species_richness", "Species richness") +
  theme(strip.background = element_blank())

r2 <- make_city_plot(all_cities, "shannon_diversity", "Species diversity") +
  ylab(NULL) +
  theme(strip.background = element_blank())

r3 <- make_city_plot(all_cities, "FDis", "Functional diversity") +
  theme(strip.background = element_blank())

r4 <- make_city_plot(all_cities, "MPD", "Phylogenetic diversity") +
  ylab(NULL) +
  theme(strip.background = element_blank())

# ---- combined multipanel ----

alt_plot <- ((r1 | r2) / (r3 | r4)) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(family = "Garamond", size = 12),
    legend.margin = margin(t = -4, b = 0),
    strip.background = element_blank(),
    plot.tag = element_text(family = "Garamond", size = 16, face = "bold"),
    plot.margin = margin(t = 0, r = 2, b = 1, l = 2)
  )

alt_plot

# ---- save figure ----

ggsave(
  filename = here("output","figures","fig_diversity_stability_context_region.png"),
  plot = alt_plot,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)
