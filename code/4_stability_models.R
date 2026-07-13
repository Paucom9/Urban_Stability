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
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
library(patchwork)
library(broom)
library(here)
library(stringr)
library(readr)
library(sf)
library(adespatial)
library(spdep)
library(tibble)


#2. Define project folders ####

input_dir  <- here("input", "data")
output_dir <- here("output")


#3. Helper functions ####

clean_bcn_id <- function(x) {
  
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  # Keep missing values as missing
  x[is.na(x) | x == ""] <- NA_character_
  
  # Remove duplicated CBMS prefixes
  x <- gsub("^ES-CTBMS\\.", "", x)
  
  # Keep uBMS IDs unchanged; add CBMS prefix only to CBMS numeric IDs
  x <- dplyr::case_when(
    is.na(x) ~ NA_character_,
    grepl("^ES_uBMS", x) ~ x,
    TRUE ~ paste0("ES-CTBMS.", x)
  )
  
  x
}


clean_site_id <- function(x, region) {
  
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  if (region == "BCN") {
    x <- clean_bcn_id(x)
  }
  
  x
}


prep_landscape <- function(built, landdiv, region) {
  
  built_clean <- built %>%
    mutate(
      SITE_ID = clean_site_id(SITE_ID, region),
      CONTEXT = stringr::str_trim(as.character(CONTEXT))
    ) %>%
    filter(!is.na(SITE_ID)) %>%
    distinct(SITE_ID, .keep_all = TRUE) %>%
    rename(CONTEXT_built = CONTEXT)
  
  land_clean <- landdiv %>%
    mutate(
      SITE_ID = clean_site_id(SITE_ID, region),
      CONTEXT = stringr::str_trim(as.character(CONTEXT))
    ) %>%
    filter(!is.na(SITE_ID)) %>%
    distinct(SITE_ID, .keep_all = TRUE) %>%
    rename(CONTEXT_land = CONTEXT)
  
  # Important:
  # Use full_join so sites present in landdiv but missing from built
  # are not lost. This prevents uBMS sites from being silently dropped.
  full_join(
    built_clean,
    land_clean,
    by = "SITE_ID"
  ) %>%
    mutate(
      CONTEXT = dplyr::coalesce(CONTEXT_built, CONTEXT_land),
      urban_context = as.numeric(CONTEXT == "inside")
    ) %>%
    select(
      -CONTEXT_built,
      -CONTEXT_land,
      -CONTEXT
    )
}


join_region_data <- function(ds_region, landscape_region, region) {
  
  out <- ds_region %>%
    left_join(landscape_region, by = "SITE_ID")
  
  out %>%
    filter(!is.na(urban_context))
}


#4. Load diversity–stability datasets ####

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


#5. Clean SITE_IDs ####

ds$BCN <- ds$BCN %>%
  mutate(SITE_ID = clean_site_id(SITE_ID, "BCN"))

ds$LND <- ds$LND %>%
  mutate(SITE_ID = clean_site_id(SITE_ID, "LND"))

ds$RND <- ds$RND %>%
  mutate(SITE_ID = clean_site_id(SITE_ID, "RND"))


#6. Load and prepare landscape variables ####

landscape <- list(
  
  BCN = prep_landscape(
    read.csv(file.path(output_dir, "data", "built_up_bcn.csv")),
    read.csv(file.path(output_dir, "data", "land_diversity_bcn.csv")),
    region = "BCN"
  ),
  
  LND = prep_landscape(
    read.csv(file.path(output_dir, "data", "built_up_lnd.csv")),
    read.csv(file.path(output_dir, "data", "land_diversity_lnd.csv")),
    region = "LND"
  ),
  
  RND = prep_landscape(
    read.csv(file.path(output_dir, "data", "built_up_rnd.csv")),
    read.csv(file.path(output_dir, "data", "land_diversity_rnd.csv")),
    region = "RND"
  )
)



#7. Merge diversity–stability and landscape data ####

results <- list(
  
  BCN = join_region_data(ds$BCN, landscape$BCN, "BCN"),
  
  LND = join_region_data(ds$LND, landscape$LND, "LND"),
  
  RND = join_region_data(ds$RND, landscape$RND, "RND")
)


#8. Spatially corrected GLMs using dbMEM ####

# -------------------------------------------------------------------------
# This section fits the 12 GLMs requested by the reviewer:
#   3 responses:
#     - community_stability
#     - species_asynchrony
#     - wm_population_stability
#
#   4 biodiversity predictors:
#     - shannon_diversity
#     - species_richness
#     - MPD
#     - FDis
#
# Model structure:
#   response ~ diversity * urban_context + city + dbMEM filters
#
# Spatial filtering:
#   dbMEMs are calculated separately by city.
#   Moran's I is then assessed on model residuals separately within each
#   city × urban_context group.
#   For each model, the minimum number of city-level dbMEMs is added until
#   no significant positive residual spatial autocorrelation remains.
# -------------------------------------------------------------------------


#8.1 Settings ####

response_vars <- c(
  "community_stability",
  "species_asynchrony",
  "wm_population_stability"
)

biodiv_vars <- c(
  "shannon_diversity",
  "species_richness",
  "MPD",
  "FDis"
)

set.seed(123)

max_dbmem <- 30
alpha_moran <- 0.05
moran_nsim <- 9999
min_group_sites <- 10

dir.create(
  file.path(output_dir, "results"),
  recursive = TRUE,
  showWarnings = FALSE
)


#8.2 Prepare GLM dataset before adding coordinates ####

all_cities_no_coords <- all_cities_check %>%
  mutate(
    SITE_ID = stringr::str_trim(as.character(SITE_ID)),
    city = factor(city, levels = c("LND", "RND", "BCN")),
    urban_context = factor(urban_context, levels = c("Urban", "Rural"))
  )


#8.3 Load coordinates from eBMS and uBMS ####

# ---- eBMS coordinates: UK, NL and BCN-CBMS ----
# transect_lon / transect_lat are already projected coordinates.

ebms_coords_raw <- read.csv(
  file.path(input_dir, "ebms_transect_coord.csv")
)

ebms_coords <- ebms_coords_raw %>%
  transmute(
    SITE_ID = stringr::str_trim(as.character(transect_id)),
    bms_id = stringr::str_trim(as.character(bms_id)),
    x_3035 = as.numeric(transect_lon),
    y_3035 = as.numeric(transect_lat)
  ) %>%
  filter(!is.na(SITE_ID), !is.na(x_3035), !is.na(y_3035)) %>%
  mutate(
    city = case_when(
      bms_id == "UKBMS" ~ "LND",
      bms_id == "NLBMS" ~ "RND",
      bms_id == "ES-CTBMS" ~ "BCN",
      TRUE ~ NA_character_
    ),
    city = factor(city, levels = c("LND", "RND", "BCN"))
  ) %>%
  filter(!is.na(city)) %>%
  distinct(city, SITE_ID, .keep_all = TRUE)


# ---- uBMS coordinates: same cleaning as original BCN script ----
# Remove _A sections and remove _T suffix to match SITE_ID in the model.

ubms_raw <- read.csv(
  file.path(input_dir, "ubms_sites.csv"),
  sep = ";",
  header = TRUE,
  stringsAsFactors = FALSE
)

ubms_coords <- ubms_raw %>%
  filter(
    !is.na(transect_longitude),
    !is.na(transect_latitude)
  ) %>%
  filter(
    !stringr::str_detect(transect_id, "_A$")
  ) %>%
  mutate(
    SITE_ID = stringr::str_remove(transect_id, "_T$"),
    SITE_ID = stringr::str_trim(as.character(SITE_ID)),
    longitude = as.numeric(transect_longitude),
    latitude  = as.numeric(transect_latitude)
  ) %>%
  filter(
    !is.na(SITE_ID),
    SITE_ID != "",
    !is.na(longitude),
    !is.na(latitude)
  )

ubms_sf <- sf::st_as_sf(
  ubms_coords,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

ubms_xy <- sf::st_coordinates(
  sf::st_transform(ubms_sf, 3035)
)

ubms_coords <- ubms_coords %>%
  mutate(
    city = factor("BCN", levels = c("LND", "RND", "BCN")),
    x_3035 = ubms_xy[, 1],
    y_3035 = ubms_xy[, 2]
  ) %>%
  select(city, SITE_ID, x_3035, y_3035) %>%
  distinct(city, SITE_ID, .keep_all = TRUE)


# ---- Combine coordinate sources ----

site_coords <- bind_rows(
  ebms_coords %>% select(city, SITE_ID, x_3035, y_3035),
  ubms_coords
) %>%
  distinct(city, SITE_ID, .keep_all = TRUE)


#8.4 Check coordinate merge ####

missing_coord_sites <- all_cities_no_coords %>%
  distinct(city, SITE_ID, urban_context) %>%
  left_join(site_coords, by = c("city", "SITE_ID")) %>%
  filter(is.na(x_3035) | is.na(y_3035)) %>%
  arrange(city, SITE_ID)

message("\n==== Sites still missing coordinates ====")
message("n = ", nrow(missing_coord_sites))
print(tibble::as_tibble(missing_coord_sites), n = 200)


all_cities <- all_cities_no_coords %>%
  select(-any_of(c("longitude", "latitude", "x_3035", "y_3035", "x_dbmem", "y_dbmem"))) %>%
  left_join(site_coords, by = c("city", "SITE_ID"))


# Standardize response and biodiversity variables within city
vars_to_standardize <- c(
  "community_stability",
  "species_asynchrony",
  "wm_population_stability",
  "species_richness",
  "shannon_diversity",
  "FDis",
  "MPD"
)

all_cities <- all_cities %>%
  group_by(city) %>%
  mutate(
    across(
      all_of(vars_to_standardize),
      ~ as.numeric(scale(.))
    )
  ) %>%
  ungroup()


message("\n==== Final GLM/dbMEM dataset check ====")

all_cities %>%
  mutate(missing_coords = is.na(x_3035) | is.na(y_3035)) %>%
  count(city, urban_context, missing_coords) %>%
  tibble::as_tibble() %>%
  print(n = 100)


message("\n==== Standardization check within city ====")

all_cities %>%
  group_by(city) %>%
  summarise(
    across(
      all_of(vars_to_standardize),
      list(
        mean = ~ round(mean(.x, na.rm = TRUE), 6),
        sd = ~ round(sd(.x, na.rm = TRUE), 6)
      )
    ),
    .groups = "drop"
  ) %>%
  tibble::as_tibble() %>%
  print(n = 100)

#8.5 Helper: jitter duplicated coordinates ####

jitter_duplicate_coords <- function(
    df,
    group_cols = c("city"),
    xcol = "x_3035",
    ycol = "y_3035",
    amount = 1
) {
  
  df %>%
    group_by(across(all_of(c(group_cols, xcol, ycol)))) %>%
    mutate(
      .dup_id = row_number(),
      .dup_n = n(),
      x_dbmem = .data[[xcol]] + if_else(.dup_n > 1, (.dup_id - 1) * amount, 0),
      y_dbmem = .data[[ycol]] + if_else(.dup_n > 1, (.dup_id - 1) * amount, 0)
    ) %>%
    ungroup() %>%
    select(-.dup_id, -.dup_n)
}


#8.6 Helper: compute dbMEMs by city ####

compute_group_dbmem <- function(
    df,
    group_cols = c("city"),
    coords = c("x_dbmem", "y_dbmem"),
    nvec = 1,
    min_group_sites = 10
) {
  
  if (nvec == 0) {
    return(matrix(numeric(0), nrow = nrow(df), ncol = 0))
  }
  
  n <- nrow(df)
  
  group_keys <- df %>%
    group_by(across(all_of(group_cols))) %>%
    group_keys()
  
  group_data <- df %>%
    mutate(.row_id = row_number()) %>%
    group_by(across(all_of(group_cols))) %>%
    group_split()
  
  ev_list <- list()
  
  for (i in seq_along(group_data)) {
    
    sub <- group_data[[i]]
    idx <- sub$.row_id
    
    if (nrow(sub) < min_group_sites) next
    
    xy <- as.matrix(sub[, coords])
    
    if (nrow(unique(as.data.frame(xy))) < min_group_sites) next
    
    mem <- tryCatch(
      adespatial::dbmem(
        xy,
        MEM.autocor = "positive",
        silent = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(mem)) next
    
    mem_mat <- as.matrix(mem)
    
    if (ncol(mem_mat) == 0) next
    
    take <- min(nvec, ncol(mem_mat))
    mem_mat <- mem_mat[, seq_len(take), drop = FALSE]
    
    group_label <- paste(as.character(group_keys[i, ]), collapse = "_")
    group_label <- make.names(group_label)
    
    padded <- matrix(0, nrow = n, ncol = ncol(mem_mat))
    padded[idx, ] <- mem_mat
    
    colnames(padded) <- paste0(
      "dbMEM_",
      group_label,
      "_",
      seq_len(ncol(mem_mat))
    )
    
    ev_list[[length(ev_list) + 1]] <- padded
  }
  
  if (length(ev_list) == 0) {
    return(matrix(numeric(0), nrow = n, ncol = 0))
  }
  
  do.call(cbind, ev_list)
}


#8.7 Helper: Moran's I by city × urban_context ####

test_moran_groups <- function(
    df,
    residuals_vec,
    group_cols = c("city", "urban_context"),
    coords = c("x_dbmem", "y_dbmem"),
    nsim = 9999,
    k = 10,
    min_group_sites = 10
) {
  
  df$.resid <- residuals_vec
  
  group_keys <- df %>%
    group_by(across(all_of(group_cols))) %>%
    group_keys()
  
  group_data <- df %>%
    group_by(across(all_of(group_cols))) %>%
    group_split()
  
  out <- vector("list", length(group_data))
  
  for (i in seq_along(group_data)) {
    
    sub <- group_data[[i]]
    key <- group_keys[i, , drop = FALSE]
    
    if (nrow(sub) < min_group_sites) {
      
      out[[i]] <- bind_cols(
        key,
        tibble(
          n_sites = nrow(sub),
          moran_I = NA_real_,
          moran_p = NA_real_,
          moran_error = "Too few sites"
        )
      )
      
      next
    }
    
    xy <- as.matrix(sub[, coords])
    
    if (nrow(unique(as.data.frame(xy))) < min_group_sites) {
      
      out[[i]] <- bind_cols(
        key,
        tibble(
          n_sites = nrow(sub),
          moran_I = NA_real_,
          moran_p = NA_real_,
          moran_error = "Too few unique coordinates"
        )
      )
      
      next
    }
    
    moran_res <- tryCatch({
      
      k_use <- min(
        k,
        floor((nrow(sub) - 1) / 3),
        nrow(sub) - 1
      )
      
      k_use <- max(k_use, 1)
      
      nb <- spdep::knn2nb(
        spdep::knearneigh(
          xy,
          k = k_use
        )
      )
      
      lw <- spdep::nb2listw(
        nb,
        style = "W",
        zero.policy = TRUE
      )
      
      mc <- spdep::moran.mc(
        sub$.resid,
        lw,
        nsim = nsim,
        zero.policy = TRUE
      )
      
      tibble(
        n_sites = nrow(sub),
        moran_I = as.numeric(mc$statistic),
        moran_p = mc$p.value,
        moran_error = NA_character_
      )
      
    }, error = function(e) {
      
      tibble(
        n_sites = nrow(sub),
        moran_I = NA_real_,
        moran_p = NA_real_,
        moran_error = e$message
      )
    })
    
    out[[i]] <- bind_cols(key, moran_res)
  }
  
  bind_rows(out)
}


#8.8 Helper: fit one dbMEM-corrected model ####

fit_one_dbmem_model <- function(
    data,
    response,
    diversity,
    max_dbmem = 30,
    alpha_moran = 0.05,
    moran_nsim = 9999,
    min_group_sites = 10
) {
  
  model_vars <- c(
    "SITE_ID",
    "city",
    "urban_context",
    response,
    diversity,
    "x_3035",
    "y_3035"
  )
  
  df_model <- data %>%
    select(all_of(model_vars)) %>%
    filter(
      complete.cases(
        across(all_of(c(response, diversity, "city", "urban_context", "x_3035", "y_3035")))
      )
    ) %>%
    droplevels() %>%
    jitter_duplicate_coords(
      group_cols = c("city")
    )
  
  base_formula <- as.formula(
    paste(response, "~", diversity, "* urban_context + city")
  )
  
  final_fit <- NULL
  final_data <- NULL
  final_formula <- NULL
  final_moran <- NULL
  selected_n <- max_dbmem
  
  for (nv in 0:max_dbmem) {
    
    message(
      "\nTesting model: ",
      response,
      " ~ ",
      diversity,
      " * urban_context + city + city-level dbMEMs; n dbMEM = ",
      nv
    )
    
    if (nv == 0) {
      
      model_data <- df_model
      model_formula <- base_formula
      
    } else {
      
      mem_mat <- compute_group_dbmem(
        df = df_model,
        group_cols = c("city"),
        coords = c("x_dbmem", "y_dbmem"),
        nvec = nv,
        min_group_sites = min_group_sites
      )
      
      if (ncol(mem_mat) > 0) {
        
        mem_df <- as.data.frame(mem_mat)
        model_data <- bind_cols(df_model, mem_df)
        
        mem_terms <- names(mem_df)
        
        model_formula <- as.formula(
          paste(
            response,
            "~",
            diversity,
            "* urban_context + city +",
            paste(mem_terms, collapse = " + ")
          )
        )
        
      } else {
        
        model_data <- df_model
        model_formula <- base_formula
      }
    }
    
    fit <- lm(model_formula, data = model_data)
    
    moran_tbl <- test_moran_groups(
      df = model_data,
      residuals_vec = residuals(fit),
      group_cols = c("city", "urban_context"),
      coords = c("x_dbmem", "y_dbmem"),
      nsim = moran_nsim,
      k = 10,
      min_group_sites = min_group_sites
    ) %>%
      mutate(n_dbmem = nv)
    
    print(tibble::as_tibble(moran_tbl), n = 100)
    
    if (all(is.na(moran_tbl$moran_p))) {
      stop(
        "Moran's I could not be calculated for any group in model: ",
        response, " ~ ", diversity,
        ". Check coordinates and spatial weights."
      )
    }
    
    remaining_autocorrelation <- moran_tbl %>%
      filter(
        is.na(moran_error),
        !is.na(moran_p),
        moran_I > 0,
        moran_p <= alpha_moran
      ) %>%
      nrow()
    
    final_fit <- fit
    final_data <- model_data
    final_formula <- model_formula
    final_moran <- moran_tbl
    selected_n <- nv
    
    if (remaining_autocorrelation == 0) {
      message("Selected ", nv, " city-level dbMEM vectors.")
      break
    }
  }
  
  final_significant_moran <- final_moran %>%
    filter(
      is.na(moran_error),
      !is.na(moran_p),
      moran_I > 0,
      moran_p <= alpha_moran
    )
  
  moran_corrected <- nrow(final_significant_moran) == 0
  
  if (!moran_corrected) {
    message(
      "WARNING: residual positive spatial autocorrelation remains for model: ",
      response, " ~ ", diversity,
      " after ", selected_n, " city-level dbMEMs."
    )
  }
  
  positive_moran <- final_moran$moran_I[final_moran$moran_I > 0]
  
  coef_tbl <- broom::tidy(final_fit) %>%
    mutate(
      response = response,
      diversity = diversity,
      n_dbmem = selected_n,
      dbmem_grouping = "city",
      moran_grouping = "city_x_urban_context",
      formula = paste(deparse(final_formula), collapse = " "),
      is_spatial_filter = stringr::str_detect(term, "^dbMEM_")
    ) %>%
    select(
      response,
      diversity,
      n_dbmem,
      dbmem_grouping,
      moran_grouping,
      term,
      estimate,
      std.error,
      statistic,
      p.value,
      is_spatial_filter,
      formula
    )
  
  summary_tbl <- tibble(
    response = response,
    diversity = diversity,
    n_dbmem = selected_n,
    n_obs = nobs(final_fit),
    r.squared = summary(final_fit)$r.squared,
    adj.r.squared = summary(final_fit)$adj.r.squared,
    dbmem_grouping = "city",
    moran_grouping = "city_x_urban_context",
    moran_corrected = moran_corrected,
    n_significant_positive_moran = nrow(final_significant_moran),
    significant_moran_groups = paste(
      final_significant_moran$city,
      final_significant_moran$urban_context,
      collapse = "; "
    ),
    max_positive_moran_I = ifelse(
      length(positive_moran) == 0,
      NA_real_,
      max(positive_moran, na.rm = TRUE)
    ),
    min_moran_p = suppressWarnings(
      min(final_moran$moran_p, na.rm = TRUE)
    )
  )
  
  moran_tbl <- final_moran %>%
    mutate(
      response = response,
      diversity = diversity,
      dbmem_grouping = "city",
      moran_grouping = "city_x_urban_context"
    ) %>%
    select(
      response,
      diversity,
      n_dbmem,
      dbmem_grouping,
      moran_grouping,
      city,
      urban_context,
      n_sites,
      moran_I,
      moran_p,
      moran_error
    )
  
  list(
    model = final_fit,
    coefficients = coef_tbl,
    summary = summary_tbl,
    moran = moran_tbl,
    data = final_data
  )
}


#8.9 Run the 12 models ####

glm_dbmem_fits <- tidyr::expand_grid(
  response = response_vars,
  diversity = biodiv_vars
) %>%
  mutate(
    fit = purrr::map2(
      response,
      diversity,
      ~ fit_one_dbmem_model(
        data = all_cities,
        response = .x,
        diversity = .y,
        max_dbmem = max_dbmem,
        alpha_moran = alpha_moran,
        moran_nsim = moran_nsim,
        min_group_sites = min_group_sites
      )
    )
  )

#8.10 Context-specific simple slopes ####

# -------------------------------------------------------------------------
# This section estimates diversity–stability slopes separately for
# Urban and Rural sites using the final dbMEM-corrected models.
#
# Interpretation:
# - The diversity × urban_context interaction tests whether Urban and Rural
#   slopes differ from each other.
# - These simple slopes test whether the diversity–stability relationship
#   is significantly different from zero within Urban and Rural sites.
# -------------------------------------------------------------------------


# ---- Significance code ----

sig_code <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.1  ~ "*",
    TRUE ~ "ns"
  )
}


# ---- Helper: extract simple slopes from one model ----

get_context_slopes <- function(response, diversity, fit) {
  
  mod <- fit$model
  
  em <- emmeans::emtrends(
    mod,
    specs = "urban_context",
    var = diversity
  )
  
  out <- summary(
    em,
    infer = c(TRUE, TRUE)
  ) %>%
    tibble::as_tibble()
  
  # emtrends creates a column named, for example:
  # shannon_diversity.trend, species_richness.trend, etc.
  trend_col <- grep("\\.trend$", names(out), value = TRUE)
  
  if (length(trend_col) != 1) {
    stop(
      "Could not identify trend column for model: ",
      response,
      " ~ ",
      diversity
    )
  }
  
  out %>%
    transmute(
      response = .env$response,
      diversity = .env$diversity,
      urban_context = urban_context,
      slope = .data[[trend_col]],
      std.error = SE,
      df = df,
      conf.low = lower.CL,
      conf.high = upper.CL,
      statistic = t.ratio,
      p.value = p.value,
      sig = sig_code(p.value)
    )
}


# ---- Calculate simple slopes for all dbMEM models ----

glm_dbmem_simple_slopes <- glm_dbmem_fits %>%
  transmute(
    slopes = purrr::pmap(
      list(response, diversity, fit),
      get_context_slopes
    )
  ) %>%
  tidyr::unnest(slopes) %>%
  left_join(
    glm_dbmem_summary %>%
      select(
        response,
        diversity,
        n_dbmem,
        moran_corrected,
        n_significant_positive_moran
      ),
    by = c("response", "diversity")
  )


# ---- Wide version for easier inspection ----

glm_dbmem_simple_slopes_wide <- glm_dbmem_simple_slopes %>%
  select(
    response,
    diversity,
    urban_context,
    slope,
    std.error,
    p.value,
    sig
  ) %>%
  mutate(
    urban_context = as.character(urban_context)
  ) %>%
  tidyr::pivot_wider(
    names_from = urban_context,
    values_from = c(
      slope,
      std.error,
      p.value,
      sig
    ),
    names_glue = "{.value}_{urban_context}"
  )


# ---- Print results ----

message("\n==== Context-specific simple slopes ====")

glm_dbmem_simple_slopes %>%
  mutate(
    slope = round(slope, 4),
    std.error = round(std.error, 4),
    conf.low = round(conf.low, 4),
    conf.high = round(conf.high, 4),
    statistic = round(statistic, 3),
    p.value = signif(p.value, 3)
  ) %>%
  arrange(response, diversity, urban_context) %>%
  print(n = 100)


message("\n==== Context-specific simple slopes, wide format ====")

glm_dbmem_simple_slopes_wide %>%
  mutate(
    across(
      where(is.numeric),
      ~ signif(.x, 3)
    )
  ) %>%
  arrange(response, diversity) %>%
  print(n = 100)

#8.11 ANOVA-style tables for dbMEM models ####

# -------------------------------------------------------------------------
# This section creates ANOVA-style tables for the dbMEM-corrected GLMs.
# It reports Df, Sum Sq, Mean Sq, F value and p-value for the focal terms.
#
# Because the models include interactions and dbMEM spatial filters, we use
# Type III tests from car::Anova after refitting models with sum contrasts.
# dbMEM terms are included in the models but not reported individually,
# because they are spatial filters rather than ecological predictors.
# -------------------------------------------------------------------------

library(car)

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

bold_sig <- function(x, p) {
  ifelse(
    !is.na(p) & p < 0.05,
    paste0("**", x, "**"),
    x
  )
}

get_anova_dbmem <- function(response, diversity, fit) {
  
  mod_formula <- formula(fit$model)
  model_data  <- fit$data
  
  # Refit with sum contrasts for type III tests
  mod_refit <- lm(
    mod_formula,
    data = model_data,
    contrasts = list(
      urban_context = contr.sum,
      city = contr.sum
    )
  )
  
  anova_tbl <- car::Anova(
    mod_refit,
    type = 3
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term_raw") %>%
    tibble::as_tibble()
  
  interaction_term <- paste0(diversity, ":urban_context")
  
  out <- anova_tbl %>%
    filter(
      term_raw %in% c(
        diversity,
        "urban_context",
        "city",
        interaction_term
      )
    ) %>%
    mutate(
      response = response,
      diversity = diversity,
      term = case_when(
        term_raw == diversity ~ diversity,
        term_raw == "urban_context" ~ "Urban context",
        term_raw == "city" ~ "Study system",
        term_raw == interaction_term ~ paste0(diversity, " × urban context"),
        TRUE ~ term_raw
      ),
      Df = Df,
      Sum_Sq = `Sum Sq`,
      Mean_Sq = `Sum Sq` / Df,
      F_value = `F value`,
      p_value = `Pr(>F)`
    ) %>%
    select(
      response,
      diversity,
      term,
      Df,
      Sum_Sq,
      Mean_Sq,
      F_value,
      p_value
    )
  
  resid_df <- df.residual(mod_refit)
  resid_ss <- sum(residuals(mod_refit)^2)
  resid_ms <- resid_ss / resid_df
  
  resid_row <- tibble(
    response = response,
    diversity = diversity,
    term = "Residuals",
    Df = resid_df,
    Sum_Sq = resid_ss,
    Mean_Sq = resid_ms,
    F_value = NA_real_,
    p_value = NA_real_
  )
  
  bind_rows(out, resid_row)
}

glm_dbmem_anova_tables <- glm_dbmem_fits %>%
  mutate(
    anova = purrr::pmap(
      list(response, diversity, fit),
      get_anova_dbmem
    )
  ) %>%
  select(anova) %>%
  tidyr::unnest(anova) %>%
  left_join(
    glm_dbmem_summary %>%
      select(
        response,
        diversity,
        n_dbmem,
        n_obs,
        r.squared,
        adj.r.squared,
        moran_corrected,
        min_moran_p
      ),
    by = c("response", "diversity")
  )

message("\n==== ANOVA-style tables for dbMEM models ====")

glm_dbmem_anova_tables %>%
  mutate(
    Sum_Sq = round(Sum_Sq, 3),
    Mean_Sq = round(Mean_Sq, 3),
    F_value = round(F_value, 3),
    p_value = format_p(p_value)
  ) %>%
  arrange(response, diversity, match(term, c(
    diversity,
    "Urban context",
    "Study system",
    paste0(diversity, " × urban context"),
    "Residuals"
  ))) %>%
  print(n = 200)



#8.12 Extract and save results ####

glm_dbmem_summary <- purrr::map_dfr(
  glm_dbmem_fits$fit,
  \(x) x$summary
)

glm_dbmem_coefficients <- purrr::map_dfr(
  glm_dbmem_fits$fit,
  \(x) x$coefficients
) %>%
  left_join(
    glm_dbmem_summary %>%
      select(
        response,
        diversity,
        moran_corrected,
        n_significant_positive_moran,
        significant_moran_groups
      ),
    by = c("response", "diversity")
  )

glm_dbmem_focal_terms <- glm_dbmem_coefficients %>%
  filter(!is_spatial_filter)

glm_dbmem_focal_terms_corrected <- glm_dbmem_focal_terms %>%
  filter(moran_corrected)

glm_dbmem_interactions_corrected <- glm_dbmem_focal_terms_corrected %>%
  filter(stringr::str_detect(term, ":urban_context")) %>%
  mutate(
    significant = p.value < 0.05
  )

glm_dbmem_moran <- purrr::map_dfr(
  glm_dbmem_fits$fit,
  \(x) x$moran
)


# ---- Print final checks ----

message("\n==== Final dbMEM model summary ====")

glm_dbmem_summary %>%
  select(
    response,
    diversity,
    n_dbmem,
    moran_corrected,
    n_significant_positive_moran,
    significant_moran_groups,
    min_moran_p
  ) %>%
  arrange(response, diversity) %>%
  print(n = 100)

message("\n==== Final dbMEM interaction results ====")

glm_dbmem_interactions_corrected %>%
  select(
    response,
    diversity,
    n_dbmem,
    term,
    estimate,
    std.error,
    statistic,
    p.value,
    moran_corrected,
    significant
  ) %>%
  mutate(
    estimate = round(estimate, 4),
    std.error = round(std.error, 4),
    statistic = round(statistic, 3),
    p.value = signif(p.value, 3)
  ) %>%
  arrange(response, p.value) %>%
  print(n = 100)


# ---- Save outputs ----

write.csv(
  glm_dbmem_coefficients,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_diversity_stability_coefficients_all_terms.csv"),
  row.names = FALSE
)

write.csv(
  glm_dbmem_focal_terms,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_diversity_stability_focal_terms.csv"),
  row.names = FALSE
)

write.csv(
  glm_dbmem_focal_terms_corrected,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_focal_terms_only_moran_corrected.csv"),
  row.names = FALSE
)

write.csv(
  glm_dbmem_interactions_corrected,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_interactions_only_moran_corrected.csv"),
  row.names = FALSE
)

write.csv(
  glm_dbmem_summary,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_diversity_stability_model_summary.csv"),
  row.names = FALSE
)

write.csv(
  glm_dbmem_moran,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_diversity_stability_moran_by_group.csv"),
  row.names = FALSE
)

write.csv(
  glm_dbmem_simple_slopes,
  file.path(
    output_dir,
    "results",
    "GLM_dbMEM_by_city_simple_slopes_by_context.csv"
  ),
  row.names = FALSE
)

write.csv(
  glm_dbmem_simple_slopes_wide,
  file.path(
    output_dir,
    "results",
    "GLM_dbMEM_by_city_simple_slopes_by_context_wide.csv"
  ),
  row.names = FALSE
)

write.csv(
  glm_dbmem_anova_tables,
  file.path(
    output_dir,
    "results",
    "GLM_dbMEM_by_city_ANOVA_style_tables.csv"
  ),
  row.names = FALSE
)

saveRDS(
  glm_dbmem_fits,
  file.path(output_dir, "results", "GLM_dbMEM_by_city_diversity_stability_model_objects.rds")
)

message("\n==== dbMEM GLM analysis complete ====")
message("Spatial filters calculated by city.")
message("Moran's I assessed by city × urban_context.")
message("Saved outputs in output/results/")


#9. Standardized pooled diversity–stability plots using dbMEM model predictions ####

#9.1 Dataset

dir.create(
  file.path(output_dir, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

all_cities_std_pooled <- all_cities %>%
  filter(!is.na(urban_context))


#9.2 Plot settings

plot_colors <- c(
  "Urban" = "#D55E00",
  "Rural" = "#0072B2"
)

#9.3 Helper: extract fitted model object

get_dbmem_fit <- function(response_name, diversity_name) {
  
  fit_obj <- glm_dbmem_fits %>%
    filter(
      .data$response == .env$response_name,
      .data$diversity == .env$diversity_name
    ) %>%
    pull(fit)
  
  if (length(fit_obj) != 1) {
    stop(
      "Could not find unique model for: ",
      response_name,
      " ~ ",
      diversity_name
    )
  }
  
  fit_obj[[1]]
}

#9.4 Helper: city-averaged predictions from dbMEM model

predict_city_average <- function(fit_object, diversity, n_points = 100) {
  
  mod <- fit_object$model
  model_data <- fit_object$data
  
  x_range <- range(model_data[[diversity]], na.rm = TRUE)
  x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
  
  city_levels <- levels(droplevels(model_data$city))
  context_levels <- levels(droplevels(model_data$urban_context))
  
  city_weights <- tibble(
    city = factor(city_levels, levels = levels(model_data$city)),
    weight = 1 / length(city_levels)
  )
  
  dbmem_terms <- grep("^dbMEM_", names(model_data), value = TRUE)
  
  pred_grid <- expand.grid(
    x = x_seq,
    urban_context = context_levels,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
  
  beta <- coef(mod)
  V <- vcov(mod)
  tt <- delete.response(terms(mod))
  
  pred_out <- purrr::map_dfr(seq_len(nrow(pred_grid)), function(i) {
    
    this_x <- pred_grid$x[i]
    this_context <- pred_grid$urban_context[i]
    
    nd <- tibble(
      city = factor(city_levels, levels = levels(model_data$city)),
      urban_context = factor(
        this_context,
        levels = levels(model_data$urban_context)
      )
    )
    
    nd[[diversity]] <- this_x
    
    # Spatial filters are set to 0 for the plotted marginal prediction.
    if (length(dbmem_terms) > 0) {
      for (z in dbmem_terms) {
        nd[[z]] <- 0
      }
    }
    
    X <- model.matrix(
      tt,
      data = nd,
      contrasts.arg = mod$contrasts
    )
    
    missing_cols <- setdiff(names(beta), colnames(X))
    
    if (length(missing_cols) > 0) {
      X_missing <- matrix(
        0,
        nrow = nrow(X),
        ncol = length(missing_cols)
      )
      colnames(X_missing) <- missing_cols
      X <- cbind(X, X_missing)
    }
    
    X <- X[, names(beta), drop = FALSE]
    
    weights <- city_weights$weight[
      match(as.character(nd$city), as.character(city_weights$city))
    ]
    
    X_avg <- matrix(
      colSums(X * weights),
      nrow = 1
    )
    
    colnames(X_avg) <- colnames(X)
    
    predicted <- as.numeric(X_avg %*% beta)
    se <- sqrt(as.numeric(X_avg %*% V %*% t(X_avg)))
    
    tibble(
      x = this_x,
      urban_context = this_context,
      predicted = predicted,
      std.error = se,
      conf.low = predicted - 1.96 * se,
      conf.high = predicted + 1.96 * se
    )
  })
  
  pred_out
}

#9.5 Helper: one panel with dbMEM model predictions

make_std_modelpred_panel <- function(
    data,
    response,
    diversity,
    xlab,
    ylab,
    show_x_axis = TRUE,
    show_y_axis = TRUE
) {
  
  fit_object <- get_dbmem_fit(
    response_name = response,
    diversity_name = diversity
  )
  
  pred <- predict_city_average(
    fit_object = fit_object,
    diversity = diversity,
    n_points = 100
  )
  
  plot_data <- data %>%
    filter(
      !is.na(.data[[response]]),
      !is.na(.data[[diversity]]),
      !is.na(urban_context)
    )
  
  p <- ggplot() +
    
    geom_point(
      data = plot_data,
      aes(
        x = .data[[diversity]],
        y = .data[[response]],
        color = urban_context
      ),
      alpha = 0.16,
      size = 0.45
    ) +
    
    geom_ribbon(
      data = pred,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = urban_context
      ),
      alpha = 0.13,
      color = NA
    ) +
    
    geom_line(
      data = pred,
      aes(
        x = x,
        y = predicted,
        color = urban_context
      ),
      linewidth = 1.0
    ) +
    
    scale_x_continuous(
      breaks = function(x) {
        seq(
          floor(min(x, na.rm = TRUE) / 2) * 2,
          ceiling(max(x, na.rm = TRUE) / 2) * 2,
          by = 2
        )
      }
    ) +
    
    scale_y_continuous(
      breaks = seq(-2, 2, by = 2)
    ) +
    
    coord_cartesian(
      ylim = c(-2, 2)
    ) +
    
    scale_color_manual(
      values = plot_colors,
      name = "Urban context"
    ) +
    
    scale_fill_manual(
      values = plot_colors,
      guide = "none"
    ) +
    
    labs(
      x = if (show_x_axis) xlab else NULL,
      y = if (show_y_axis) ylab else NULL
    ) +
    
    theme_classic(
      base_family = "Garamond",
      base_size = 15
    ) +
    
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.75
      ),
      axis.line = element_blank(),
      legend.position = "bottom",
      
      axis.title.x = if (show_x_axis) {
        element_text(family = "Garamond", size = 15)
      } else {
        element_blank()
      },
      
      axis.title.y = if (show_y_axis) {
        element_text(family = "Garamond", size = 15)
      } else {
        element_blank()
      },
      
      axis.text = element_text(
        family = "Garamond",
        size = 11
      ),
      
      aspect.ratio = 1,
      plot.margin = margin(t = 1, r = 2, b = 1, l = 2, unit = "pt")
    )
  
  if (!show_x_axis) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  if (!show_y_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  p
}


#9.6 Generate panels

# Row 1: Community stability

p_cs_sr <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "community_stability",
  diversity = "species_richness",
  xlab = "Species richness",
  ylab = "Community stability",
  show_x_axis = FALSE,
  show_y_axis = TRUE
)

p_cs_sd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "community_stability",
  diversity = "shannon_diversity",
  xlab = "Species diversity",
  ylab = "Community stability",
  show_x_axis = FALSE,
  show_y_axis = FALSE
)

p_cs_fd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "community_stability",
  diversity = "FDis",
  xlab = "Functional diversity",
  ylab = "Community stability",
  show_x_axis = FALSE,
  show_y_axis = FALSE
)

p_cs_mpd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "community_stability",
  diversity = "MPD",
  xlab = "Phylogenetic diversity",
  ylab = "Community stability",
  show_x_axis = FALSE,
  show_y_axis = FALSE
)

# Row 2: Species asynchrony

p_sa_sr <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "species_asynchrony",
  diversity = "species_richness",
  xlab = "Species richness",
  ylab = "Species asynchrony",
  show_x_axis = FALSE,
  show_y_axis = TRUE
)

p_sa_sd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "species_asynchrony",
  diversity = "shannon_diversity",
  xlab = "Species diversity",
  ylab = "Species asynchrony",
  show_x_axis = FALSE,
  show_y_axis = FALSE
)

p_sa_fd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "species_asynchrony",
  diversity = "FDis",
  xlab = "Functional diversity",
  ylab = "Species asynchrony",
  show_x_axis = FALSE,
  show_y_axis = FALSE
)

p_sa_mpd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "species_asynchrony",
  diversity = "MPD",
  xlab = "Phylogenetic diversity",
  ylab = "Species asynchrony",
  show_x_axis = FALSE,
  show_y_axis = FALSE
)

# Row 3: Population stability

p_ps_sr <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "wm_population_stability",
  diversity = "species_richness",
  xlab = "Species richness",
  ylab = "Population stability",
  show_x_axis = TRUE,
  show_y_axis = TRUE
)

p_ps_sd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "wm_population_stability",
  diversity = "shannon_diversity",
  xlab = "Species diversity",
  ylab = "Population stability",
  show_x_axis = TRUE,
  show_y_axis = FALSE
)

p_ps_fd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "wm_population_stability",
  diversity = "FDis",
  xlab = "Functional diversity",
  ylab = "Population stability",
  show_x_axis = TRUE,
  show_y_axis = FALSE
)

p_ps_mpd <- make_std_modelpred_panel(
  all_cities_std_pooled,
  response = "wm_population_stability",
  diversity = "MPD",
  xlab = "Phylogenetic diversity",
  ylab = "Population stability",
  show_x_axis = TRUE,
  show_y_axis = FALSE
)


#9.7 Combine figure

std_modelpred_3x4_plot <- (
  (p_cs_sr | p_cs_sd | p_cs_fd | p_cs_mpd) /
    (p_sa_sr | p_sa_sd | p_sa_fd | p_sa_mpd) /
    (p_ps_sr | p_ps_sd | p_ps_fd | p_ps_mpd)
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Garamond",
      size = 15
    ),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
  )

std_modelpred_3x4_plot


#9.8 Save figure

ggsave(
  filename = file.path(
    output_dir,
    "figures",
    "fig_standardized_modelpred_3stability_x_4diversity.png"
  ),
  plot = std_modelpred_3x4_plot,
  width = 8,
  height = 7,
  dpi = 600,
  bg = "white"
)


#10 Supplementary region-specific raw relationships ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

dir.create(
  file.path(output_dir, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

plot_colors <- c(
  "Urban" = "#D55E00",
  "Rural" = "#0072B2"
)

stability_vars <- c(
  "community_stability",
  "species_asynchrony",
  "wm_population_stability"
)

diversity_vars <- c(
  "species_richness",
  "shannon_diversity",
  "FDis",
  "MPD"
)

stability_labels <- c(
  community_stability = "Community stability",
  species_asynchrony = "Species asynchrony",
  wm_population_stability = "Population stability"
)

diversity_labels <- c(
  species_richness = "Species richness",
  shannon_diversity = "Species diversity",
  FDis = "Functional diversity",
  MPD = "Phylogenetic diversity"
)

city_levels <- c("LND", "RND", "BCN")

# RAW VALUES
plot_data_s18 <- all_cities_check %>%
  select(
    city,
    urban_context,
    all_of(stability_vars),
    all_of(diversity_vars)
  ) %>%
  pivot_longer(
    cols = all_of(stability_vars),
    names_to = "stability_metric",
    values_to = "stability_value"
  ) %>%
  pivot_longer(
    cols = all_of(diversity_vars),
    names_to = "diversity_metric",
    values_to = "diversity_value"
  ) %>%
  mutate(
    city = factor(city, levels = city_levels),
    urban_context = factor(urban_context, levels = c("Urban", "Rural"))
  ) %>%
  filter(
    !is.na(stability_value),
    !is.na(diversity_value),
    !is.na(city),
    !is.na(urban_context)
  )

make_s18_panel <- function(stab, div, city_i, show_y, show_x, show_title) {
  
  df <- plot_data_s18 %>%
    filter(
      stability_metric == stab,
      diversity_metric == div,
      city == city_i
    )
  
  ggplot(
    df,
    aes(
      x = diversity_value,
      y = stability_value,
      color = urban_context,
      fill = urban_context
    )
  ) +
    geom_point(
      alpha = 0.18,
      size = 0.45
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.7,
      alpha = 0.14
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_x_continuous(
      breaks = scales::breaks_width(2)
    ) +
    scale_y_continuous(
      breaks = scales::breaks_width(2)
    ) +
    labs(
      title = if (show_title) as.character(city_i) else NULL,
      x = if (show_x) diversity_labels[[div]] else NULL,
      y = if (show_y) stability_labels[[stab]] else NULL
    ) +
    theme_classic(
      base_family = "Garamond",
      base_size = 12
    ) +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.4
      ),
      axis.line = element_blank(),
      
      plot.title = element_text(
        family = "Garamond",
        size = 16,
        hjust = 0.5,
        margin = margin(b = 4)
      ),
      
      axis.title.x = element_text(
        family = "Garamond",
        size = 15,
        margin = margin(t = 6)
      ),
      axis.title.y = element_text(
        family = "Garamond",
        size = 15,
        margin = margin(r = 6)
      ),
      
      axis.text = element_text(
        family = "Garamond",
        size = 10
      ),
      axis.ticks = element_line(linewidth = 0.25),
      
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_header_plot <- function(label = NULL) {
  
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      family = "Garamond",
      size = 5
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
}

make_empty_plot <- function() {
  ggplot() + theme_void()
}

# Header row: city names below legend
header_list <- list()

for (div in diversity_vars) {
  for (city_i in city_levels) {
    header_list[[length(header_list) + 1]] <- make_header_plot(city_i)
  }
  if (div != diversity_vars[length(diversity_vars)]) {
    header_list[[length(header_list) + 1]] <- make_empty_plot()
  }
}

# Main panels with spacer columns between diversity metrics
panel_list <- list()

for (stab_i in seq_along(stability_vars)) {
  
  stab <- stability_vars[stab_i]
  
  for (div in diversity_vars) {
    
    for (city_i in city_levels) {
      
      panel_list[[length(panel_list) + 1]] <- make_s18_panel(
        stab = stab,
        div = div,
        city_i = city_i,
        show_y = div == diversity_vars[1] & city_i == city_levels[1],
        show_x = stab == stability_vars[length(stability_vars)] & city_i == "RND",
        show_title = FALSE
      )
    }
    
    if (div != diversity_vars[length(diversity_vars)]) {
      panel_list[[length(panel_list) + 1]] <- make_empty_plot()
    }
  }
}

legend_plot <- ggplot(
  plot_data_s18,
  aes(
    x = diversity_value,
    y = stability_value,
    color = urban_context,
    fill = urban_context
  )
) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.7,
    alpha = 0.14
  ) +
  scale_color_manual(values = plot_colors) +
  scale_fill_manual(values = plot_colors) +
  theme_void(base_family = "Garamond") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Garamond",
      size = 12
    )
  )

legend_only <- cowplot::get_legend(legend_plot)

header_row <- wrap_plots(
  header_list,
  ncol = 15,
  widths = c(1, 1, 1, 0.25, 1, 1, 1, 0.25, 1, 1, 1, 0.25, 1, 1, 1)
)

panel_grid <- wrap_plots(
  panel_list,
  ncol = 15,
  widths = c(1, 1, 1, 0.25, 1, 1, 1, 0.25, 1, 1, 1, 0.25, 1, 1, 1)
)

fig_s6_raw_region_relationships <- 
  cowplot::ggdraw(legend_only) /
  header_row /
  panel_grid +
  plot_layout(
    heights = c(0.35, 0.25, 6)
  )

fig_s6_raw_region_relationships

ggsave(
  filename = file.path(
    output_dir,
    "figures",
    "Fig_S6.png"
  ),
  plot = fig_s6_raw_region_relationships,
  width = 16,
  height = 8,
  dpi = 600,
  bg = "white"
)
