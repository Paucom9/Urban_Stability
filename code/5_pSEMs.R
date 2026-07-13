# ============================================================================ #
#   5_pSEMs.R                                                                  #
#   Author: Pau Colom                                                          #
#                                                                              #
#   It runs spatially corrected piecewise structural equation models for       #
#   London, Randstad, and Barcelona at 1-, 2-, and 5-km landscape buffers      # 
#                                                                              #
# ============================================================================ #


# --------------------------------------------------------------------------- #
# 1. Libraries
# --------------------------------------------------------------------------- #

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(here)
library(sf)
library(piecewiseSEM)
library(spdep)
library(adespatial)
library(nlme)

# --------------------------------------------------------------------------- #
# 2. Folders and settings
# --------------------------------------------------------------------------- #

input_dir   <- here("input", "data")
output_dir  <- here("output")
data_dir    <- file.path(output_dir, "data")
results_dir <- file.path(output_dir, "results", "SEMs")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

regions <- c("LND", "RND", "BCN")
buffers <- c(1000, 2000, 5000)

max_mem     <- 120
k_neigh     <- 8
alpha_moran <- 0.05
max_select  <- 90

# use GLS spatial correlation for the species richness equation.
# This is intended to avoid over-parameterizing the Randstad richness equation
# with many MEMs while still modelling the ecological path BS/LD -> SR.
use_gls_for_richness <- TRUE

# --------------------------------------------------------------------------- #
# 3. Helper functions: IDs, data loading, coordinates
# --------------------------------------------------------------------------- #

clean_bcn_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[is.na(x) | x == ""] <- NA_character_
  x <- gsub("^ES-CTBMS\\.", "", x)
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    grepl("^ES_uBMS", x) ~ x,
    TRUE ~ paste0("ES-CTBMS.", x)
  )
}

clean_site_id <- function(x, region) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[is.na(x) | x == ""] <- NA_character_
  if (region == "BCN") {
    x <- clean_bcn_id(x)
  }
  x
}

prep_landscape <- function(built, landdiv, region) {
  
  built_clean <- built %>%
    mutate(
      SITE_ID = clean_site_id(SITE_ID, region),
      CONTEXT = stringr::str_to_lower(stringr::str_trim(as.character(CONTEXT)))
    ) %>%
    filter(!is.na(SITE_ID)) %>%
    distinct(SITE_ID, .keep_all = TRUE) %>%
    rename(CONTEXT_built = CONTEXT)
  
  land_clean <- landdiv %>%
    mutate(
      SITE_ID = clean_site_id(SITE_ID, region),
      CONTEXT = stringr::str_to_lower(stringr::str_trim(as.character(CONTEXT)))
    ) %>%
    filter(!is.na(SITE_ID)) %>%
    distinct(SITE_ID, .keep_all = TRUE) %>%
    rename(CONTEXT_land = CONTEXT)
  
  full_join(
    built_clean,
    land_clean,
    by = "SITE_ID"
  ) %>%
    mutate(
      CONTEXT = dplyr::coalesce(CONTEXT_built, CONTEXT_land),
      urban_context = as.numeric(CONTEXT == "inside")
    ) %>%
    dplyr::select(
      -CONTEXT_built,
      -CONTEXT_land,
      -CONTEXT
    )
}

join_region_data <- function(ds_region, landscape_region, region) {
  
  out <- ds_region %>%
    left_join(landscape_region, by = "SITE_ID")
  
  message("\n---- ", region, " merge check ----")
  message("Sites in diversity/stability data: ", n_distinct(ds_region$SITE_ID))
  message("Sites after join: ", n_distinct(out$SITE_ID))
  message("Sites without urban_context: ", sum(is.na(out$urban_context)))
  
  if (region == "BCN") {
    message(
      "uBMS in diversity/stability data: ",
      ds_region %>% filter(grepl("^ES_uBMS", SITE_ID)) %>% nrow()
    )
    message(
      "uBMS after join with urban_context: ",
      out %>% filter(grepl("^ES_uBMS", SITE_ID), !is.na(urban_context)) %>% nrow()
    )
  }
  
  out %>%
    filter(!is.na(urban_context))
}

load_site_coordinates <- function(input_dir) {
  
  # ---- eBMS coordinates: UK, NL and BCN-CBMS ----
  # transect_lon / transect_lat are already projected coordinates in EPSG:3035.
  ebms_coords_raw <- read.csv(
    file.path(input_dir, "ebms_transect_coord.csv")
  )
  
  ebms_coords <- ebms_coords_raw %>%
    transmute(
      SITE_ID_raw = stringr::str_trim(as.character(transect_id)),
      bms_id = stringr::str_trim(as.character(bms_id)),
      x = as.numeric(transect_lon),
      y = as.numeric(transect_lat)
    ) %>%
    filter(!is.na(SITE_ID_raw), !is.na(x), !is.na(y)) %>%
    mutate(
      city = case_when(
        bms_id == "UKBMS" ~ "LND",
        bms_id == "NLBMS" ~ "RND",
        bms_id == "ES-CTBMS" ~ "BCN",
        TRUE ~ NA_character_
      ),
      SITE_ID = case_when(
        city == "BCN" ~ clean_site_id(SITE_ID_raw, "BCN"),
        city == "LND" ~ clean_site_id(SITE_ID_raw, "LND"),
        city == "RND" ~ clean_site_id(SITE_ID_raw, "RND"),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(city), !is.na(SITE_ID)) %>%
    dplyr::select(city, SITE_ID, x, y) %>%
    distinct(city, SITE_ID, .keep_all = TRUE)
  
  # ---- uBMS coordinates for Barcelona ----
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
      SITE_ID = clean_site_id(SITE_ID, "BCN"),
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
      city = "BCN",
      x = ubms_xy[, 1],
      y = ubms_xy[, 2]
    ) %>%
    dplyr::select(city, SITE_ID, x, y) %>%
    distinct(city, SITE_ID, .keep_all = TRUE)
  
  bind_rows(
    ebms_coords,
    ubms_coords
  ) %>%
    distinct(city, SITE_ID, .keep_all = TRUE)
}

load_region_sem_raw <- function(region, data_dir, site_coords) {
  
  suffix <- tolower(region)
  
  ds_region <- read.csv(
    file.path(data_dir, paste0("diversity_stability_", suffix, ".csv"))
  ) %>%
    mutate(
      SITE_ID = clean_site_id(SITE_ID, region)
    )
  
  landscape_region <- prep_landscape(
    built = read.csv(file.path(data_dir, paste0("built_up_", suffix, ".csv"))),
    landdiv = read.csv(file.path(data_dir, paste0("land_diversity_", suffix, ".csv"))),
    region = region
  )
  
  out <- join_region_data(
    ds_region = ds_region,
    landscape_region = landscape_region,
    region = region
  ) %>%
    mutate(city = region) %>%
    left_join(site_coords, by = c("city", "SITE_ID"))
  
  missing_coords <- out %>%
    filter(is.na(x) | is.na(y)) %>%
    distinct(city, SITE_ID, urban_context)
  
  message("\n==== ", region, " coordinate check ====")
  message("Sites missing coordinates: ", nrow(missing_coords))
  if (nrow(missing_coords) > 0) {
    print(tibble::as_tibble(missing_coords), n = 100)
  }
  
  out <- out %>%
    filter(!is.na(x), !is.na(y))
  
  message("\n==== ", region, " data after coordinate join ====")
  message("N sites: ", nrow(out))
  print(table(out$urban_context, useNA = "ifany"))
  
  out
}

# --------------------------------------------------------------------------- #
# 4. Helper functions: formulas, VIF, dbMEM, Moran
# --------------------------------------------------------------------------- #

make_formula <- function(response, predictors, env = parent.frame()) {
  f <- reformulate(
    termlabels = predictors,
    response = response
  )
  environment(f) <- env
  f
}

clean_empty_names <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  names(x)[names(x) == ""] <- "signif"
  tibble::as_tibble(x, .name_repair = "unique")
}

split_mem_string <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") {
    character(0)
  } else {
    stringr::str_split(x, " \\+ ")[[1]]
  }
}

manual_vif <- function(df, vars) {
  purrr::map_dfr(
    vars,
    function(v) {
      other_vars <- setdiff(vars, v)
      if (length(other_vars) == 0) {
        return(
          tibble::tibble(
            variable = v,
            R2 = NA_real_,
            VIF = NA_real_
          )
        )
      }
      f <- make_formula(v, other_vars)
      mod <- lm(f, data = df)
      r2 <- summary(mod)$r.squared
      tibble::tibble(
        variable = v,
        R2 = r2,
        VIF = 1 / (1 - r2)
      )
    }
  )
}


# Diagnose which sites are lost before fitting the pSEM.
# This is important for BCN, where the merged dataset can contain 64 sites,
# but complete-case filtering for the SEM may retain fewer sites if any of the
# required stability, diversity, landscape or coordinate variables are missing.
diagnose_psem_complete_cases <- function(df_input, region, buffer, results_dir) {
  diagnostic_vars <- setdiff(
    names(df_input),
    c("SITE_ID", "city")
  )
  
  missing_matrix <- is.na(df_input[, diagnostic_vars, drop = FALSE])
  complete_case <- rowSums(missing_matrix) == 0
  
  site_retention <- df_input %>%
    dplyr::mutate(
      region = region,
      buffer = buffer,
      complete_case = complete_case,
      .before = 1
    ) %>%
    dplyr::count(region, buffer, urban_context, complete_case, name = "n_sites") %>%
    dplyr::arrange(region, buffer, urban_context, complete_case)
  
  missing_by_variable <- tibble::tibble(
    region = region,
    buffer = buffer,
    variable = diagnostic_vars,
    n_missing = colSums(missing_matrix),
    n_available = nrow(df_input),
    prop_missing = n_missing / n_available
  ) %>%
    dplyr::arrange(dplyr::desc(n_missing), variable)
  
  dropped_sites <- df_input %>%
    dplyr::mutate(
      region = region,
      buffer = buffer,
      n_missing = rowSums(missing_matrix),
      missing_vars = apply(
        missing_matrix,
        1,
        function(z) paste(diagnostic_vars[z], collapse = "; ")
      ),
      .before = 1
    ) %>%
    dplyr::filter(n_missing > 0) %>%
    dplyr::select(
      region,
      buffer,
      SITE_ID,
      city,
      urban_context,
      n_missing,
      missing_vars
    ) %>%
    dplyr::arrange(region, buffer, urban_context, SITE_ID)
  
  message("Sites available before complete-case filter: ", nrow(df_input))
  message("Complete-case sites available for pSEM: ", sum(complete_case))
  message("Sites dropped by complete-case filter: ", sum(!complete_case))
  
  if (any(missing_by_variable$n_missing > 0)) {
    message("Variables causing site loss:")
    print(
      missing_by_variable %>%
        dplyr::filter(n_missing > 0),
      n = 100
    )
  }
  
  if (nrow(dropped_sites) > 0) {
    message("Dropped sites by urban_context:")
    print(
      dropped_sites %>%
        dplyr::count(region, buffer, urban_context, name = "n_dropped")
    )
    
    utils::write.csv(
      dropped_sites,
      file.path(
        results_dir,
        paste0("Dropped_sites_", region, "_buffer_", buffer, ".csv")
      ),
      row.names = FALSE
    )
  }
  
  list(
    site_retention = site_retention,
    missing_by_variable = missing_by_variable,
    dropped_sites = dropped_sites
  )
}


signif_code <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

extract_lm_path_coefs <- function(mod, model_data, response, region, buffer,
                                  built_var, landdiv_var) {
  
  # Works for both lm() and nlme::gls() component models.
  # For gls(), coefficients are stored in summary(mod)$tTable.
  # For lm(), coefficients are stored in summary(mod)$coefficients.
  if (inherits(mod, "gls")) {
    coef_mat <- summary(mod)$tTable
    coef_tbl <- as.data.frame(coef_mat, check.names = FALSE) %>%
      tibble::rownames_to_column("Predictor") %>%
      dplyr::rename(
        Estimate = Value,
        Std.Error = Std.Error,
        Crit.Value = `t-value`,
        P.Value = `p-value`
      ) %>%
      dplyr::mutate(
        DF = mod$dims$N - length(stats::coef(mod))
      )
  } else {
    coef_mat <- summary(mod)$coefficients
    coef_tbl <- as.data.frame(coef_mat, check.names = FALSE) %>%
      tibble::rownames_to_column("Predictor") %>%
      dplyr::rename(
        Estimate = Estimate,
        Std.Error = `Std. Error`,
        Crit.Value = `t value`,
        P.Value = `Pr(>|t|)`
      ) %>%
      dplyr::mutate(
        DF = stats::df.residual(mod)
      )
  }
  
  if (is.null(coef_tbl) || nrow(coef_tbl) == 0) {
    return(tibble::tibble())
  }
  
  coef_tbl <- coef_tbl %>%
    dplyr::filter(Predictor != "(Intercept)")
  
  if (nrow(coef_tbl) == 0) {
    return(tibble::tibble())
  }
  
  y_sd <- stats::sd(model_data[[response]], na.rm = TRUE)
  
  coef_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      .x_sd = if (Predictor %in% names(model_data)) {
        stats::sd(model_data[[Predictor]], na.rm = TRUE)
      } else {
        NA_real_
      },
      Std.Estimate = ifelse(
        !is.na(.x_sd) && !is.na(y_sd) && y_sd > 0,
        Estimate * .x_sd / y_sd,
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      region = region,
      buffer = buffer,
      Response = response,
      built_var = built_var,
      landdiv_var = landdiv_var,
      signif = signif_code(P.Value),
      extraction = ifelse(inherits(mod, "gls"), "manual_gls", "manual_lm")
    ) %>%
    dplyr::select(
      region,
      buffer,
      Response,
      Predictor,
      Estimate,
      Std.Error,
      DF,
      Crit.Value,
      P.Value,
      Std.Estimate,
      signif,
      built_var,
      landdiv_var,
      extraction
    )
}

extract_correlated_error <- function(mod1, mod2, response1, response2,
                                     region, buffer, built_var, landdiv_var) {
  
  res1 <- stats::residuals(mod1)
  res2 <- stats::residuals(mod2)
  ok <- stats::complete.cases(res1, res2)
  
  if (sum(ok) < 4) {
    return(tibble::tibble())
  }
  
  ct <- stats::cor.test(res1[ok], res2[ok])
  est <- unname(ct$estimate)
  pval <- ct$p.value
  
  tibble::tibble(
    region = region,
    buffer = buffer,
    Response = paste0("~~", response1),
    Predictor = paste0("~~", response2),
    Estimate = est,
    Std.Error = NA_real_,
    DF = unname(ct$parameter),
    Crit.Value = unname(ct$statistic),
    P.Value = pval,
    Std.Estimate = est,
    signif = signif_code(pval),
    built_var = built_var,
    landdiv_var = landdiv_var,
    extraction = "manual_residual_correlation"
  )
}

extract_raw_correlation <- function(model_data, var1, var2,
                                    region, buffer, built_var, landdiv_var) {
  ok <- stats::complete.cases(model_data[[var1]], model_data[[var2]])
  
  if (sum(ok) < 4) {
    return(tibble::tibble())
  }
  
  ct <- stats::cor.test(model_data[[var1]][ok], model_data[[var2]][ok])
  est <- unname(ct$estimate)
  pval <- ct$p.value
  
  tibble::tibble(
    region = region,
    buffer = buffer,
    Response = paste0("~~", var1),
    Predictor = paste0("~~", var2),
    Estimate = est,
    Std.Error = NA_real_,
    DF = unname(ct$parameter),
    Crit.Value = unname(ct$statistic),
    P.Value = pval,
    Std.Estimate = est,
    signif = signif_code(pval),
    built_var = built_var,
    landdiv_var = landdiv_var,
    extraction = "manual_raw_correlation"
  )
}

extract_component_path_coefs <- function(component_models, model_data, region,
                                         buffer, built_var, landdiv_var) {
  
  path_coefs <- purrr::imap_dfr(
    component_models,
    function(mod, response) {
      extract_lm_path_coefs(
        mod = mod,
        model_data = model_data,
        response = response,
        region = region,
        buffer = buffer,
        built_var = built_var,
        landdiv_var = landdiv_var
      )
    }
  )
  
  correlated_errors <- dplyr::bind_rows(
    extract_correlated_error(
      mod1 = component_models$species_asynchrony,
      mod2 = component_models$wm_population_stability,
      response1 = "species_asynchrony",
      response2 = "wm_population_stability",
      region = region,
      buffer = buffer,
      built_var = built_var,
      landdiv_var = landdiv_var
    ),
    extract_raw_correlation(
      model_data = model_data,
      var1 = "built_buffer",
      var2 = "landdiv_buffer",
      region = region,
      buffer = buffer,
      built_var = built_var,
      landdiv_var = landdiv_var
    )
  )
  
  dplyr::bind_rows(path_coefs, correlated_errors)
}

jitter_duplicate_xy <- function(df, x_col = "x", y_col = "y", amount = 1) {
  df %>%
    group_by(across(all_of(c(x_col, y_col)))) %>%
    mutate(
      .dup_id = row_number(),
      .dup_n = n(),
      x_dbmem = .data[[x_col]] + if_else(.dup_n > 1, (.dup_id - 1) * amount, 0),
      y_dbmem = .data[[y_col]] + if_else(.dup_n > 1, (.dup_id - 1) * amount, 0)
    ) %>%
    ungroup() %>%
    dplyr::select(-.dup_id, -.dup_n)
}

add_dbmems <- function(df, max_mem = 40, coords = c("x_dbmem", "y_dbmem")) {
  
  xy <- df %>%
    dplyr::select(all_of(coords)) %>%
    as.matrix()
  
  mem_obj <- adespatial::dbmem(
    xy,
    MEM.autocor = "positive",
    silent = TRUE
  )
  
  mem_mat <- as.data.frame(as.matrix(mem_obj))
  
  if (ncol(mem_mat) == 0) {
    warning("No positive dbMEMs were generated.")
    return(df)
  }
  
  n_keep <- min(max_mem, ncol(mem_mat))
  mem_mat <- mem_mat[, seq_len(n_keep), drop = FALSE]
  names(mem_mat) <- paste0("MEM", seq_len(n_keep))
  
  bind_cols(df, mem_mat)
}

make_listw <- function(df, k = 8, coords = c("x_dbmem", "y_dbmem")) {
  
  xy <- df %>%
    dplyr::select(all_of(coords)) %>%
    as.matrix()
  
  if (nrow(xy) < 3) {
    stop("Not enough sites to build spatial weights.")
  }
  
  k_use <- min(k, nrow(xy) - 1)
  
  knn <- spdep::knearneigh(xy, k = k_use)
  nb  <- spdep::knn2nb(knn)
  
  spdep::nb2listw(
    nb,
    style = "W",
    zero.policy = TRUE
  )
}

moran_lm <- function(df, response, predictors, listw) {
  
  f <- make_formula(response, predictors)
  mod <- lm(f, data = df)
  res <- residuals(mod)
  
  out <- tryCatch({
    mt <- spdep::moran.test(
      res,
      listw,
      zero.policy = TRUE,
      alternative = "greater"
    )
    tibble::tibble(
      response = response,
      moran_I = unname(mt$estimate[["Moran I statistic"]]),
      moran_p = mt$p.value,
      n = length(res)
    )
  }, error = function(e) {
    tibble::tibble(
      response = response,
      moran_I = NA_real_,
      moran_p = NA_real_,
      n = length(res)
    )
  })
  
  out
}


fit_richness_gls <- function(df, predictors = c("built_buffer", "landdiv_buffer")) {
  
  f <- make_formula("species_richness", predictors)
  
  # Use scaled coordinates to avoid numerical issues with projected metres.
  if (!all(c("x_gls", "y_gls") %in% names(df))) {
    df <- df %>%
      dplyr::mutate(
        x_gls = as.numeric(scale(x_dbmem)),
        y_gls = as.numeric(scale(y_dbmem))
      )
  }
  
  fit <- tryCatch(
    nlme::gls(
      f,
      data = df,
      correlation = nlme::corExp(
        form = ~ x_gls + y_gls,
        nugget = TRUE
      ),
      method = "ML",
      control = nlme::glsControl(
        msMaxIter = 200,
        msVerbose = FALSE,
        returnObject = TRUE
      )
    ),
    error = function(e1) {
      warning(
        "corExp GLS failed for species_richness: ", e1$message,
        ". Trying corGaus."
      )
      tryCatch(
        nlme::gls(
          f,
          data = df,
          correlation = nlme::corGaus(
            form = ~ x_gls + y_gls,
            nugget = TRUE
          ),
          method = "ML",
          control = nlme::glsControl(
            msMaxIter = 200,
            msVerbose = FALSE,
            returnObject = TRUE
          )
        ),
        error = function(e2) {
          warning(
            "corGaus GLS also failed for species_richness: ", e2$message,
            ". Falling back to lm()."
          )
          stats::lm(f, data = df)
        }
      )
    }
  )
  
  fit
}

moran_model_residuals <- function(mod, listw, response = NA_character_) {
  
  res <- tryCatch(
    {
      if (inherits(mod, "gls")) {
        stats::residuals(mod, type = "normalized")
      } else {
        stats::residuals(mod)
      }
    },
    error = function(e) stats::residuals(mod)
  )
  
  out <- tryCatch({
    mt <- spdep::moran.test(
      res,
      listw,
      zero.policy = TRUE,
      alternative = "greater"
    )
    tibble::tibble(
      response = response,
      moran_I = unname(mt$estimate[["Moran I statistic"]]),
      moran_p = mt$p.value,
      n = length(res)
    )
  }, error = function(e) {
    tibble::tibble(
      response = response,
      moran_I = NA_real_,
      moran_p = NA_real_,
      n = length(res)
    )
  })
  
  out
}

moran_gls_richness <- function(df, predictors, listw) {
  mod <- fit_richness_gls(df, predictors = predictors)
  moran_model_residuals(
    mod = mod,
    listw = listw,
    response = "species_richness"
  )
}

select_mems_for_equation <- function(df, response, predictors, mem_candidates,
                                     listw, alpha = 0.05, max_select = 30) {
  
  selected <- character(0)
  remaining <- mem_candidates
  
  initial <- moran_lm(
    df = df,
    response = response,
    predictors = predictors,
    listw = listw
  )
  
  current <- initial
  
  if (is.na(current$moran_p) || current$moran_p >= alpha) {
    return(
      tibble::tibble(
        response = response,
        selected_mems = "",
        n_mems = 0,
        initial_moran_I = initial$moran_I,
        initial_moran_p = initial$moran_p,
        final_moran_I = current$moran_I,
        final_moran_p = current$moran_p
      )
    )
  }
  
  while (
    current$moran_p < alpha &&
    length(selected) < max_select &&
    length(remaining) > 0
  ) {
    
    candidates_list <- list()
    
    for (mem in remaining) {
      
      preds_test <- c(predictors, selected, mem)
      
      out <- try(
        moran_lm(
          df = df,
          response = response,
          predictors = preds_test,
          listw = listw
        ),
        silent = TRUE
      )
      
      if (!inherits(out, "try-error")) {
        candidates_list[[mem]] <- out %>%
          dplyr::mutate(mem = mem) %>%
          dplyr::select(mem, moran_I, moran_p)
      }
    }
    
    candidates_test <- dplyr::bind_rows(candidates_list) %>%
      dplyr::filter(!is.na(moran_p))
    
    if (nrow(candidates_test) == 0) break
    
    best_mem <- candidates_test %>%
      dplyr::arrange(desc(moran_p), abs(moran_I)) %>%
      dplyr::slice(1) %>%
      dplyr::pull(mem)
    
    selected <- c(selected, best_mem)
    remaining <- setdiff(remaining, best_mem)
    
    current <- moran_lm(
      df = df,
      response = response,
      predictors = c(predictors, selected),
      listw = listw
    )
  }
  
  tibble::tibble(
    response = response,
    selected_mems = paste(selected, collapse = " + "),
    n_mems = length(selected),
    initial_moran_I = initial$moran_I,
    initial_moran_p = initial$moran_p,
    final_moran_I = current$moran_I,
    final_moran_p = current$moran_p
  )
}

# --------------------------------------------------------------------------- #
# 5. Prepare and residualize SEM data by region
# --------------------------------------------------------------------------- #

sem_base_vars <- c(
  "SITE_ID",
  "city",
  "species_asynchrony",
  "wm_population_stability",
  "shannon_diversity",
  "FDis",
  "MPD",
  "species_richness",
  "urban_context",
  "built1000", "built2000", "built5000",
  "landdiv1000", "landdiv2000", "landdiv5000",
  "x", "y"
)

# Residualization logic:
#   species_richness = base diversity axis
#   shannon_resid    = Shannon diversity independent of species richness
#   FDis_resid       = functional diversity independent of richness and residual Shannon
#   MPD_resid        = phylogenetic diversity independent of richness, residual Shannon,
#                      and residual functional diversity
#
# This keeps richness as an explicit diversity predictor of stability, while
# the three residualized diversity variables represent non-richness components.
# In v10, richness is modelled as an endogenous response to landscape structure
# to retain the ecological pathway landscape -> richness -> stability.

residualize_diversity_vars <- function(df) {
  
  df %>%
    mutate(
      shannon_resid = residuals(
        lm(
          shannon_diversity ~ species_richness,
          data = .,
          na.action = na.exclude
        )
      ),
      FDis_resid = residuals(
        lm(
          FDis ~ species_richness + shannon_resid,
          data = .,
          na.action = na.exclude
        )
      ),
      MPD_resid = residuals(
        lm(
          MPD ~ species_richness + shannon_resid + FDis_resid,
          data = .,
          na.action = na.exclude
        )
      )
    ) %>%
    mutate(
      across(
        c(shannon_resid, FDis_resid, MPD_resid),
        ~ as.numeric(scale(.))
      )
    )
}

prepare_region_sem_data <- function(region, sem_raw) {
  
  missing_vars <- setdiff(sem_base_vars, names(sem_raw))
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variables in ", region, " SEM data: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  scale_vars <- setdiff(
    sem_base_vars,
    c("SITE_ID", "city", "urban_context", "x", "y")
  )
  
  sem_scaled <- sem_raw %>%
    dplyr::select(dplyr::all_of(sem_base_vars)) %>%
    mutate(
      city = region,
      urban_context = as.numeric(urban_context),
      across(
        all_of(scale_vars),
        ~ as.numeric(scale(.))
      )
    )
  
  sem_resid <- residualize_diversity_vars(sem_scaled)
  
  resid_vars <- c(
    "species_richness",
    "shannon_resid",
    "FDis_resid",
    "MPD_resid"
  )
  
  cor_resid_vars <- cor(
    sem_resid %>% select(all_of(resid_vars)),
    use = "pairwise.complete.obs"
  )
  
  cor_resid_long <- as.data.frame(cor_resid_vars) %>%
    tibble::rownames_to_column("var1") %>%
    tidyr::pivot_longer(
      cols = -var1,
      names_to = "var2",
      values_to = "correlation"
    ) %>%
    mutate(region = region, .before = 1)
  
  vif_resid_vars <- manual_vif(
    sem_resid %>% drop_na(all_of(resid_vars)),
    resid_vars
  ) %>%
    mutate(region = region, .before = 1)
  
  message("\n================ ", region, " CORRELATIONS AFTER RESIDUALIZATION ================")
  print(cor_resid_vars)
  
  message("\n================ ", region, " VIF AFTER RESIDUALIZATION ================")
  print(vif_resid_vars)
  
  list(
    data = sem_resid,
    correlations = cor_resid_long,
    vif = vif_resid_vars
  )
}


# --------------------------------------------------------------------------- #
# 6. Helper functions: d-separation and ecological Fisher C
# --------------------------------------------------------------------------- #

calc_fisher_from_pvalues <- function(p_values) {
  p_values <- suppressWarnings(as.numeric(p_values))
  p_values <- p_values[!is.na(p_values) & p_values > 0]
  
  if (length(p_values) == 0) {
    return(
      tibble::tibble(
        Fisher.C = NA_real_,
        df = 0,
        P.Value = NA_real_,
        n_claims = 0
      )
    )
  }
  
  fisher_c <- -2 * sum(log(p_values))
  fisher_df <- 2 * length(p_values)
  fisher_p <- stats::pchisq(fisher_c, df = fisher_df, lower.tail = FALSE)
  
  tibble::tibble(
    Fisher.C = fisher_c,
    df = fisher_df,
    P.Value = fisher_p,
    n_claims = length(p_values)
  )
}

extract_dsep_tables <- function(psem_fit, region, buffer) {
  
  dsep_full <- tryCatch(
    piecewiseSEM::dSep(psem_fit) %>%
      clean_empty_names() %>%
      dplyr::mutate(region = region, buffer = buffer, .before = 1),
    error = function(e) {
      warning(
        "dSep failed for ", region, " buffer ", buffer,
        ": ", e$message
      )
      tibble::tibble(
        region = region,
        buffer = buffer,
        dsep_error = e$message
      )
    }
  )
  
  if (!"Independ.Claim" %in% names(dsep_full)) {
    return(
      list(
        dsep_full = dsep_full,
        dsep_ecological = tibble::tibble(region = region, buffer = buffer),
        fisher_ecological = tibble::tibble(
          region = region,
          buffer = buffer,
          Fisher.C = NA_real_,
          df = NA_real_,
          P.Value = NA_real_,
          n_claims = NA_integer_,
          fisher_type = "ecological_dSep_excluding_MEM_claims"
        )
      )
    )
  }
  
  # MEMs are nuisance spatial covariates. For an ecological fit diagnostic,
  # we remove d-separation claims where the tested independence relationship
  # directly involves a MEM variable. The remaining claims still come from the
  # spatially corrected component models, so MEMs remain included as controls
  # when they are part of the relevant component equation.
  dsep_ecological <- dsep_full %>%
    dplyr::filter(
      !stringr::str_detect(Independ.Claim, "\\bMEM[0-9]+\\b")
    ) %>%
    dplyr::mutate(fisher_type = "ecological_dSep_excluding_MEM_claims")
  
  fisher_ecological <- calc_fisher_from_pvalues(dsep_ecological$P.Value) %>%
    dplyr::mutate(
      region = region,
      buffer = buffer,
      fisher_type = "ecological_dSep_excluding_MEM_claims",
      .before = 1
    )
  
  list(
    dsep_full = dsep_full,
    dsep_ecological = dsep_ecological,
    fisher_ecological = fisher_ecological
  )
}

# --------------------------------------------------------------------------- #
# 7. Function to run pSEM for one region and one buffer
# --------------------------------------------------------------------------- #

run_region_psem_buffer <- function(region, buffer, data_resid,
                                   max_mem = 40, k_neigh = 8,
                                   alpha = 0.05, max_select = 30,
                                   small_n_mem_fraction = 0.25,
                                   min_resid_df = 8) {
  
  message("\n============================================================")
  message("Running pSEM for region: ", region, " | buffer: ", buffer, " m")
  message("============================================================")
  
  built_var   <- paste0("built", buffer)
  landdiv_var <- paste0("landdiv", buffer)
  
  message("Using landscape variables: ", built_var, " and ", landdiv_var)
  
  required_vars <- c(
    "SITE_ID",
    "city",
    "species_asynchrony",
    "wm_population_stability",
    "species_richness",
    "shannon_resid",
    "FDis_resid",
    "MPD_resid",
    "urban_context",
    built_var,
    landdiv_var,
    "x",
    "y"
  )
  
  missing_input <- setdiff(required_vars, names(data_resid))
  
  if (length(missing_input) > 0) {
    stop(
      "Missing variables in data_resid for ", region, " buffer ", buffer, ": ",
      paste(missing_input, collapse = ", ")
    )
  }
  
  df_psem_input <- data_resid %>%
    dplyr::select(dplyr::all_of(required_vars)) %>%
    dplyr::mutate(
      built_buffer = .data[[built_var]],
      landdiv_buffer = .data[[landdiv_var]]
    ) %>%
    dplyr::select(-dplyr::all_of(c(built_var, landdiv_var)))
  
  site_diagnostics <- diagnose_psem_complete_cases(
    df_input = df_psem_input,
    region = region,
    buffer = buffer,
    results_dir = results_dir
  )
  
  df_psem <- df_psem_input %>%
    tidyr::drop_na() %>%
    jitter_duplicate_xy(
      x_col = "x",
      y_col = "y",
      amount = 1
    ) %>%
    dplyr::mutate(
      x_gls = as.numeric(scale(x_dbmem)),
      y_gls = as.numeric(scale(y_dbmem))
    )
  
  message("N complete-case sites used in pSEM: ", nrow(df_psem))
  
  if (nrow(df_psem) < 10) {
    warning("Very few sites for ", region, " buffer ", buffer, ". Results may be unstable.")
  }
  
  # Small regions such as BCN have few sites. Allowing too many MEMs can
  # quasi-saturate the component lm() models and make piecewiseSEM fail when
  # it builds the d-separation tests. This cap keeps spatial filtering useful
  # but prevents overfitting. In v10 the richness equation may need more MEMs,
  # especially in RND, so the fraction cap is slightly more permissive than in v9.
  sample_size_cap <- max(
    0,
    floor(nrow(df_psem) * small_n_mem_fraction)
  )
  
  effective_max_select_global <- min(
    max_select,
    sample_size_cap
  )
  
  message(
    "Maximum MEMs allowed per equation after sample-size cap: ",
    effective_max_select_global
  )
  
  df_psem <- add_dbmems(
    df = df_psem,
    max_mem = max_mem,
    coords = c("x_dbmem", "y_dbmem")
  )
  
  mem_candidates <- names(df_psem)[grepl("^MEM", names(df_psem))]
  
  if (length(mem_candidates) == 0) {
    stop("No MEM variables were created for ", region, " buffer ", buffer, ".")
  }
  
  listw <- make_listw(
    df = df_psem,
    k = k_neigh,
    coords = c("x_dbmem", "y_dbmem")
  )
  
  # Ecological pSEM structure.
  # urban_context is kept in the data only for site diagnostics, but it is not
  # included in the pSEM. The urban/rural contrast is tested in the main GLMs.
  # Here, the pSEM focuses on continuous landscape-mediated pathways.
  # In v10, richness is modelled as a response to landscape structure.
  eqs <- list(
    species_asynchrony = c(
      "species_richness",
      "shannon_resid",
      "FDis_resid",
      "MPD_resid",
      "built_buffer",
      "landdiv_buffer"
    ),
    
    wm_population_stability = c(
      "species_richness",
      "shannon_resid",
      "FDis_resid",
      "MPD_resid",
      "built_buffer",
      "landdiv_buffer"
    ),
    
    species_richness = c(
      "built_buffer",
      "landdiv_buffer"
    ),
    
    shannon_resid = c(
      "built_buffer",
      "landdiv_buffer"
    ),
    
    FDis_resid = c(
      "built_buffer",
      "landdiv_buffer"
    ),
    
    MPD_resid = c(
      "built_buffer",
      "landdiv_buffer"
    )
  )
  
  all_model_vars <- unique(c(names(eqs), unlist(eqs)))
  missing_model_vars <- setdiff(all_model_vars, names(df_psem))
  
  if (length(missing_model_vars) > 0) {
    stop(
      "Variables used in SEM but missing from df_psem: ",
      paste(missing_model_vars, collapse = ", ")
    )
  }
  
  max_select_for_response <- function(response_name) {
    
    n_base_terms <- length(eqs[[response_name]])
    
    # lm residual df = n - intercept - base terms - selected MEMs.
    # Keep at least min_resid_df residual df in each component equation.
    df_cap <- nrow(df_psem) - 1 - n_base_terms - min_resid_df
    df_cap <- max(0, df_cap)
    
    response_fraction_cap <- sample_size_cap
    
    min(
      max_select,
      response_fraction_cap,
      df_cap,
      length(mem_candidates)
    )
  }
  
  # ------------------------------------------------------------------------- #
  # Moran before
  # ------------------------------------------------------------------------- #
  
  moran_before_list <- list()
  
  for (resp in names(eqs)) {
    preds <- eqs[[resp]]
    moran_before_list[[resp]] <- moran_lm(
      df = df_psem,
      response = resp,
      predictors = preds,
      listw = listw
    )
  }
  
  moran_before <- dplyr::bind_rows(moran_before_list) %>%
    dplyr::mutate(region = region, buffer = buffer, stage = "before", .before = 1)
  
  # ------------------------------------------------------------------------- #
  # Select MEMs
  # ------------------------------------------------------------------------- #
  
  mem_selection_list <- list()
  
  for (resp in names(eqs)) {
    preds <- eqs[[resp]]
    this_max_select <- max_select_for_response(resp)
    
    if (resp == "species_richness" && use_gls_for_richness) {
      message(
        "Using spatial GLS for ", resp,
        " | no MEM selection for this equation"
      )
      
      initial_richness <- moran_lm(
        df = df_psem,
        response = resp,
        predictors = preds,
        listw = listw
      )
      
      final_richness <- moran_gls_richness(
        df = df_psem,
        predictors = preds,
        listw = listw
      )
      
      mem_selection_list[[resp]] <- tibble::tibble(
        response = resp,
        selected_mems = "",
        n_mems = 0,
        initial_moran_I = initial_richness$moran_I,
        initial_moran_p = initial_richness$moran_p,
        final_moran_I = final_richness$moran_I,
        final_moran_p = final_richness$moran_p,
        spatial_method = "GLS_corExp_or_fallback",
        max_mems_allowed = 0
      )
      next
    }
    
    message(
      "Selecting MEMs for ", resp,
      " | max allowed = ", this_max_select
    )
    
    mem_selection_list[[resp]] <- select_mems_for_equation(
      df = df_psem,
      response = resp,
      predictors = preds,
      mem_candidates = mem_candidates,
      listw = listw,
      alpha = alpha,
      max_select = this_max_select
    ) %>%
      dplyr::mutate(
        spatial_method = "MEM",
        max_mems_allowed = this_max_select
      )
  }
  
  mem_selection <- dplyr::bind_rows(mem_selection_list) %>%
    dplyr::mutate(region = region, buffer = buffer, .before = 1)
  
  # ------------------------------------------------------------------------- #
  # Build spatial equations
  # ------------------------------------------------------------------------- #
  
  spatial_terms_list <- vector("list", nrow(mem_selection))
  names(spatial_terms_list) <- mem_selection$response
  
  for (i in seq_len(nrow(mem_selection))) {
    spatial_terms_list[[mem_selection$response[i]]] <- split_mem_string(
      mem_selection$selected_mems[i]
    )
  }
  
  eqs_spatial <- list()
  
  for (resp in names(eqs)) {
    preds <- eqs[[resp]]
    extra <- spatial_terms_list[[resp]]
    if (is.null(extra)) extra <- character(0)
    eqs_spatial[[resp]] <- c(preds, extra)
  }
  
  missing_spatial_vars <- lapply(
    eqs_spatial,
    function(x) setdiff(x, names(df_psem))
  )
  
  if (any(lengths(missing_spatial_vars) > 0)) {
    print(missing_spatial_vars)
    stop("Some selected spatial terms are missing from df_psem.")
  }
  
  # ------------------------------------------------------------------------- #
  # Fit component models
  # ------------------------------------------------------------------------- #
  
  model_env <- environment()
  
  m_async <- lm(
    make_formula("species_asynchrony", eqs_spatial$species_asynchrony, env = model_env),
    data = df_psem
  )
  
  m_pop <- lm(
    make_formula("wm_population_stability", eqs_spatial$wm_population_stability, env = model_env),
    data = df_psem
  )
  
  if (use_gls_for_richness) {
    m_richness <- fit_richness_gls(
      df = df_psem,
      predictors = eqs$species_richness
    )
    
    # IMPORTANT: fit_richness_gls() is a helper function, so the original
    # gls call stores local helper objects (e.g. data = df, model = f).
    # piecewiseSEM later re-evaluates the model call to recover model data.
    # A local object such as `df_psem` is not visible when pSEM/dSep/coefs
    # later re-evaluate the gls call, so we store a unique copy in .GlobalEnv
    # and patch the gls call to point to that stable object.
    gls_data_name <- paste0(".psem_gls_data_", region, "_", buffer)
    assign(gls_data_name, df_psem, envir = .GlobalEnv)
    
    m_richness$call$model <- make_formula(
      "species_richness",
      eqs$species_richness,
      env = .GlobalEnv
    )
    m_richness$call$data <- as.name(gls_data_name)
  } else {
    m_richness <- lm(
      make_formula("species_richness", eqs_spatial$species_richness, env = model_env),
      data = df_psem
    )
  }
  
  m_shannon <- lm(
    make_formula("shannon_resid", eqs_spatial$shannon_resid, env = model_env),
    data = df_psem
  )
  
  m_fdis <- lm(
    make_formula("FDis_resid", eqs_spatial$FDis_resid, env = model_env),
    data = df_psem
  )
  
  m_mpd <- lm(
    make_formula("MPD_resid", eqs_spatial$MPD_resid, env = model_env),
    data = df_psem
  )
  
  # ------------------------------------------------------------------------- #
  # Build pSEM
  # ------------------------------------------------------------------------- #
  
  psem_fit <- piecewiseSEM::psem(
    m_async,
    m_pop,
    m_richness,
    m_shannon,
    m_fdis,
    m_mpd,
    species_asynchrony %~~% wm_population_stability,
    built_buffer %~~% landdiv_buffer,
    data = df_psem
  )
  
  # ------------------------------------------------------------------------- #
  # Moran after
  # ------------------------------------------------------------------------- #
  
  moran_after_list <- list()
  
  for (resp in names(eqs)) {
    preds <- eqs[[resp]]
    extra <- spatial_terms_list[[resp]]
    if (is.null(extra)) extra <- character(0)
    
    if (resp == "species_richness" && use_gls_for_richness) {
      moran_after_list[[resp]] <- moran_model_residuals(
        mod = m_richness,
        listw = listw,
        response = resp
      )
    } else {
      moran_after_list[[resp]] <- moran_lm(
        df = df_psem,
        response = resp,
        predictors = c(preds, extra),
        listw = listw
      )
    }
  }
  
  moran_after <- dplyr::bind_rows(moran_after_list) %>%
    dplyr::mutate(region = region, buffer = buffer, stage = "after", .before = 1)
  
  moran_all <- dplyr::bind_rows(moran_before, moran_after)
  
  # ------------------------------------------------------------------------- #
  # Extract outputs
  # ------------------------------------------------------------------------- #
  
  psem_summary <- tryCatch(
    summary(psem_fit),
    error = function(e) {
      warning(
        "pSEM summary failed for ", region, " buffer ", buffer,
        ": ", e$message
      )
      NULL
    }
  )
  
  component_models <- list(
    species_asynchrony = m_async,
    wm_population_stability = m_pop,
    species_richness = m_richness,
    shannon_resid = m_shannon,
    FDis_resid = m_fdis,
    MPD_resid = m_mpd
  )
  
  # Use manual extraction from the component lm() models for the path table.
  # This is more robust for small regional SEMs such as BCN, where
  # piecewiseSEM::coefs() can fail while the component models are valid.
  psem_coefs <- extract_component_path_coefs(
    component_models = component_models,
    model_data = df_psem,
    region = region,
    buffer = buffer,
    built_var = built_var,
    landdiv_var = landdiv_var
  )
  
  psem_coefs_clean <- psem_coefs %>%
    dplyr::filter(!grepl("^MEM", Predictor))
  
  psem_coefs_piecewise <- tryCatch(
    piecewiseSEM::coefs(
      psem_fit,
      standardize = "scale"
    ) %>%
      clean_empty_names() %>%
      dplyr::mutate(
        region = region,
        buffer = buffer,
        built_var = built_var,
        landdiv_var = landdiv_var,
        extraction = "piecewiseSEM",
        .before = 1
      ),
    error = function(e) {
      warning(
        "pSEM coefficients failed for ", region, " buffer ", buffer,
        ": ", e$message,
        ". Manual lm() coefficients were still extracted."
      )
      tibble::tibble(
        region = region,
        buffer = buffer,
        error = e$message,
        extraction = "piecewiseSEM_failed"
      )
    }
  )
  
  psem_r2 <- tryCatch(
    piecewiseSEM::rsquared(psem_fit) %>%
      clean_empty_names() %>%
      dplyr::mutate(region = region, buffer = buffer, .before = 1),
    error = function(e) {
      warning(
        "pSEM R2 failed for ", region, " buffer ", buffer,
        ": ", e$message
      )
      tibble::tibble(region = region, buffer = buffer, error = e$message)
    }
  )
  
  psem_fisher <- tryCatch(
    piecewiseSEM::fisherC(psem_fit) %>%
      clean_empty_names() %>%
      dplyr::mutate(region = region, buffer = buffer, .before = 1),
    error = function(e) {
      warning(
        "Fisher's C failed for ", region, " buffer ", buffer,
        ": ", e$message
      )
      tibble::tibble(
        region = region,
        buffer = buffer,
        Fisher.C = NA_real_,
        df = NA_real_,
        P.Value = NA_real_,
        fisher_error = e$message
      )
    }
  )
  
  dsep_outputs <- extract_dsep_tables(
    psem_fit = psem_fit,
    region = region,
    buffer = buffer
  )
  
  psem_dsep_full <- dsep_outputs$dsep_full
  psem_dsep_ecological <- dsep_outputs$dsep_ecological
  psem_fisher_ecological <- dsep_outputs$fisher_ecological
  
  psem_aic <- tryCatch(
    AIC(psem_fit) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(region = region, buffer = buffer, .before = 1),
    error = function(e) {
      warning(
        "pSEM AIC failed for ", region, " buffer ", buffer,
        ": ", e$message
      )
      tibble::tibble(region = region, buffer = buffer, error = e$message)
    }
  )
  
  list(
    region = region,
    buffer = buffer,
    built_var = built_var,
    landdiv_var = landdiv_var,
    data = df_psem,
    eqs = eqs,
    eqs_spatial = eqs_spatial,
    mem_selection = mem_selection,
    moran_all = moran_all,
    site_retention = site_diagnostics$site_retention,
    missing_by_variable = site_diagnostics$missing_by_variable,
    dropped_sites = site_diagnostics$dropped_sites,
    psem = psem_fit,
    summary = psem_summary,
    coefs = psem_coefs,
    coefs_clean = psem_coefs_clean,
    coefs_piecewise = psem_coefs_piecewise,
    r2 = psem_r2,
    fisher = psem_fisher,
    fisher_ecological = psem_fisher_ecological,
    dsep_full = psem_dsep_full,
    dsep_ecological = psem_dsep_ecological,
    aic = psem_aic
  )
}

# --------------------------------------------------------------------------- #
# 7. Load data for all regions
# --------------------------------------------------------------------------- #

message("\n================ LOAD COORDINATES ================")
site_coords <- load_site_coordinates(input_dir)

message("\n================ LOAD REGION DATA ================")
sem_raw <- purrr::set_names(regions) %>%
  purrr::map(
    function(r) {
      load_region_sem_raw(
        region = r,
        data_dir = data_dir,
        site_coords = site_coords
      )
    }
  )

all_regions_check <- purrr::imap_dfr(
  sem_raw,
  function(df, r) {
    df %>%
      mutate(region = r) %>%
      count(region, urban_context)
  }
)

message("\n==== Final region urban/rural counts ====")
print(all_regions_check)

# --------------------------------------------------------------------------- #
# 8. Standardize and residualize by region
# --------------------------------------------------------------------------- #

sem_prepared <- purrr::imap(
  sem_raw,
  function(df, r) {
    prepare_region_sem_data(
      region = r,
      sem_raw = df
    )
  }
)

sem_resid <- purrr::map(sem_prepared, "data")

resid_correlations_all <- purrr::map_dfr(sem_prepared, "correlations")
resid_vif_all <- purrr::map_dfr(sem_prepared, "vif")

write.csv(
  resid_correlations_all,
  file.path(results_dir, "AllRegions_residualized_diversity_correlations.csv"),
  row.names = FALSE
)

write.csv(
  resid_vif_all,
  file.path(results_dir, "AllRegions_residualized_diversity_VIF.csv"),
  row.names = FALSE
)

# --------------------------------------------------------------------------- #
# 9. Run pSEMs for all regions and buffers
# --------------------------------------------------------------------------- #

psem_results <- list()

for (r in regions) {
  
  psem_results[[r]] <- list()
  
  for (b in buffers) {
    
    message("\n\n########## START REGION ", r, " | BUFFER ", b, " ##########")
    
    psem_results[[r]][[as.character(b)]] <- run_region_psem_buffer(
      region = r,
      buffer = b,
      data_resid = sem_resid[[r]],
      max_mem = max_mem,
      k_neigh = k_neigh,
      alpha = alpha_moran,
      max_select = max_select
    )
    
    message("########## END REGION ", r, " | BUFFER ", b, " ##########\n")
  }
}

# --------------------------------------------------------------------------- #
# 10. Combine outputs
# --------------------------------------------------------------------------- #

flatten_results <- function(results_nested, element) {
  purrr::map_dfr(
    results_nested,
    function(region_list) {
      purrr::map_dfr(region_list, element)
    }
  )
}

psem_fisher_all <- flatten_results(psem_results, "fisher")
psem_fisher_ecological_all <- flatten_results(psem_results, "fisher_ecological")
psem_dsep_full_all <- flatten_results(psem_results, "dsep_full")
psem_dsep_ecological_all <- flatten_results(psem_results, "dsep_ecological")
psem_aic_all <- flatten_results(psem_results, "aic")
psem_r2_all <- flatten_results(psem_results, "r2")
psem_coefs_all <- flatten_results(psem_results, "coefs")
psem_coefs_clean_all <- flatten_results(psem_results, "coefs_clean")
psem_coefs_piecewise_all <- flatten_results(psem_results, "coefs_piecewise")
psem_moran_all <- flatten_results(psem_results, "moran_all")
psem_mem_selection_all <- flatten_results(psem_results, "mem_selection")
psem_site_retention_all <- flatten_results(psem_results, "site_retention")
psem_missing_by_variable_all <- flatten_results(psem_results, "missing_by_variable")
psem_dropped_sites_all <- flatten_results(psem_results, "dropped_sites")

# --------------------------------------------------------------------------- #
# 11. Spatial correction check
# --------------------------------------------------------------------------- #

psem_moran_check <- psem_moran_all %>%
  filter(stage == "after") %>%
  mutate(moran_corrected = moran_p >= 0.05) %>%
  arrange(region, buffer, moran_p)

# --------------------------------------------------------------------------- #
# 12. Coefficient tables across regions and buffers
# --------------------------------------------------------------------------- #

psem_dsep_ecological_problematic <- psem_dsep_ecological_all %>%
  dplyr::filter(
    "P.Value" %in% names(.),
    !is.na(P.Value),
    P.Value < 0.05
  ) %>%
  dplyr::arrange(region, buffer, P.Value)

psem_key_paths <- psem_coefs_clean_all %>%
  filter(P.Value < 0.1) %>%
  dplyr::select(
    region,
    buffer,
    Response,
    Predictor,
    Estimate,
    Std.Estimate,
    P.Value,
    signif
  ) %>%
  arrange(region, Response, Predictor, buffer)

psem_paths_wide <- psem_coefs_clean_all %>%
  mutate(
    path = paste(Response, "<-", Predictor),
    estimate_label = paste0(
      round(Std.Estimate, 3),
      " ",
      .data$signif,
      " (p = ",
      base::signif(P.Value, 3),
      ")"
    )
  ) %>%
  dplyr::select(region, path, buffer, estimate_label) %>%
  pivot_wider(
    names_from = buffer,
    values_from = estimate_label,
    names_prefix = "buffer_"
  ) %>%
  arrange(region, path)

main_paths <- c(
  "species_richness <- built_buffer",
  "species_richness <- landdiv_buffer",
  "shannon_resid <- built_buffer",
  "shannon_resid <- landdiv_buffer",
  "FDis_resid <- built_buffer",
  "FDis_resid <- landdiv_buffer",
  "MPD_resid <- built_buffer",
  "MPD_resid <- landdiv_buffer",
  "species_asynchrony <- species_richness",
  "species_asynchrony <- shannon_resid",
  "species_asynchrony <- FDis_resid",
  "species_asynchrony <- MPD_resid",
  "species_asynchrony <- built_buffer",
  "species_asynchrony <- landdiv_buffer",
  "wm_population_stability <- species_richness",
  "wm_population_stability <- shannon_resid",
  "wm_population_stability <- FDis_resid",
  "wm_population_stability <- MPD_resid",
  "wm_population_stability <- built_buffer",
  "wm_population_stability <- landdiv_buffer",
  "~~species_asynchrony <- ~~wm_population_stability",
  "~~built_buffer <- ~~landdiv_buffer"
)

psem_main_paths_wide <- psem_paths_wide %>%
  filter(path %in% main_paths)

# --------------------------------------------------------------------------- #
# 13. Save outputs
# --------------------------------------------------------------------------- #

write.csv(
  psem_fisher_all,
  file.path(results_dir, "AllRegions_pSEM_fisher_all_buffers_full_with_MEMs.csv"),
  row.names = FALSE
)

write.csv(
  psem_fisher_ecological_all,
  file.path(results_dir, "AllRegions_pSEM_fisher_ecological_dSep_excluding_MEM_claims.csv"),
  row.names = FALSE
)

write.csv(
  psem_dsep_full_all,
  file.path(results_dir, "AllRegions_pSEM_dSep_all_buffers_full_with_MEMs.csv"),
  row.names = FALSE
)

write.csv(
  psem_dsep_ecological_all,
  file.path(results_dir, "AllRegions_pSEM_dSep_ecological_excluding_MEM_claims.csv"),
  row.names = FALSE
)

write.csv(
  psem_dsep_ecological_problematic,
  file.path(results_dir, "AllRegions_pSEM_dSep_ecological_problematic_p_lt_0.05.csv"),
  row.names = FALSE
)

write.csv(
  psem_aic_all,
  file.path(results_dir, "AllRegions_pSEM_AIC_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_r2_all,
  file.path(results_dir, "AllRegions_pSEM_R2_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_coefs_all,
  file.path(results_dir, "AllRegions_pSEM_coefficients_all_buffers_with_MEMs.csv"),
  row.names = FALSE
)

write.csv(
  psem_coefs_clean_all,
  file.path(results_dir, "AllRegions_pSEM_coefficients_all_buffers_clean.csv"),
  row.names = FALSE
)

write.csv(
  psem_coefs_piecewise_all,
  file.path(results_dir, "AllRegions_pSEM_coefficients_piecewiseSEM_raw.csv"),
  row.names = FALSE
)

write.csv(
  psem_key_paths,
  file.path(results_dir, "AllRegions_pSEM_key_paths_p_lt_0.1.csv"),
  row.names = FALSE
)

write.csv(
  psem_paths_wide,
  file.path(results_dir, "AllRegions_pSEM_paths_wide_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_main_paths_wide,
  file.path(results_dir, "AllRegions_pSEM_main_paths_wide_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_moran_all,
  file.path(results_dir, "AllRegions_pSEM_Moran_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_moran_check,
  file.path(results_dir, "AllRegions_pSEM_Moran_after_check.csv"),
  row.names = FALSE
)

write.csv(
  psem_mem_selection_all,
  file.path(results_dir, "AllRegions_pSEM_MEM_selection_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_site_retention_all,
  file.path(results_dir, "AllRegions_pSEM_site_retention_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_missing_by_variable_all,
  file.path(results_dir, "AllRegions_pSEM_missing_by_variable_all_buffers.csv"),
  row.names = FALSE
)

write.csv(
  psem_dropped_sites_all,
  file.path(results_dir, "AllRegions_pSEM_dropped_sites_all_buffers.csv"),
  row.names = FALSE
)

saveRDS(
  psem_results,
  file.path(results_dir, "AllRegions_pSEM_full_results_all_buffers.rds")
)

# --------------------------------------------------------------------------- #
# 14. Print key outputs
# --------------------------------------------------------------------------- #

message("\n================ MODEL FIT: FULL FISHER C WITH MEMs ================")
print(psem_fisher_all, n = 100)

message("\n================ MODEL FIT: ECOLOGICAL FISHER C EXCLUDING MEM CLAIMS ================")
print(psem_fisher_ecological_all, n = 100)

message("\n================ ECOLOGICAL dSEP CLAIMS p < 0.05 ================")
print(psem_dsep_ecological_problematic, n = 300)

message("\n================ AIC ================")
print(psem_aic_all, n = 100)

message("\n================ R2 ================")
print(psem_r2_all, n = 300)

message("\n================ MORAN AFTER CHECK ================")
print(psem_moran_check, n = 300)

message("\n================ KEY PATHS p < 0.1 ================")
print(psem_key_paths, n = 300)

message("\n================ MAIN PATHS ACROSS BUFFERS ================")
print(psem_main_paths_wide, n = 300)

message("\n==== All-region pSEM analysis complete ====")
message("Outputs saved in: ", results_dir)
