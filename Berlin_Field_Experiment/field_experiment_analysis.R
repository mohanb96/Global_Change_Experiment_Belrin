
############################################################
# Field experiment analysis for The manuscript Figure 2
# and field-experiment supplementary figures/tables.
#
# Reorganized from Code_field_study.txt.
# Main data: Field_experiment/data_5years.csv
# Main figure produced: Figure 2
# Supplementary outputs produced: raw ridge plots, NMDS/PERMANOVA,
# null-model deviation/SSD plots, and random-forest R2 summaries.
############################################################

########################
# 0. Packages/config   #
########################

required_pkgs <- c(
  "dplyr", "ggplot2", "ggridges", "ggepi", "patchwork",
  "cowplot", "vegan", "caret", "party", "grid", "splines"
)

check_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Please install missing packages first: install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "),
      "))"
    )
  }
}

check_required_packages(required_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(ggepi)
  library(patchwork)
  library(cowplot)
  library(vegan)
  library(caret)
  library(party)
  library(grid)
  library(splines)
})

set.seed(100)

cfg <- list(
  data_path = Sys.getenv(
    "FIELD_DATA_PATH",
    unset = "data_5years.csv"
  ),
  output_dir = Sys.getenv(
    "FIELD_OUTPUT_DIR",
    unset = "Field_experiment/reorganized_outputs"
  ),
  drop_first_column = TRUE,
  save_outputs = TRUE,
  run_main_figure = TRUE,
  run_supp_raw_ridge = TRUE,
  run_supp_nmds = TRUE,
  run_supp_deviation = TRUE,
  run_supp_random_forest = TRUE,
  run_supp_spline_factor_response = TRUE,
  run_lm_factor_number_comparison = TRUE,
  n_boot_mean = 1000,
  n_boot_effect = 1000,
  n_null_boot_expected = 1000,
  # Original field code used 100 bootstrap draws per treatment combination
  # inside the expected-effect null model, even when n_iter was 1000.
  # Keeping this at 100 makes Figure 2 calculations and heatmap RNG order
  # closer to the original script.
  n_expected_boot_per_combo_main = 100,
  n_null_boot_deviation = 100,
  n_boot_slope = 1000,
  # Random-forest settings for Supplementary Fig. S9.
  # Use n_rf_iter = 1000 for final manuscript output; smaller values are only for testing.
  n_rf_iter = 1000,
  n_rf_trees = 50,
  n_eachlv_rf = 12,
  n_rf_expected_boot_per_combo = 100,
  rf_seed = 20260507,
  rf_include_pretreatment = FALSE,
  rf_years = NULL,
  nmds_trymax = 500,
  nmds_maxit = 999,
  permanova_permutations = 999,
  heatmap_p_adjust_method = "BH",
  heatmap_p_two_sided = TRUE,
  # TRUE restores the original random-number order: expected null effects are
  # computed before heatmap treatment-vs-control bootstraps.
  heatmap_legacy_rng_order = TRUE,
  spline_df = 3,
  spline_include_pretreatment = TRUE,
  lm_comparison_include_pretreatment = FALSE,
  lm_comparison_expected_boot_per_combo = 100,
  # Keep FALSE to reproduce the existing supplementary interaction classification.
  # Set TRUE for a conventional two-sided bootstrap tail probability.
  interaction_p_two_sided = FALSE
)

if (cfg$save_outputs) {
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
}

########################
# 1. Study definitions #
########################

############################################################
# Response-control panel
#
# Edit this table to analyze any numeric response variable in
# data_5years.csv. The rest of the script uses this table.
#
# null_model options: "additive", "multiplicative", "dominative"
# direction options : "negative" or "positive"
#   negative = lower response values indicate stronger stress effects
#   positive = higher response values indicate stronger stress effects
############################################################

response_settings <- data.frame(
  response = c("richness_all", "biomass", "shannon_all", "evenness_all"),
  label = c("Plant richness", "Plant biomass", "Shannon diversity", "Plant evenness"),
  null_model = c("multiplicative", "dominative", "multiplicative", "dominative"),
  direction = c("negative", "negative", "negative", "negative"),
  stringsAsFactors = FALSE
)

#response_settings <- data.frame(
#  response = c("biomass_jun", "biomass"),
#  label = c("Biomass June", "Biomass June+October"),
#  null_model = c("dominative", "dominative"),
#  direction = c("negative", "negative"),
#  stringsAsFactors = FALSE
#)

responses <- response_settings$response

response_labels <- stats::setNames(response_settings$label, response_settings$response)

response_label <- function(response_name) {
  lab <- response_labels[[response_name]]
  if (is.null(lab) || is.na(lab) || !nzchar(lab)) {
    return(gsub("_", " ", response_name))
  }
  lab
}

valid_null_models <- c("additive", "multiplicative", "dominative")

null_model_to_E <- function(null_model) {
  null_model <- tolower(null_model)
  out <- ifelse(
    null_model == "additive", "E1",
    ifelse(null_model == "multiplicative", "E2",
           ifelse(null_model == "dominative", "E3", NA_character_))
  )
  out
}

final_null_model_by_response <- stats::setNames(
  tolower(response_settings$null_model),
  response_settings$response
)

if (any(!final_null_model_by_response %in% valid_null_models)) {
  stop("Invalid null_model in response_settings. Use one of: ", paste(valid_null_models, collapse = ", "))
}

expected_effect_column_by_response <- stats::setNames(
  null_model_to_E(final_null_model_by_response),
  names(final_null_model_by_response)
)

response_direction <- stats::setNames(
  tolower(response_settings$direction),
  response_settings$response
)

if (any(!response_direction %in% c("negative", "positive"))) {
  stop("Invalid direction in response_settings. Use 'negative' or 'positive'.")
}

stressors <- c("W", "N", "D", "HM", "MP", "S", "A", "I", "AF", "SF")

stressor_labels <- c(
  W = "Warming",
  N = "N deposition",
  D = "Drought",
  HM = "Heavy metal",
  MP = "Microplastic",
  S = "Salinity",
  A = "Antibiotic",
  I = "Insecticide",
  AF = "Antifungal agent",
  SF = "Surfactant"
)

year_lookup <- c("0" = "2021", "1" = "2022", "2" = "2023", "3" = "2024", "4" = "2025")
lv_levels <- c("1", "2", "4", "6", "8", "10")
h_levels <- lv_levels[-1]
remark_order <- c("CT", stressors, lv_levels)
treat_order_16 <- c(stressors, lv_levels)

lv_cols <- c(
  "0" = "#6D6E71",
  "1" = "#FDE333",
  "2" = "#7ED357",
  "4" = "#00B28A",
  "6" = "#008298",
  "8" = "#274983",
  "10" = "#4B0055"
)

remark_cols <- c(
  "CT" = "#6D6E71",
  "W" = "#726658",
  "N" = "#726658",
  "D" = "#726658",
  "HM" = "#726658",
  "MP" = "#726658",
  "S" = "#726658",
  "A" = "#726658",
  "I" = "#726658",
  "AF" = "#726658",
  "SF" = "#726658",
  "1" = "#FDE333",
  "2" = "#7ED357",
  "4" = "#00B28A",
  "6" = "#008298",
  "8" = "#274983",
  "10" = "#4B0055"
)

########################
# 2. General utilities #
########################

assert_columns <- function(data, cols, context = "data") {
  missing <- setdiff(cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing columns in ", context, ": ", paste(missing, collapse = ", "))
  }
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

read_field_data <- function(data_path, drop_first_column = TRUE) {
  dat <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (drop_first_column) dat <- dat[, -1, drop = FALSE]

  assert_columns(dat, unique(c("time", "Lv", "remark", stressors, responses)), "field data")

  dat$time <- as.character(dat$time)
  dat$Lv <- safe_numeric(dat$Lv)
  dat$remark <- as.character(dat$remark)
  dat$year <- unname(year_lookup[as.character(dat$time)])

  for (s in stressors) dat[[s]] <- safe_numeric(dat[[s]])
  for (r in responses) dat[[r]] <- safe_numeric(dat[[r]])

  dat
}

make_year_list <- function(dat) {
  dat <- dat[!is.na(dat$year), , drop = FALSE]
  out <- split(dat, dat$year)
  out[sort(names(out))]
}

# Single-factor rows are kept twice: once with their identity (W, N, ...),
# and once as factor-number level "1". This allows heatmaps to show both
# individual GCF identity and the one-factor level.
duplicate_single_factor_rows <- function(dataset) {
  lv_1 <- dataset[dataset$Lv == 1 & !is.na(dataset$Lv), , drop = FALSE]
  if (nrow(lv_1) > 0) {
    lv_1$remark <- "1"
    dataset <- rbind(dataset, lv_1)
  }
  dataset
}

prep_year_df <- function(df_tmp) {
  df <- duplicate_single_factor_rows(df_tmp)
  df$.__rowid <- seq_len(nrow(df))
  df$Lv <- safe_numeric(df$Lv)
  df$remark <- as.character(df$remark)
  present <- remark_order[remark_order %in% unique(df$remark)]
  if (length(present) == 0) present <- unique(df$remark)
  df$remark <- factor(df$remark, levels = present)
  for (s in stressors) df[[s]] <- safe_numeric(df[[s]])
  for (r in responses) df[[r]] <- safe_numeric(df[[r]])
  df
}

make_stressor_key <- function(comb_df) {
  if (nrow(comb_df) == 0) return(character(0))
  apply(comb_df, 1, function(x) paste0(as.integer(safe_numeric(x)), collapse = ""))
}

unique_stressor_patterns <- function(df) {
  out <- df[, stressors, drop = FALSE]
  out[] <- lapply(out, safe_numeric)
  out <- unique(out)
  rownames(out) <- NULL
  out
}

safe_tail_p <- function(diff_vec, two_sided = FALSE) {
  diff_vec <- diff_vec[is.finite(diff_vec)]
  if (length(diff_vec) == 0) return(NA_real_)
  p <- mean(diff_vec > 0)
  p <- min(p, 1 - p)
  if (two_sided) p <- min(1, 2 * p)
  p
}

sig_stars <- function(p) {
  ifelse(
    !is.finite(p), "",
    ifelse(p <= 0.001, "***", ifelse(p <= 0.01, "**", ifelse(p <= 0.05, "*", "")))
  )
}

##############################
# 3. Bootstrap summaries     #
##############################

bootstrap_mean_ci <- function(x, n_boot = cfg$n_boot_mean, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) < 1) return(c(ci_low = NA_real_, mean = NA_real_, ci_high = NA_real_))
  bs <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  c(
    ci_low = unname(quantile(bs, probs[1], na.rm = TRUE)),
    mean = mean(bs, na.rm = TRUE),
    ci_high = unname(quantile(bs, probs[2], na.rm = TRUE))
  )
}

bootstrap_group_means <- function(response, data, targets = remark_order, n_boot = cfg$n_boot_mean) {
  out <- vector("list", length(targets))
  names(out) <- targets
  for (trt in targets) {
    population <- data[data$remark == trt, response]
    out[[trt]] <- bootstrap_mean_ci(population, n_boot = n_boot)
  }
  out <- as.data.frame(do.call(rbind, out))
  out$target <- targets
  out <- out[, c("target", "ci_low", "mean", "ci_high")]
  colnames(out) <- c("target", "X2.5.", "mean", "X97.5.")
  out
}

bootstrap_effect_vs_control <- function(response, data, targets = treat_order_16, n_boot = cfg$n_boot_effect) {
  population_CT <- data[data$remark == "CT", response]
  population_CT <- population_CT[is.finite(population_CT)]
  size_CT <- length(population_CT)

  out <- data.frame(
    treatment = targets,
    ES_mean = NA_real_,
    ES_2.5 = NA_real_,
    ES_97.5 = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (size_CT < 2) return(out)

  for (k in seq_along(targets)) {
    trt <- targets[k]

    if (trt == "1") {
      population_TR <- data[data$remark %in% stressors, response]
    } else {
      population_TR <- data[data$remark == trt, response]
    }

    population_TR <- population_TR[is.finite(population_TR)]
    size_TR <- length(population_TR)
    if (size_TR < 2) next

    bs <- replicate(n_boot, {
      k_CT <- mean(sample(population_CT, size_CT, replace = TRUE), na.rm = TRUE)
      k_TR <- mean(sample(population_TR, size_TR, replace = TRUE), na.rm = TRUE)
      k_TR - k_CT
    })

    out$ES_mean[k] <- mean(bs, na.rm = TRUE)
    out$ES_2.5[k] <- unname(quantile(bs, 0.025, na.rm = TRUE))
    out$ES_97.5[k] <- unname(quantile(bs, 0.975, na.rm = TRUE))
    out$p_value[k] <- safe_tail_p(bs, two_sided = cfg$heatmap_p_two_sided)
  }

  out
}

######################################
# 4. Null models and factor-number   #
######################################

simulate_null_effects <- function(response, data, selected_factors, n_boot = 100, return_response = FALSE) {
  null_models <- c("additive", "multiplicative", "dominative")
  output <- setNames(vector("list", length(null_models)), null_models)

  population_CT <- data[data$remark == "CT", response]
  population_CT <- population_CT[is.finite(population_CT)]
  size_CT <- length(population_CT)
  CT <- mean(population_CT, na.rm = TRUE)

  selected_factors <- selected_factors[selected_factors %in% as.character(unique(data$remark))]

  if (size_CT < 2 || !is.finite(CT) || length(selected_factors) == 0) {
    return(lapply(output, function(x) numeric(0)))
  }

  for (type in null_models) {
    bs <- numeric(0)
    for (id in seq_len(n_boot)) {
      k_CT <- mean(sample(population_CT, size_CT, replace = TRUE), na.rm = TRUE)
      if (!is.finite(k_CT)) next

      each_effect <- numeric(0)
      for (trt in selected_factors) {
        population_TR <- data[data$remark == trt, response]
        population_TR <- population_TR[is.finite(population_TR)]
        size_TR <- length(population_TR)
        if (size_TR < 1) next

        k_TR <- mean(sample(population_TR, size_TR, replace = TRUE), na.rm = TRUE)
        if (!is.finite(k_TR)) next

        if (type == "additive") each_effect <- c(each_effect, k_TR - k_CT)
        if (type == "multiplicative") {
          if (abs(k_CT) < .Machine$double.eps) next
          each_effect <- c(each_effect, (k_TR - k_CT) / k_CT)
        }
        if (type == "dominative") each_effect <- c(each_effect, k_TR - k_CT)
      }

      if (length(each_effect) == 0) next

      if (type == "additive") joint_effect <- sum(each_effect)
      if (type == "multiplicative") joint_effect <- (prod(1 + each_effect) - 1) * k_CT
      if (type == "dominative") joint_effect <- each_effect[which.max(abs(each_effect))]

      value <- if (return_response) CT + joint_effect else joint_effect
      if (is.finite(value)) bs <- c(bs, value)
    }
    output[[type]] <- bs
  }

  output
}

calc_expected_effects_for_rows <- function(response, data, n_boot = cfg$n_null_boot_expected) {
  df_lv <- data[data$remark %in% lv_levels, , drop = FALSE]
  if (nrow(df_lv) == 0) {
    return(data.frame(.__rowid = integer(0), E1 = numeric(0), E2 = numeric(0), E3 = numeric(0)))
  }

  effect_map <- list()

  for (Lv in lv_levels) {
    if (Lv == "1") {
      comb <- unique_stressor_patterns(data[data$remark %in% stressors, , drop = FALSE])
    } else {
      comb <- unique_stressor_patterns(data[data$remark == Lv, , drop = FALSE])
    }
    if (nrow(comb) == 0) next

    keys <- make_stressor_key(comb)
    effect_map[[Lv]] <- list()

    for (j in seq_len(nrow(comb))) {
      selected <- stressors[which(safe_numeric(comb[j, stressors]) == 1)]
      sim <- simulate_null_effects(
        response = response,
        data = data,
        selected_factors = selected,
        n_boot = n_boot,
        return_response = FALSE
      )
      effect_map[[Lv]][[keys[j]]] <- c(
        E1 = mean(sim$additive, na.rm = TRUE),
        E2 = mean(sim$multiplicative, na.rm = TRUE),
        E3 = mean(sim$dominative, na.rm = TRUE)
      )
    }
  }

  row_keys <- make_stressor_key(df_lv[, stressors, drop = FALSE])
  row_lv <- as.character(df_lv$Lv)
  row_lv[!(row_lv %in% lv_levels)] <- as.character(df_lv$remark[!(row_lv %in% lv_levels)])

  E1 <- E2 <- E3 <- rep(NA_real_, nrow(df_lv))

  for (i in seq_len(nrow(df_lv))) {
    Lv <- row_lv[i]
    key <- row_keys[i]
    if (!Lv %in% names(effect_map)) next
    if (!key %in% names(effect_map[[Lv]])) next
    vals <- effect_map[[Lv]][[key]]
    E1[i] <- vals["E1"]
    E2[i] <- vals["E2"]
    E3[i] <- vals["E3"]
  }

  data.frame(.__rowid = df_lv$.__rowid, E1 = E1, E2 = E2, E3 = E3)
}

safe_relative_deviation <- function(actual_mean, null_mean) {
  if (!is.finite(actual_mean) || !is.finite(null_mean)) return(NA_real_)
  if (abs(null_mean) < .Machine$double.eps) return(NA_real_)
  (actual_mean - null_mean) / null_mean
}

classify_interaction <- function(deviation, p_value, model_name, direction = "negative", alpha = 0.05) {
  if (!is.finite(deviation) || !is.finite(p_value) || p_value >= alpha) return(model_name)
  if (direction == "positive") {
    if (deviation > 0) return("Synergistic")
    if (deviation < 0) return("Antagonistic")
  }
  if (direction == "negative") {
    if (deviation > 0) return("Antagonistic")
    if (deviation < 0) return("Synergistic")
  }
  model_name
}

calc_deviation_all_models <- function(response, data, year_label, n_boot = cfg$n_null_boot_deviation) {
  null_models <- c("additive", "multiplicative", "dominative")
  df_H <- data[data$remark %in% h_levels, , drop = FALSE]
  df_CK <- data[data$remark == "CT", , drop = FALSE]
  df_H <- df_H[is.finite(df_H[[response]]), , drop = FALSE]
  df_CK <- df_CK[is.finite(df_CK[[response]]), , drop = FALSE]

  ck_mean <- mean(df_CK[[response]], na.rm = TRUE)

  base <- data.frame(
    year = year_label,
    response = response,
    null_model = NA_character_,
    row_id = df_H$.__rowid,
    treatment = as.character(df_H$remark),
    Lv = df_H$Lv,
    deviation = NA_real_,
    squared_deviation = NA_real_,
    P = NA_real_,
    actual_ES = NA_real_,
    null_ES = NA_real_,
    selected_factors = NA_character_,
    interaction_type = NA_character_,
    stringsAsFactors = FALSE
  )

  out <- setNames(lapply(null_models, function(m) {
    tmp <- base
    tmp$null_model <- m
    tmp
  }), null_models)

  for (i in seq_len(nrow(df_H))) {
    comb_values <- safe_numeric(unlist(df_H[i, stressors, drop = FALSE], use.names = FALSE))
    selected <- stressors[which(comb_values == 1)]
    selected_label <- paste(selected, collapse = "+")

    sim <- simulate_null_effects(
      response = response,
      data = data,
      selected_factors = selected,
      n_boot = n_boot,
      return_response = TRUE
    )

    actual_mean <- df_H[[response]][i]

    for (m in null_models) {
      nd <- sim[[m]]
      null_mean <- mean(nd, na.rm = TRUE)
      dev <- safe_relative_deviation(actual_mean, null_mean)
      p <- safe_tail_p(actual_mean - nd, two_sided = cfg$interaction_p_two_sided)

      out[[m]]$deviation[i] <- dev
      out[[m]]$squared_deviation[i] <- dev^2
      out[[m]]$P[i] <- p
      out[[m]]$actual_ES[i] <- actual_mean - ck_mean
      out[[m]]$null_ES[i] <- null_mean - ck_mean
      out[[m]]$selected_factors[i] <- selected_label
      out[[m]]$interaction_type[i] <- classify_interaction(
        deviation = dev,
        p_value = p,
        model_name = m,
        direction = response_direction[[response]]
      )
    }
  }

  levels_i <- c("Antagonistic", "additive", "multiplicative", "dominative", "Synergistic")
  for (m in null_models) out[[m]]$interaction_type <- factor(out[[m]]$interaction_type, levels = levels_i)
  out
}

bootstrap_deviation_slope <- function(df_dv, n_boot = cfg$n_boot_slope, probs = c(0.025, 0.975)) {
  d <- df_dv %>% dplyr::filter(is.finite(deviation), is.finite(Lv))
  if (nrow(d) < 5 || length(unique(d$Lv)) < 2) {
    return(data.frame(n = nrow(d), slope = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  fit <- lm(deviation ~ Lv, data = d)
  slopes <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(d)), size = nrow(d), replace = TRUE)
    dd <- d[idx, , drop = FALSE]
    if (length(unique(dd$Lv)) < 2) return(NA_real_)
    unname(coef(lm(deviation ~ Lv, data = dd))["Lv"])
  })
  slopes <- slopes[is.finite(slopes)]

  data.frame(
    n = nrow(d),
    slope = unname(coef(fit)["Lv"]),
    ci_low = unname(quantile(slopes, probs[1], na.rm = TRUE)),
    ci_high = unname(quantile(slopes, probs[2], na.rm = TRUE))
  )
}

calc_ssd_summary <- function(deviation_all) {
  deviation_all %>%
    dplyr::filter(is.finite(deviation)) %>%
    dplyr::group_by(year, response, null_model) %>%
    dplyr::summarise(
      n_deviation = dplyr::n(),
      SSD = sum(deviation^2, na.rm = TRUE),
      mean_squared_deviation = mean(deviation^2, na.rm = TRUE),
      root_mean_squared_deviation = sqrt(mean(deviation^2, na.rm = TRUE)),
      mean_abs_deviation = mean(abs(deviation), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(year, response) %>%
    dplyr::mutate(
      SSD_rank = rank(SSD, ties.method = "min", na.last = "keep"),
      best_model_by_SSD = SSD_rank == 1
    ) %>%
    dplyr::ungroup()
}

calc_deviation_lm_summary <- function(deviation_all) {
  split_list <- split(deviation_all, list(deviation_all$year, deviation_all$response, deviation_all$null_model), drop = TRUE)
  out <- lapply(split_list, function(d) {
    yr <- unique(d$year)[1]
    resp <- unique(d$response)[1]
    mdl <- unique(d$null_model)[1]
    d <- d %>% dplyr::filter(is.finite(deviation), is.finite(Lv))
    if (nrow(d) < 3 || length(unique(d$Lv)) < 2) {
      return(data.frame(year = yr, response = resp, null_model = mdl, n = nrow(d),
                        intercept = NA_real_, slope = NA_real_, slope_se = NA_real_,
                        slope_t = NA_real_, slope_p = NA_real_, slope_ci_low = NA_real_,
                        slope_ci_high = NA_real_, r_squared = NA_real_, adj_r_squared = NA_real_,
                        model_p = NA_real_))
    }
    fit <- lm(deviation ~ Lv, data = d)
    sm <- summary(fit)
    cf <- sm$coefficients
    ci <- confint(fit, level = 0.95)
    fstat <- sm$fstatistic
    model_p <- if (!is.null(fstat)) pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) else NA_real_
    data.frame(
      year = yr,
      response = resp,
      null_model = mdl,
      n = nrow(d),
      intercept = unname(coef(fit)["(Intercept)"]),
      slope = unname(coef(fit)["Lv"]),
      slope_se = cf["Lv", "Std. Error"],
      slope_t = cf["Lv", "t value"],
      slope_p = cf["Lv", "Pr(>|t|)"],
      slope_ci_low = ci["Lv", 1],
      slope_ci_high = ci["Lv", 2],
      r_squared = sm$r.squared,
      adj_r_squared = sm$adj.r.squared,
      model_p = model_p
    )
  })
  dplyr::bind_rows(out)
}

##############################
# 5. Plotting functions      #
##############################

plot_es_heatmap <- function(df_hm, response_name) {
  df_hm$treatment <- factor(df_hm$treatment, levels = treat_order_16)
  df_hm$year <- factor(df_hm$year, levels = sort(unique(df_hm$year)))

  ggplot(df_hm, aes(x = treatment, y = year, fill = ES_mean)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = stars), size = 3) +
    scale_y_discrete(limits = rev(levels(df_hm$year))) +
    scale_x_discrete(labels = c(stressor_labels, "1" = "1", "2" = "2", "4" = "4", "6" = "6", "8" = "8", "10" = "10")) +
    scale_fill_gradient2(
      low = "#8C510A",
      mid = "#F2EFE2",
      high = "#01665E",
      midpoint = 0,
      na.value = "grey90",
      name = "Mean ES\n(TR - CT)"
    ) +
    labs(x = NULL, y = NULL, title = response_label(response_name)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 10)
    )
}

plot_factor_number_effects <- function(tbl, response_name) {
  model_i <- unique(tbl$interaction_model)
  ggplot(tbl, aes(x = year, y = slope, group = 1)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.6) +
    geom_point(size = 2.4) +
    scale_x_continuous(breaks = sort(unique(tbl$year))) +
    labs(
      x = NULL,
      y = paste0("Standardized effect of\nfactor number on\n", tolower(response_label(response_name))),
      title = paste0("Calculated by deviations from ", model_i, " model")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 9),
      axis.text.x = element_text(angle = 0)
    )
}

make_raw_plot <- function(df, response_mean, response_name) {
  fill_cols <- remark_cols
  names(fill_cols) <- names(remark_cols)
  fill_cols <- paste0(fill_cols, "80")

  ct_mean <- response_mean$mean[response_mean$target == "CT"][1]

  ggplot() +
    theme_bw(base_size = 9) +
    coord_flip() +
    stat_density_ridges(
      data = df,
      aes_string(x = response_name, y = "remark", fill = "remark"),
      geom = "density_ridges_gradient",
      rel_min_height = 0.01,
      jittered_points = TRUE,
      color = "#00000000",
      position = position_points_jitter(height = 0.2, yoffset = 0.15),
      point_size = 0.8,
      point_alpha = 0.25,
      scale = 0.5
    ) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    geom_estci(
      data = response_mean,
      aes(x = mean, y = target, xmin = X2.5., xmax = X97.5., xintercept = ct_mean, color = target),
      center.linecolour = "black",
      size = 0.5,
      ci.linesize = 0.4,
      position = position_nudge(y = -0.15)
    ) +
    scale_color_manual(values = remark_cols, drop = FALSE) +
    labs(x = response_label(response_name), y = NULL) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "#4D728510", color = "white"),
      panel.grid.major = element_line(color = "white", linewidth = 0.1),
      panel.grid.minor = element_line(color = "white", linewidth = 0.2)
    )
}

make_deviation_plot <- function(deviation_df, ssd_row = NULL, lm_row = NULL) {
  best_model <- FALSE
  ssd_value <- NA_real_
  if (!is.null(ssd_row) && nrow(ssd_row) > 0) {
    best_model <- isTRUE(ssd_row$best_model_by_SSD[1])
    ssd_value <- ssd_row$SSD[1]
  }

  slope_txt <- "NA"
  p_txt <- "NA"
  r2_txt <- "NA"
  if (!is.null(lm_row) && nrow(lm_row) > 0) {
    slope_txt <- ifelse(is.finite(lm_row$slope[1]), format(round(lm_row$slope[1], 4), nsmall = 4), "NA")
    p_txt <- ifelse(is.finite(lm_row$slope_p[1]), ifelse(lm_row$slope_p[1] < 0.001, "<0.001", format(round(lm_row$slope_p[1], 4), nsmall = 4)), "NA")
    r2_txt <- ifelse(is.finite(lm_row$r_squared[1]), format(round(lm_row$r_squared[1], 3), nsmall = 3), "NA")
  }

  ggplot(deviation_df, aes(x = Lv, y = deviation)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = ifelse(best_model, "#FFF3B030", "#FFFFFF00")) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey70", alpha = 0.35, linewidth = 0.7) +
    geom_point(aes(color = interaction_type), size = 1.5, alpha = 0.85, position = position_jitter(width = 0.12, height = 0)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
    scale_color_manual(
      values = c(
        "Antagonistic" = "#5499C7",
        "additive" = "#909497",
        "multiplicative" = "#909497",
        "dominative" = "#909497",
        "Synergistic" = "#F16667"
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0(unique(deviation_df$year), " | ", unique(deviation_df$null_model)),
      subtitle = paste0("SSD=", round(ssd_value, 3), "; slope=", slope_txt, "; P=", p_txt, "; R2=", r2_txt),
      x = "Number of factors",
      y = "Deviation from null model"
    ) +
    theme_bw(base_size = 8) +
    theme(legend.position = "none", plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6))
}

##############################
# 6. Main Figure 2 analysis  #
##############################

run_main_figure2 <- function(dat) {
  year_list <- make_year_list(dat)

  yearly_es_heatmap <- data.frame()
  yearly_responses <- data.frame()
  yearly_factor_effect <- data.frame()
  raw_plots <- list()
  expected_effect_tables <- list()

  for (yr in names(year_list)) {
    message("Running main field analysis for ", yr)
    df <- prep_year_df(year_list[[yr]])
    raw_plots[[yr]] <- list()
    expected_effect_tables[[yr]] <- list()

    for (resp in responses) {
      response_mean <- bootstrap_group_means(resp, df, remark_order, n_boot = cfg$n_boot_mean)
      raw_plots[[yr]][[resp]] <- make_raw_plot(df, response_mean, resp)

      yearly_responses <- rbind(
        yearly_responses,
        data.frame(
          year = as.numeric(yr),
          response = resp,
          response_label = response_label(resp),
          treatment = response_mean$target,
          Respon_2.5 = response_mean$X2.5.,
          Respon_mean = response_mean$mean,
          Respon_97.5 = response_mean$X97.5.,
          stringsAsFactors = FALSE
        )
      )

      # In the original script, expected null effects were calculated before
      # treatment-vs-control bootstrap tests. Because both steps are stochastic,
      # changing this order can slightly alter BH-adjusted heatmap stars when
      # p-values are close to thresholds. The legacy branch restores that order.
      if (isTRUE(cfg$heatmap_legacy_rng_order)) {
        null_dist_main <- rf_null_distribution_rep(
          response = resp,
          data = df,
          sub_n_perm = cfg$n_expected_boot_per_combo_main
        )
        expected_tbl <- rf_expected_effects_for_rows(null_dist_main, df)
      } else {
        expected_tbl <- calc_expected_effects_for_rows(resp, df, n_boot = cfg$n_null_boot_expected)
      }
      expected_effect_tables[[yr]][[resp]] <- expected_tbl

      es_tbl <- bootstrap_effect_vs_control(resp, df, treat_order_16, n_boot = cfg$n_boot_effect)
      es_tbl$p_adj <- p.adjust(es_tbl$p_value, method = cfg$heatmap_p_adjust_method)
      es_tbl$stars <- sig_stars(es_tbl$p_adj)

      yearly_es_heatmap <- rbind(
        yearly_es_heatmap,
        data.frame(
          year = as.numeric(yr),
          response = resp,
          response_label = response_label(resp),
          es_tbl,
          p_adjust_method = cfg$heatmap_p_adjust_method,
          stringsAsFactors = FALSE
        )
      )

      dev_models <- calc_deviation_all_models(resp, df, year_label = yr, n_boot = cfg$n_null_boot_deviation)
      model_i <- final_null_model_by_response[[resp]]
      slope_i <- bootstrap_deviation_slope(dev_models[[model_i]], n_boot = cfg$n_boot_slope)

      yearly_factor_effect <- rbind(
        yearly_factor_effect,
        data.frame(
          year = as.numeric(yr),
          response = resp,
          response_label = response_label(resp),
          interaction_model = model_i,
          slope = slope_i$slope,
          ci_low = slope_i$ci_low,
          ci_high = slope_i$ci_high,
          n = slope_i$n,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  heatmap_plots <- lapply(responses, function(resp) {
    plot_es_heatmap(yearly_es_heatmap %>% dplyr::filter(response == resp), resp)
  })
  names(heatmap_plots) <- responses

  factor_plots <- lapply(responses, function(resp) {
    plot_factor_number_effects(yearly_factor_effect %>% dplyr::filter(response == resp), resp)
  })
  names(factor_plots) <- responses

  fig2_rows <- lapply(responses, function(resp) {
    heatmap_plots[[resp]] + factor_plots[[resp]] + plot_layout(widths = c(5, 4))
  })

  fig2 <- wrap_plots(fig2_rows, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a")

  if (cfg$save_outputs) {
    fig_height <- max(6, 3.2 * length(responses))
    ggsave(file.path(cfg$output_dir, "figures", "Figure2_field_experiment.pdf"), fig2, width = 12, height = fig_height)
    ggsave(file.path(cfg$output_dir, "figures", "Figure2_field_experiment.tiff"), fig2, width = 12, height = fig_height, dpi = 350)
    write.csv(yearly_es_heatmap, file.path(cfg$output_dir, "tables", "Figure2_heatmap_effect_sizes.csv"), row.names = FALSE)
    write.csv(yearly_responses, file.path(cfg$output_dir, "tables", "field_bootstrapped_group_means.csv"), row.names = FALSE)
    write.csv(yearly_factor_effect, file.path(cfg$output_dir, "tables", "Figure2_factor_number_effects.csv"), row.names = FALSE)
    save(raw_plots, heatmap_plots, factor_plots, expected_effect_tables,
         file = file.path(cfg$output_dir, "objects", "field_main_figure_objects.RData"))
  }

  list(
    figure2 = fig2,
    raw_plots = raw_plots,
    heatmap_plots = heatmap_plots,
    factor_plots = factor_plots,
    yearly_es_heatmap = yearly_es_heatmap,
    yearly_responses = yearly_responses,
    yearly_factor_effect = yearly_factor_effect,
    expected_effect_tables = expected_effect_tables
  )
}

#######################################
# 7. Supplementary raw ridge figures  #
#######################################

save_raw_ridge_supplements <- function(raw_plots) {
  if (is.null(raw_plots) || length(raw_plots) == 0) return(invisible(NULL))

  years <- sort(names(raw_plots))

  all_plot_list <- unlist(lapply(years, function(yr) {
    lapply(responses, function(resp) {
      raw_plots[[yr]][[resp]] + ggtitle(paste(yr, response_label(resp)))
    })
  }), recursive = FALSE)

  fig_all <- wrap_plots(all_plot_list, ncol = length(responses)) +
    plot_annotation(tag_levels = "a")

  out <- list(FigureS_raw_all_responses = fig_all)

  # Preserve the original Figure S2/S3 split when those four responses are present.
  if (all(c("richness_all", "biomass") %in% responses)) {
    fig_s2_list <- unlist(lapply(years, function(yr) {
      list(raw_plots[[yr]][["richness_all"]], raw_plots[[yr]][["biomass"]])
    }), recursive = FALSE)
    out$FigureS2 <- wrap_plots(fig_s2_list, ncol = 2) + plot_annotation(tag_levels = "a")
  }

  if (all(c("shannon_all", "evenness_all") %in% responses)) {
    fig_s3_list <- unlist(lapply(years, function(yr) {
      list(raw_plots[[yr]][["shannon_all"]], raw_plots[[yr]][["evenness_all"]])
    }), recursive = FALSE)
    out$FigureS3 <- wrap_plots(fig_s3_list, ncol = 2) + plot_annotation(tag_levels = "a")
  }

  if (cfg$save_outputs) {
    ggsave(
      file.path(cfg$output_dir, "figures", "FigureS_raw_ridge_all_selected_responses.pdf"),
      fig_all,
      width = max(8, 4.5 * length(responses)),
      height = max(6, 3 * length(years))
    )
    if (!is.null(out$FigureS2)) {
      ggsave(file.path(cfg$output_dir, "figures", "FigureS2_raw_ridge_richness_biomass.pdf"), out$FigureS2, width = 10, height = 15)
    }
    if (!is.null(out$FigureS3)) {
      ggsave(file.path(cfg$output_dir, "figures", "FigureS3_raw_ridge_shannon_evenness.pdf"), out$FigureS3, width = 10, height = 15)
    }
  }

  invisible(out)
}

#######################################
# 8. Supplementary NMDS/PERMANOVA     #
#######################################

non_species_cols <- c(
  "SampleNo", "time", stressors, "Lv", "remark", "year", ".__rowid",
  "grass_sum", "herb_sum", "legume_sum", "all_abundance",
  "grass_prop", "herb_prop", "legume_prop",
  "richness_all", "diversity_all", "diversity_grass", "diversity_herb", "diversity_legume",
  "shannon_all", "evenness_all", "simpson_all",
  "shannon_grass", "evenness_grass", "simpson_grass",
  "shannon_herb", "evenness_herb", "simpson_herb",
  "shannon_legume", "evenness_legume", "simpson_legume", "biomass"
)

identify_species_cols <- function(dat) {
  candidates <- setdiff(colnames(dat), unique(c(non_species_cols, responses)))
  candidates <- candidates[vapply(dat[, candidates, drop = FALSE], function(x) {
    xx <- safe_numeric(x)
    all(is.na(xx) | is.finite(xx))
  }, logical(1))]
  candidates
}

prepare_nmds_data <- function(data, selected_year, species_cols, transform_to_relative = TRUE) {
  df <- data[data$year == selected_year, , drop = FALSE]
  if (nrow(df) == 0) stop("No rows found for year ", selected_year)

  community <- df[, species_cols, drop = FALSE]
  community[] <- lapply(community, safe_numeric)
  community[is.na(community)] <- 0

  keep_samples <- rowSums(community, na.rm = TRUE) > 0
  community <- community[keep_samples, , drop = FALSE]
  groups <- df[keep_samples, , drop = FALSE]

  keep_species <- colSums(community, na.rm = TRUE) > 0
  community <- community[, keep_species, drop = FALSE]

  if (nrow(community) < 3) stop("Fewer than 3 non-empty samples for year ", selected_year)
  if (ncol(community) < 2) stop("Fewer than 2 non-zero species columns for year ", selected_year)

  if (transform_to_relative) community <- vegan::decostand(community, method = "total")

  list(community = community, groups = groups, species_used = colnames(community))
}

run_nmds <- function(community_matrix) {
  set.seed(123)
  vegan::metaMDS(
    community_matrix,
    distance = "bray",
    k = 2,
    maxit = cfg$nmds_maxit,
    trymax = cfg$nmds_trymax,
    wascores = TRUE,
    autotransform = FALSE,
    trace = FALSE
  )
}

extract_nmds_scores <- function(nmds_fit, groups) {
  site_scores <- as.data.frame(vegan::scores(nmds_fit, display = "sites"))
  site_scores$SampleNo <- groups$SampleNo
  site_scores$time <- groups$time
  site_scores$year <- groups$year
  site_scores$Lv <- as.character(groups$Lv)
  site_scores$remark <- as.character(groups$remark)
  site_scores
}

plot_nmds_by_lv <- function(nmds_fit, groups, year_label) {
  site_scores <- extract_nmds_scores(nmds_fit, groups)
  site_scores$Lv <- factor(site_scores$Lv, levels = names(lv_cols))

  centroids <- site_scores %>%
    dplyr::group_by(Lv) %>%
    dplyr::summarise(cNMDS1 = mean(NMDS1, na.rm = TRUE), cNMDS2 = mean(NMDS2, na.rm = TRUE), .groups = "drop")

  site_scores <- site_scores %>% dplyr::left_join(centroids, by = "Lv")

  ggplot(site_scores, aes(x = NMDS1, y = NMDS2, color = Lv)) +
    geom_segment(aes(xend = cNMDS1, yend = cNMDS2), alpha = 0.35, linewidth = 0.25) +
    geom_point(size = 1.7, alpha = 0.85) +
    geom_point(data = centroids, aes(x = cNMDS1, y = cNMDS2, color = Lv), shape = 15, size = 3, inherit.aes = FALSE) +
    scale_color_manual(values = lv_cols, drop = FALSE) +
    labs(
      title = year_label,
      subtitle = paste0("R2 and P from PERMANOVA; NMDS stress = ", round(nmds_fit$stress, 3)),
      x = "NMDS1",
      y = "NMDS2",
      color = "Number of factors"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}

run_permanova_by_lv <- function(community_matrix, groups) {
  groups$Lv <- factor(groups$Lv)
  vegan::adonis2(community_matrix ~ Lv, data = groups, method = "bray", permutations = cfg$permanova_permutations)
}

run_permanova_full <- function(community_matrix, groups) {
  groups$Lv <- factor(groups$Lv)
  vegan::adonis2(
    community_matrix ~ Lv + W + N + D + HM + MP + S + A + I + AF + SF,
    data = groups,
    method = "bray",
    permutations = cfg$permanova_permutations
  )
}

permanova_to_df <- function(x, year_label, test_name) {
  out <- as.data.frame(x)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out$year <- year_label
  out$test <- test_name
  out[, c("year", "test", "term", setdiff(colnames(out), c("year", "test", "term")))]
}

run_nmds_supplement <- function(dat) {
  species_cols <- identify_species_cols(dat)
  message("Detected ", length(species_cols), " species columns for NMDS.")

  years <- sort(unique(na.omit(dat$year)))
  nmds_scores <- list()
  nmds_plots <- list()
  permanova_Lv <- list()
  permanova_full <- list()

  for (yr in years) {
    message("Running NMDS/PERMANOVA for ", yr)
    inp <- prepare_nmds_data(dat, selected_year = yr, species_cols = species_cols)
    fit <- run_nmds(inp$community)
    nmds_scores[[yr]] <- extract_nmds_scores(fit, inp$groups)
    permanova_Lv[[yr]] <- run_permanova_by_lv(inp$community, inp$groups)
    permanova_full[[yr]] <- run_permanova_full(inp$community, inp$groups)
    p <- plot_nmds_by_lv(fit, inp$groups, yr)

    pval <- as.data.frame(permanova_Lv[[yr]])
    r2 <- pval$R2[1]
    pv <- pval$`Pr(>F)`[1]
    p <- p + labs(subtitle = paste0("R2 = ", round(r2, 3), ", P = ", signif(pv, 3), "; stress = ", round(fit$stress, 3)))
    nmds_plots[[yr]] <- p
  }

  fig_s4 <- wrap_plots(nmds_plots, ncol = 3) + plot_annotation(tag_levels = "a")

  permanova_Lv_table <- dplyr::bind_rows(lapply(names(permanova_Lv), function(yr) permanova_to_df(permanova_Lv[[yr]], yr, "Lv")))
  permanova_full_table <- dplyr::bind_rows(lapply(names(permanova_full), function(yr) permanova_to_df(permanova_full[[yr]], yr, "Lv_plus_stressors")))
  nmds_scores_table <- dplyr::bind_rows(nmds_scores)

  if (cfg$save_outputs) {
    ggsave(file.path(cfg$output_dir, "figures", "FigureS4_NMDS_PERMANOVA.pdf"), fig_s4, width = 12, height = 8)
    write.csv(nmds_scores_table, file.path(cfg$output_dir, "tables", "NMDS_site_scores_all_years.csv"), row.names = FALSE)
    write.csv(permanova_Lv_table, file.path(cfg$output_dir, "tables", "PERMANOVA_by_Lv_all_years.csv"), row.names = FALSE)
    write.csv(permanova_full_table, file.path(cfg$output_dir, "tables", "PERMANOVA_Lv_plus_stressors_all_years.csv"), row.names = FALSE)
    save(nmds_plots, nmds_scores, permanova_Lv, permanova_full,
         file = file.path(cfg$output_dir, "objects", "NMDS_PERMANOVA_results.RData"))
  }

  list(figureS4 = fig_s4, nmds_scores = nmds_scores_table, permanova_Lv = permanova_Lv_table, permanova_full = permanova_full_table)
}

##########################################
# 9. Supplementary deviation/SSD plots   #
##########################################

run_deviation_supplement <- function(dat) {
  years <- sort(unique(na.omit(dat$year)))
  deviation_all <- data.frame()

  for (yr in years) {
    message("Running supplementary deviation analysis for ", yr)
    df <- prep_year_df(dat[dat$year == yr, , drop = FALSE])
    for (resp in responses) {
      dev_models <- calc_deviation_all_models(resp, df, year_label = yr, n_boot = cfg$n_null_boot_deviation)
      deviation_all <- rbind(deviation_all, dplyr::bind_rows(dev_models))
    }
  }

  ssd_summary <- calc_ssd_summary(deviation_all)
  lm_summary <- calc_deviation_lm_summary(deviation_all)

  best_model <- ssd_summary %>%
    dplyr::filter(best_model_by_SSD) %>%
    dplyr::select(year, response, null_model, SSD, mean_squared_deviation, root_mean_squared_deviation, mean_abs_deviation)

  null_models <- c("additive", "multiplicative", "dominative")
  deviation_plots <- list()

  for (resp in responses) {
    plot_list <- list()
    for (yr in years) {
      for (mdl in null_models) {
        d <- deviation_all %>% dplyr::filter(response == resp, year == yr, null_model == mdl)
        ssd_row <- ssd_summary %>% dplyr::filter(response == resp, year == yr, null_model == mdl)
        lm_row <- lm_summary %>% dplyr::filter(response == resp, year == yr, null_model == mdl)
        plot_list[[paste(yr, mdl, sep = "_")]] <- make_deviation_plot(d, ssd_row, lm_row)
      }
    }
    deviation_plots[[resp]] <- wrap_plots(plot_list, ncol = 3) + plot_annotation(title = response_label(resp))
  }

  if (cfg$save_outputs) {
    for (resp in responses) {
      ggsave(
        file.path(cfg$output_dir, "figures", paste0("FigureS_deviation_SSD_", resp, ".pdf")),
        deviation_plots[[resp]], width = 12, height = 15
      )
    }
    write.csv(deviation_all, file.path(cfg$output_dir, "tables", "deviation_all_years.csv"), row.names = FALSE)
    write.csv(ssd_summary, file.path(cfg$output_dir, "tables", "ssd_summary_all_years.csv"), row.names = FALSE)
    write.csv(lm_summary, file.path(cfg$output_dir, "tables", "deviation_lm_summary_all_years.csv"), row.names = FALSE)
    write.csv(best_model, file.path(cfg$output_dir, "tables", "best_null_model_by_year_response.csv"), row.names = FALSE)
    save(deviation_plots, deviation_all, ssd_summary, lm_summary, best_model,
         file = file.path(cfg$output_dir, "objects", "deviation_SSD_results.RData"))
  }

  list(deviation_all = deviation_all, ssd_summary = ssd_summary, lm_summary = lm_summary, best_model = best_model, plots = deviation_plots)
}

##########################################
# 10. Supplementary random forest R2     #
##########################################

#   1. Expected factor-identity effects are calculated per treatment combination.
#   2. Each RF iteration first builds a balanced bootstrap dataset:
#        - n_eachlv one-factor rows sampled with replacement
#        - all multi-factor rows sampled with replacement
#   3. R2 is calculated from party::cforest out-of-bag predictions on that
#      bootstrap dataset: predict(fit, OOB = TRUE).
#   4. Lv, ES, LvES, and All models are evaluated on the same bootstrap rows
#      within each iteration.

bt <- function(x) paste0("`", x, "`")

make_formula <- function(response, predictors) {
  as.formula(paste(bt(response), "~", paste(bt(predictors), collapse = " + ")))
}

make_rf_controls <- function(n_tree) {
  ctrl <- tryCatch(
    party::cforest_control(ntree = n_tree, minsplit = 5, minbucket = 2),
    error = function(e) NULL
  )
  if (!is.null(ctrl)) return(ctrl)

  ctrl <- tryCatch(
    party::cforest_unbiased(ntree = n_tree, minsplit = 5, minbucket = 2),
    error = function(e) NULL
  )
  if (!is.null(ctrl)) return(ctrl)

  party::cforest_unbiased(ntree = n_tree)
}

safe_r2 <- function(pred, obs) {
  pred <- safe_numeric(pred)
  obs <- safe_numeric(obs)

  keep <- is.finite(pred) & is.finite(obs)
  if (sum(keep) < 3) return(NA_real_)
  if (length(unique(obs[keep])) < 2) return(NA_real_)
  if (length(unique(pred[keep])) < 2) return(NA_real_)

  as.numeric(
    tryCatch(
      caret::postResample(pred = pred[keep], obs = obs[keep])[["Rsquared"]],
      error = function(e) NA_real_
    )
  )
}

fit_rf_oob_r2 <- function(fml, rdf, response_name, n_tree) {
  needed <- all.vars(fml)
  rdf <- rdf[stats::complete.cases(rdf[, needed, drop = FALSE]), , drop = FALSE]

  if (nrow(rdf) < 5) return(NA_real_)
  if (length(unique(rdf[[response_name]])) < 2) return(NA_real_)

  fit <- tryCatch(
    party::cforest(
      formula = fml,
      data = rdf,
      controls = make_rf_controls(n_tree)
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)

  pred <- tryCatch(
    predict(fit, OOB = TRUE),
    error = function(e) rep(NA_real_, nrow(rdf))
  )

  safe_r2(pred, rdf[[response_name]])
}

summary_ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(CI.low = NA_real_, Mean = NA_real_, CI.high = NA_real_))
  }

  c(
    CI.low = unname(stats::quantile(x, 0.025, na.rm = TRUE)),
    Mean = unname(stats::quantile(x, 0.500, na.rm = TRUE)),
    CI.high = unname(stats::quantile(x, 0.975, na.rm = TRUE))
  )
}

rf_null_distribution_rep <- function(response, data, sub_n_perm = cfg$n_rf_expected_boot_per_combo) {
  output <- list()

  population_CT <- data[data$remark == "CT", response]
  population_CT <- population_CT[is.finite(population_CT)]
  size_CT <- length(population_CT)

  for (Lv in lv_levels) {
    combination <- data[data$remark == Lv, stressors, drop = FALSE]

    resampled <- list(
      Additive = vector("list", nrow(combination)),
      Multiplicative = vector("list", nrow(combination)),
      Dominative = vector("list", nrow(combination))
    )

    if (nrow(combination) == 0 || size_CT < 1) {
      output[[Lv]] <- resampled
      next
    }

    for (type in c("Additive", "Multiplicative", "Dominative")) {
      for (j in seq_len(nrow(combination))) {
        bs <- numeric(0)
        comb_row <- safe_numeric(unlist(combination[j, ], use.names = FALSE))
        selected_stressors <- stressors[which(comb_row == 1)]

        if (length(selected_stressors) == 0) {
          resampled[[type]][[j]] <- bs
          next
        }

        for (id in seq_len(sub_n_perm)) {
          k_CT <- mean(sample(population_CT, size_CT, replace = TRUE), na.rm = TRUE)
          if (!is.finite(k_CT)) next

          each_effect <- numeric(0)

          for (trt in selected_stressors) {
            population_TR <- data[data$remark == trt, response]
            population_TR <- population_TR[is.finite(population_TR)]
            size_TR <- length(population_TR)
            if (size_TR < 1) next

            k_TR <- mean(sample(population_TR, size_TR, replace = TRUE), na.rm = TRUE)
            if (!is.finite(k_TR)) next

            if (type == "Additive") each_effect <- c(each_effect, k_TR - k_CT)
            if (type == "Multiplicative") {
              if (abs(k_CT) < .Machine$double.eps) next
              each_effect <- c(each_effect, (k_TR - k_CT) / k_CT)
            }
            if (type == "Dominative") each_effect <- c(each_effect, k_TR - k_CT)
          }

          if (length(each_effect) == 0) next

          if (type == "Additive") joint_effect <- sum(each_effect)
          if (type == "Multiplicative") joint_effect <- (prod(1 + each_effect) - 1) * k_CT
          if (type == "Dominative") joint_effect <- each_effect[which.max(abs(each_effect))]

          if (is.finite(joint_effect)) bs <- c(bs, joint_effect)
        }

        resampled[[type]][[j]] <- bs
      }
    }

    output[[Lv]] <- resampled
  }

  output
}

rf_expected_effects_for_rows <- function(null_dist, data) {
  df_lv <- data[data$remark %in% lv_levels, , drop = FALSE]

  if (nrow(df_lv) == 0) {
    return(data.frame(.__rowid = integer(0), E1 = numeric(0), E2 = numeric(0), E3 = numeric(0)))
  }

  key_map <- list()

  for (Lv in lv_levels) {
    if (Lv == "1") {
      comb <- data[data$remark %in% stressors, stressors, drop = FALSE]
    } else {
      comb <- data[data$remark == Lv, stressors, drop = FALSE]
    }

    if (nrow(comb) == 0) next

    comb_key <- make_stressor_key(comb)
    key_map[[Lv]] <- stats::setNames(seq_along(comb_key), comb_key)
  }

  row_key <- make_stressor_key(df_lv[, stressors, drop = FALSE])
  row_lv <- as.character(df_lv$Lv)
  row_lv[!(row_lv %in% lv_levels)] <- as.character(df_lv$remark[!(row_lv %in% lv_levels)])

  E1 <- E2 <- E3 <- rep(NA_real_, nrow(df_lv))

  for (i in seq_len(nrow(df_lv))) {
    Lv <- row_lv[i]
    if (!Lv %in% names(key_map)) next
    if (!Lv %in% names(null_dist)) next

    j <- unname(key_map[[Lv]][row_key[i]])
    if (is.na(j)) next

    add <- null_dist[[Lv]][["Additive"]][[j]]
    mul <- null_dist[[Lv]][["Multiplicative"]][[j]]
    dom <- null_dist[[Lv]][["Dominative"]][[j]]

    E1[i] <- if (length(add)) mean(add, na.rm = TRUE) else NA_real_
    E2[i] <- if (length(mul)) mean(mul, na.rm = TRUE) else NA_real_
    E3[i] <- if (length(dom)) mean(dom, na.rm = TRUE) else NA_real_
  }

  data.frame(.__rowid = df_lv$.__rowid, E1 = E1, E2 = E2, E3 = E3)
}

make_rf_r2_plot <- function(rf_raw, rf_summary, response_name, year_label = NULL) {
  level_order <- c("Lv", "ES", "LvES", "All")

  rf_raw$Model <- factor(rf_raw$Model, levels = level_order)
  rf_summary$Model <- factor(rf_summary$Model, levels = level_order)

  cols <- c(
    Lv = "#E69F37",
    ES = "#E69F37",
    LvES = "#8D2F1E",
    All = "#9C9E83"
  )

  ggplot(rf_raw, aes(x = Model, y = R2, fill = Model)) +
    geom_violin(color = "#00000000", alpha = 0.5, trim = FALSE) +
    geom_pointrange(
      data = rf_summary,
      aes(y = Mean, ymin = CI.low, ymax = CI.high, color = Model),
      fatten = 3,
      linewidth = 0.3
    ) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_color_manual(values = cols, drop = FALSE) +
    scale_x_discrete(
      labels = c(
        Lv = "Number of factors",
        ES = "Factor identity effect",
        LvES = "Number of factors\n+ factor identity effect",
        All = "Factor composition"
      ),
      drop = FALSE
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = response_label(response_name),
      subtitle = paste0("ES predictor: ", final_null_model_by_response[[response_name]]),
      x = NULL,
      y = "Out-of-bag R2"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "#6D6E7130", color = "white"),
      panel.grid.major = element_line(color = "white", linewidth = 0.2),
      panel.grid.minor = element_line(color = "white", linewidth = 0.2),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

run_rf_year <- function(df, expected_effect_tables, year_label) {
  model_levels <- c("Lv", "ES", "LvES", "All")

  out_raw <- data.frame()
  out_summary <- data.frame()
  plots <- list()

  for (resp in responses) {
    df_rf <- df[df$remark %in% lv_levels, , drop = FALSE] %>%
      dplyr::left_join(expected_effect_tables[[resp]], by = ".__rowid")

    ES_col <- expected_effect_column_by_response[[resp]]
    if (is.null(ES_col) || is.na(ES_col)) stop("Missing RF ES column for response: ", resp)
    if (!ES_col %in% names(df_rf)) stop("Expected ES column not found: ", ES_col)

    df_rf$E_RF <- safe_numeric(df_rf[[ES_col]])
    df_rf$Lv <- safe_numeric(df_rf$Lv)

    for (s in stressors) df_rf[[s]] <- safe_numeric(df_rf[[s]])
    df_rf[[resp]] <- safe_numeric(df_rf[[resp]])

    needed_all <- unique(c(resp, "Lv", "E_RF", stressors))
    df_rf <- df_rf[stats::complete.cases(df_rf[, needed_all, drop = FALSE]), , drop = FALSE]
    df_rf <- df_rf[is.finite(df_rf[[resp]]) & is.finite(df_rf$Lv) & is.finite(df_rf$E_RF), , drop = FALSE]

    id_lv1 <- which(df_rf$Lv == 1)
    id_lvh <- which(df_rf$Lv > 1)

    if (length(id_lv1) < 1 || length(id_lvh) < 1 || nrow(df_rf) < 5 || length(unique(df_rf[[resp]])) < 2) {
      warning("Skipping RF for ", year_label, ", ", resp, ": insufficient usable rows.")
      rf_raw_resp <- data.frame(
        year = year_label,
        response = resp,
        Model = factor(rep(model_levels, each = cfg$n_rf_iter), levels = model_levels),
        Iter = rep(seq_len(cfg$n_rf_iter), times = length(model_levels)),
        R2 = NA_real_,
        ES_source = ES_col,
        ES_model = final_null_model_by_response[[resp]],
        stringsAsFactors = FALSE
      )

      sm <- rf_raw_resp %>%
        dplyr::group_by(response, Model, ES_source, ES_model) %>%
        dplyr::summarise(
          CI.low = summary_ci(R2)["CI.low"],
          Mean = summary_ci(R2)["Mean"],
          CI.high = summary_ci(R2)["CI.high"],
          .groups = "drop"
        )

      out_raw <- rbind(out_raw, rf_raw_resp)
      out_summary <- rbind(out_summary, sm)
      plots[[resp]] <- make_rf_r2_plot(rf_raw_resp, sm, resp, year_label)
      next
    }

    fml <- list(
      Lv = make_formula(resp, "Lv"),
      ES = make_formula(resp, "E_RF"),
      LvES = make_formula(resp, c("Lv", "E_RF")),
      All = make_formula(resp, c("Lv", "E_RF", stressors))
    )

    rf_raw_resp <- data.frame(
      year = year_label,
      response = resp,
      Iter = rep(seq_len(cfg$n_rf_iter), each = length(model_levels)),
      Model = factor(rep(model_levels, times = cfg$n_rf_iter), levels = model_levels),
      R2 = NA_real_,
      ES_source = ES_col,
      ES_model = final_null_model_by_response[[resp]],
      stringsAsFactors = FALSE
    )

    row_id <- 0

    for (iter_i in seq_len(cfg$n_rf_iter)) {
      set.seed(
        cfg$rf_seed +
          as.integer(year_label) * 10000 +
          match(resp, responses) * 1000 +
          iter_i
      )

      rid <- c(
        sample(id_lv1, cfg$n_eachlv_rf, replace = TRUE),
        sample(id_lvh, length(id_lvh), replace = TRUE)
      )

      rdf <- df_rf[rid, , drop = FALSE]

      for (model_name in model_levels) {
        row_id <- row_id + 1
        rf_raw_resp$R2[row_id] <- fit_rf_oob_r2(
          fml = fml[[model_name]],
          rdf = rdf,
          response_name = resp,
          n_tree = cfg$n_rf_trees
        )
      }

      if (iter_i %% 100 == 0 || iter_i == cfg$n_rf_iter) {
        message("  ", year_label, " | ", resp, " | finished RF iteration ", iter_i, " / ", cfg$n_rf_iter)
      }
    }

    sm <- do.call(
      rbind,
      lapply(
        model_levels,
        function(mm) {
          ss <- summary_ci(rf_raw_resp$R2[rf_raw_resp$Model == mm])
          data.frame(
            year = year_label,
            response = resp,
            Model = factor(mm, levels = model_levels),
            CI.low = ss["CI.low"],
            Mean = ss["Mean"],
            CI.high = ss["CI.high"],
            ES_source = ES_col,
            ES_model = final_null_model_by_response[[resp]],
            stringsAsFactors = FALSE
          )
        }
      )
    )

    out_raw <- rbind(out_raw, rf_raw_resp)
    out_summary <- rbind(out_summary, sm)
    plots[[resp]] <- make_rf_r2_plot(rf_raw_resp, sm, resp, year_label)
  }

  list(rf_raw = out_raw, rf_summary = out_summary, plots = plots)
}

get_rf_years <- function(dat) {
  years <- sort(unique(na.omit(dat$year)))

  if (!is.null(cfg$rf_years)) {
    years <- intersect(cfg$rf_years, years)
  }

  if (!isTRUE(cfg$rf_include_pretreatment)) {
    years <- setdiff(years, "2021")
  }

  years
}

run_rf_supplement <- function(dat, main_results = NULL) {
  years <- get_rf_years(dat)
  if (length(years) == 0) stop("No years selected for RF analysis. Check cfg$rf_years and cfg$rf_include_pretreatment.")

  rf_raw_all <- data.frame()
  rf_summary_all <- data.frame()
  rf_plots <- list()

  for (yr in years) {
    message("Running random forest supplement for ", yr)
    df <- prep_year_df(dat[dat$year == yr, , drop = FALSE])

    expected <- list()
    for (resp in responses) {
      message("  computing RF expected ES for ", resp)
      null_dist <- rf_null_distribution_rep(
        response = resp,
        data = df,
        sub_n_perm = cfg$n_rf_expected_boot_per_combo
      )
      expected[[resp]] <- rf_expected_effects_for_rows(null_dist, df)
    }

    rf_out <- run_rf_year(df, expected, year_label = yr)

    rf_raw_all <- rbind(rf_raw_all, rf_out$rf_raw)
    rf_summary_all <- rbind(rf_summary_all, rf_out$rf_summary)
    rf_plots[[yr]] <- rf_out$plots
  }

  plot_list <- list()
  for (yr in years) {
    for (resp in responses) {
      plot_list[[paste(yr, resp, sep = "_")]] <- rf_plots[[yr]][[resp]] +
        ggtitle(paste(yr, response_label(resp)))
    }
  }

  fig_s9 <- wrap_plots(plot_list, ncol = length(responses)) +
    plot_annotation(tag_levels = "a")

  if (cfg$save_outputs) {
    ggsave(file.path(cfg$output_dir, "figures", "FigureS9_random_forest_R2_corrected.pdf"), fig_s9, width = 14, height = 12)
    ggsave(file.path(cfg$output_dir, "figures", "FigureS9_random_forest_R2_corrected.tiff"), fig_s9, width = 14, height = 12, dpi = 350)
    write.csv(rf_raw_all, file.path(cfg$output_dir, "tables", "rf_r2_all_years_corrected.csv"), row.names = FALSE)
    write.csv(rf_summary_all, file.path(cfg$output_dir, "tables", "rf_r2_summary_all_years_corrected.csv"), row.names = FALSE)
    save(
      rf_plots,
      rf_raw_all,
      rf_summary_all,
      fig_s9,
      file = file.path(cfg$output_dir, "objects", "random_forest_R2_results_corrected.RData")
    )
  }

  list(figureS9 = fig_s9, rf_raw_all = rf_raw_all, rf_summary_all = rf_summary_all, rf_plots = rf_plots)
}



###############################################
# 11. Supplementary spline factor-response    #
###############################################

get_years_for_supplement <- function(dat, include_pretreatment = TRUE, years = NULL) {
  yy <- sort(unique(na.omit(dat$year)))
  if (!is.null(years)) yy <- intersect(years, yy)
  if (!isTRUE(include_pretreatment)) yy <- setdiff(yy, "2021")
  yy
}

prepare_factor_response_data <- function(dat, include_pretreatment = TRUE, years = NULL) {
  yy <- get_years_for_supplement(dat, include_pretreatment = include_pretreatment, years = years)
  base <- dat[dat$year %in% yy, , drop = FALSE]
  base$Lv_numeric <- safe_numeric(base$Lv)

  out <- list()
  k <- 0
  for (resp in responses) {
    tmp <- base[, c("year", "remark", "Lv", "Lv_numeric", resp), drop = FALSE]
    names(tmp)[names(tmp) == resp] <- "value"
    tmp$value <- safe_numeric(tmp$value)
    tmp <- tmp[is.finite(tmp$value) & is.finite(tmp$Lv_numeric), , drop = FALSE]
    if (nrow(tmp) == 0) next
    tmp$response <- resp
    tmp$response_label <- response_label(resp)
    k <- k + 1
    out[[k]] <- tmp
  }
  dplyr::bind_rows(out)
}

spline_stats_one <- function(df, response_name, year_label) {
  dd <- df %>%
    dplyr::filter(response == response_name, year == year_label, is.finite(value), is.finite(Lv_numeric))

  if (nrow(dd) < 5 || length(unique(dd$Lv_numeric)) < 3 || length(unique(dd$value)) < 2) {
    return(data.frame(
      response = response_name,
      response_label = response_label(response_name),
      year = year_label,
      n = nrow(dd),
      n_factor_levels = length(unique(dd$Lv_numeric)),
      spline_df = NA_real_,
      AIC_null = NA_real_,
      AIC_linear = NA_real_,
      AIC_spline = NA_real_,
      delta_AIC_spline_vs_linear = NA_real_,
      R2_linear = NA_real_,
      R2_spline = NA_real_,
      adj_R2_linear = NA_real_,
      adj_R2_spline = NA_real_,
      p_linear_slope = NA_real_,
      p_spline_model = NA_real_,
      p_spline_vs_linear = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  df_spline <- min(cfg$spline_df, max(1, length(unique(dd$Lv_numeric)) - 1))

  fit_null <- lm(value ~ 1, data = dd)
  fit_linear <- lm(value ~ Lv_numeric, data = dd)
  fit_spline <- lm(value ~ splines::ns(Lv_numeric, df = df_spline), data = dd)

  sm_linear <- summary(fit_linear)
  sm_spline <- summary(fit_spline)

  p_linear <- tryCatch(sm_linear$coefficients["Lv_numeric", "Pr(>|t|)"], error = function(e) NA_real_)

  p_spline_model <- tryCatch(
    stats::anova(fit_null, fit_spline)$`Pr(>F)`[2],
    error = function(e) NA_real_
  )

  p_spline_vs_linear <- tryCatch(
    stats::anova(fit_linear, fit_spline)$`Pr(>F)`[2],
    error = function(e) NA_real_
  )

  data.frame(
    response = response_name,
    response_label = response_label(response_name),
    year = year_label,
    n = nrow(dd),
    n_factor_levels = length(unique(dd$Lv_numeric)),
    spline_df = df_spline,
    AIC_null = AIC(fit_null),
    AIC_linear = AIC(fit_linear),
    AIC_spline = AIC(fit_spline),
    delta_AIC_spline_vs_linear = AIC(fit_spline) - AIC(fit_linear),
    R2_linear = sm_linear$r.squared,
    R2_spline = sm_spline$r.squared,
    adj_R2_linear = sm_linear$adj.r.squared,
    adj_R2_spline = sm_spline$adj.r.squared,
    p_linear_slope = p_linear,
    p_spline_model = p_spline_model,
    p_spline_vs_linear = p_spline_vs_linear,
    stringsAsFactors = FALSE
  )
}

plot_spline_factor_response <- function(plot_data) {
  if (nrow(plot_data) == 0) {
    return(ggplot() + theme_void() + labs(title = "No spline data available"))
  }

  plot_data$response_label <- factor(
    plot_data$response_label,
    levels = unique(vapply(responses, response_label, character(1)))
  )

  ggplot(plot_data, aes(x = Lv_numeric, y = value)) +
    geom_point(
      aes(color = factor(Lv_numeric)),
      alpha = 0.35,
      size = 1.3,
      position = position_jitter(width = 0.10, height = 0)
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ splines::ns(x, df = cfg$spline_df),
      se = TRUE,
      color = "black",
      fill = "grey70",
      linewidth = 0.7,
      alpha = 0.35
    ) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 8, 10)) +
    scale_color_manual(values = lv_cols, guide = "none") +
    facet_grid(response_label ~ year, scales = "free_y") +
    labs(
      x = "Number of factors",
      y = "Response value",
      title = "Relationships between factor number and selected response variables",
      subtitle = paste0("Spline regression with df = ", cfg$spline_df, "; grey bands are 95% confidence intervals")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
}

run_spline_factor_response_supplement <- function(dat) {
  plot_data <- prepare_factor_response_data(
    dat,
    include_pretreatment = cfg$spline_include_pretreatment,
    years = NULL
  )

  stat_rows <- list()
  k <- 0
  for (resp in responses) {
    for (yr in sort(unique(plot_data$year))) {
      k <- k + 1
      stat_rows[[k]] <- spline_stats_one(plot_data, resp, yr)
    }
  }

  stats_tbl <- dplyr::bind_rows(stat_rows)
  fig_spline <- plot_spline_factor_response(plot_data)

  if (cfg$save_outputs) {
    ggsave(
      file.path(cfg$output_dir, "figures", "FigureS_factor_number_response_spline.pdf"),
      fig_spline,
      width = max(10, 2.8 * length(unique(plot_data$year))),
      height = max(4, 2.6 * length(responses))
    )
    ggsave(
      file.path(cfg$output_dir, "figures", "FigureS_factor_number_response_spline.tiff"),
      fig_spline,
      width = max(10, 2.8 * length(unique(plot_data$year))),
      height = max(4, 2.6 * length(responses)),
      dpi = 350
    )
    write.csv(plot_data, file.path(cfg$output_dir, "tables", "factor_number_response_spline_plot_data.csv"), row.names = FALSE)
    write.csv(stats_tbl, file.path(cfg$output_dir, "tables", "factor_number_response_spline_statistics.csv"), row.names = FALSE)
    save(plot_data, stats_tbl, fig_spline, file = file.path(cfg$output_dir, "objects", "factor_number_response_spline_results.RData"))
  }

  list(figure = fig_spline, plot_data = plot_data, statistics = stats_tbl)
}

############################################################
# 12. Linear model comparison: identity vs identity + Lv    #
############################################################

lm_r2 <- function(fit) {
  out <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
  as.numeric(out)
}

lm_adj_r2 <- function(fit) {
  out <- tryCatch(summary(fit)$adj.r.squared, error = function(e) NA_real_)
  as.numeric(out)
}

fit_lm_factor_number_comparison_one <- function(df, response_name, year_label) {
  null_dist <- rf_null_distribution_rep(
    response = response_name,
    data = df,
    sub_n_perm = cfg$lm_comparison_expected_boot_per_combo
  )

  expected_tbl <- rf_expected_effects_for_rows(null_dist, df)
  es_col <- expected_effect_column_by_response[[response_name]]
  model_name <- final_null_model_by_response[[response_name]]

  d_lm <- df[df$remark %in% lv_levels, , drop = FALSE] %>%
    dplyr::left_join(expected_tbl, by = ".__rowid")

  d_lm$E_RF <- safe_numeric(d_lm[[es_col]])
  d_lm$Lv <- safe_numeric(d_lm$Lv)
  d_lm[[response_name]] <- safe_numeric(d_lm[[response_name]])

  d_lm <- d_lm[
    is.finite(d_lm[[response_name]]) & is.finite(d_lm$E_RF) & is.finite(d_lm$Lv),
    ,
    drop = FALSE
  ]

  if (nrow(d_lm) < 5 || length(unique(d_lm[[response_name]])) < 2 || length(unique(d_lm$Lv)) < 2) {
    long <- data.frame(
      Response = response_name,
      Response_label = response_label(response_name),
      Year = year_label,
      Model = c("Factor_identity_effect", "Factor_identity_effect+Factor_number"),
      AIC = NA_real_,
      DeltaAIC = NA_real_,
      R2 = NA_real_,
      Adj_R2 = NA_real_,
      ANOVA_F_added_factor_number = NA_real_,
      ANOVA_p_added_factor_number = NA_real_,
      Factor_identity_effect = paste(model_name, "model"),
      n = nrow(d_lm),
      stringsAsFactors = FALSE
    )
    pair <- data.frame(
      Response = response_name,
      Response_label = response_label(response_name),
      Year = year_label,
      Factor_identity_effect = paste(model_name, "model"),
      n = nrow(d_lm),
      AIC_identity = NA_real_,
      AIC_identity_factor_number = NA_real_,
      DeltaAIC_identity = NA_real_,
      DeltaAIC_identity_factor_number = NA_real_,
      Best_model = NA_character_,
      R2_identity = NA_real_,
      R2_identity_factor_number = NA_real_,
      Adj_R2_identity = NA_real_,
      Adj_R2_identity_factor_number = NA_real_,
      ANOVA_F_added_factor_number = NA_real_,
      ANOVA_p_added_factor_number = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(long = long, pair = pair, model_data = d_lm))
  }

  f_identity <- make_formula(response_name, "E_RF")
  f_plus <- make_formula(response_name, c("E_RF", "Lv"))

  fit_identity <- lm(f_identity, data = d_lm)
  fit_plus <- lm(f_plus, data = d_lm)

  aic_vals <- c(
    Factor_identity_effect = AIC(fit_identity),
    `Factor_identity_effect+Factor_number` = AIC(fit_plus)
  )
  delta_vals <- aic_vals - min(aic_vals, na.rm = TRUE)

  anova_tbl <- tryCatch(stats::anova(fit_identity, fit_plus), error = function(e) NULL)
  anova_F <- NA_real_
  anova_p <- NA_real_
  if (!is.null(anova_tbl) && nrow(anova_tbl) >= 2) {
    anova_F <- anova_tbl$F[2]
    anova_p <- anova_tbl$`Pr(>F)`[2]
  }

  long <- data.frame(
    Response = response_name,
    Response_label = response_label(response_name),
    Year = year_label,
    Model = names(aic_vals),
    AIC = as.numeric(aic_vals),
    DeltaAIC = as.numeric(delta_vals),
    R2 = c(lm_r2(fit_identity), lm_r2(fit_plus)),
    Adj_R2 = c(lm_adj_r2(fit_identity), lm_adj_r2(fit_plus)),
    ANOVA_F_added_factor_number = anova_F,
    ANOVA_p_added_factor_number = anova_p,
    Factor_identity_effect = paste(model_name, "model"),
    n = nrow(d_lm),
    stringsAsFactors = FALSE
  )

  pair <- data.frame(
    Response = response_name,
    Response_label = response_label(response_name),
    Year = year_label,
    Factor_identity_effect = paste(model_name, "model"),
    n = nrow(d_lm),
    AIC_identity = unname(aic_vals["Factor_identity_effect"]),
    AIC_identity_factor_number = unname(aic_vals["Factor_identity_effect+Factor_number"]),
    DeltaAIC_identity = unname(delta_vals["Factor_identity_effect"]),
    DeltaAIC_identity_factor_number = unname(delta_vals["Factor_identity_effect+Factor_number"]),
    Best_model = names(aic_vals)[which.min(aic_vals)],
    R2_identity = lm_r2(fit_identity),
    R2_identity_factor_number = lm_r2(fit_plus),
    Adj_R2_identity = lm_adj_r2(fit_identity),
    Adj_R2_identity_factor_number = lm_adj_r2(fit_plus),
    ANOVA_F_added_factor_number = anova_F,
    ANOVA_p_added_factor_number = anova_p,
    stringsAsFactors = FALSE
  )

  list(long = long, pair = pair, model_data = d_lm)
}

run_lm_factor_number_comparison <- function(dat) {
  years <- get_years_for_supplement(
    dat,
    include_pretreatment = cfg$lm_comparison_include_pretreatment,
    years = NULL
  )

  long_rows <- list()
  pair_rows <- list()
  data_rows <- list()
  k <- 0

  for (yr in years) {
    message("Running LM factor-number comparison for ", yr)
    df <- prep_year_df(dat[dat$year == yr, , drop = FALSE])

    for (resp in responses) {
      k <- k + 1
      out <- fit_lm_factor_number_comparison_one(df, resp, yr)
      long_rows[[k]] <- out$long
      pair_rows[[k]] <- out$pair
      data_rows[[k]] <- data.frame(
        year = yr,
        response = resp,
        response_label = response_label(resp),
        out$model_data,
        stringsAsFactors = FALSE
      )
    }
  }

  long_tbl <- dplyr::bind_rows(long_rows) %>%
    dplyr::arrange(Response, Year, Model)

  pair_tbl <- dplyr::bind_rows(pair_rows) %>%
    dplyr::arrange(Response, Year)

  model_data <- dplyr::bind_rows(data_rows)

  if (cfg$save_outputs) {
    write.csv(
      long_tbl,
      file.path(cfg$output_dir, "tables", "Supplementary_Data5_LM_factor_number_comparison_long.csv"),
      row.names = FALSE
    )
    write.csv(
      pair_tbl,
      file.path(cfg$output_dir, "tables", "Supplementary_Data5_LM_factor_number_comparison_pairs.csv"),
      row.names = FALSE
    )
    write.csv(
      model_data,
      file.path(cfg$output_dir, "tables", "Supplementary_Data5_LM_model_data.csv"),
      row.names = FALSE
    )
    save(long_tbl, pair_tbl, model_data, file = file.path(cfg$output_dir, "objects", "lm_factor_number_comparison_results.RData"))
  }

  list(long = long_tbl, pair = pair_tbl, model_data = model_data)
}

##############################
# 11. Run selected analyses  #
##############################

df_field <- read_field_data(cfg$data_path, cfg$drop_first_column)

main_results <- NULL
if (cfg$run_main_figure) {
  main_results <- run_main_figure2(df_field)
  print(main_results$figure2)
}

if (cfg$run_supp_raw_ridge && !is.null(main_results)) {
  raw_supp <- save_raw_ridge_supplements(main_results$raw_plots)
}

if (cfg$run_supp_nmds) {
  nmds_results <- run_nmds_supplement(df_field)
}

if (cfg$run_supp_deviation) {
  deviation_results <- run_deviation_supplement(df_field)
}

if (cfg$run_supp_random_forest) {
  rf_results <- run_rf_supplement(df_field, main_results)
}

if (cfg$run_supp_spline_factor_response) {
  spline_results <- run_spline_factor_response_supplement(df_field)
  print(spline_results$figure)
}

if (cfg$run_lm_factor_number_comparison) {
  lm_factor_number_results <- run_lm_factor_number_comparison(df_field)
}

message("Field experiment analysis complete. Outputs saved to: ", cfg$output_dir)
