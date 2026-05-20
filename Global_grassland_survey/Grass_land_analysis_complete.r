############################################################
## 0. Packages
############################################################

library(dplyr)
library(ggplot2)
library(simpleboot)
library(MuMIn)
library(relaimpo)
library(car)
library(gvlma)

############################################################
## 1. Import dataset
############################################################

df <- read.csv("data.csv")[, -1]

############################################################
## 2. Filter grassland subset
############################################################

dff <- df %>%
  dplyr::filter(Biome %in% c("Arid", "Temperate", "Continental")) %>%
  dplyr::filter(
    NDVI_2001_2020_Normalized >= 0.2,
    NDVI_2001_2020_Normalized <= 0.65
  ) %>%
  dplyr::filter(
    Latitude <= 55,
    Latitude >= -55,
    Latitude >= 20 | Latitude <= -20
  )

cat("Number of grassland subset rows:", nrow(dff), "\n")

############################################################
## 3. Define response variable
############################################################

#resp_var <- "NDVI_2001_2020_Normalized"
 resp_var <- "Plant_cover_v3"

if (!resp_var %in% names(dff)) {
  stop(paste0("Response variable not found: ", resp_var))
}

############################################################
## 4. Define stressor columns
############################################################
## Change these if your column positions differ

# historical records 1970_2000
# stressor_cols <- names(dff)[19:24]

# soil sampling date 2016_2019
stressor_cols <- names(dff)[24:29]

# climate anomalies
# stressor_cols <- names(dff)[31:36]

cat("Stressors used:\n")
print(stressor_cols)

MS_raw <- dff[, stressor_cols, drop = FALSE]

############################################################
## 5. Scale stressors within grassland subset to 0-100
############################################################

scale_01_to_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(NA_real_, length(x)))
  }
  100 * (x - rng[1]) / (rng[2] - rng[1])
}

MS <- as.data.frame(lapply(MS_raw, scale_01_to_100))

############################################################
## 6. Threshold analysis
############################################################

thres <- seq(5, 95, 1)
RES  <- list()
RES2 <- list()

for (i in seq_along(thres)) {
  
  dfi <- data.frame(
    MF = dff[[resp_var]],
    nMS = apply(MS, 1, function(rr, tt) sum(rr >= tt, na.rm = TRUE), tt = thres[i]),
    threshold = thres[i]
  )
  
  dfi <- dfi[complete.cases(dfi), ]
  
  # Skip thresholds with no variation in stressor counts
  if (nrow(dfi) < 10) next
  if (length(unique(dfi$nMS)) < 2) next
  
  RES[[length(RES) + 1]] <- dfi
  
  mdl <- lm(MF ~ nMS, data = dfi)
  resi <- simpleboot::lm.boot(mdl, 100)
  
  boot_coefs <- unlist(lapply(resi$boot.list, function(mm) mm$coef[2]))
  boot_r2    <- unlist(lapply(resi$boot.list, function(mm) {
    if (!is.null(mm$rsquare)) mm$rsquare else NA_real_
  }))
  
  RES2[[length(RES2) + 1]] <- data.frame(
    threshold = thres[i],
    coefs = boot_coefs,
    R2 = boot_r2
  )
}

RES  <- do.call(rbind, RES)
RES2 <- do.call(rbind, RES2)

############################################################
## 7. Summarise threshold results
############################################################

mytable <- RES2 %>%
  dplyr::group_by(threshold) %>%
  dplyr::summarise(
    n_boot = sum(is.finite(coefs)),
    Average.R2   = mean(R2, na.rm = TRUE),
    sd.R2        = sd(R2, na.rm = TRUE),
    Average.coef = mean(coefs, na.rm = TRUE),
    ci.025       = quantile(coefs, 0.005, na.rm = TRUE),
    ci.975       = quantile(coefs, 0.995, na.rm = TRUE),
    pval = {
      cc <- coefs[is.finite(coefs)]
      if (length(cc) == 0) NA_real_
      else if (mean(cc < 0) > 0.5) sum(cc < 0) else sum(cc > 0)
    },
    .groups = "drop"
  )

############################################################
## 8. Plot threshold analysis
############################################################

p_threshold <- ggplot(mytable, aes(x = threshold, y = Average.coef)) +
  theme_classic() +
  geom_line() +
  geom_ribbon(aes(ymin = ci.025, ymax = ci.975), fill = "grey", alpha = 0.4) +
  geom_hline(yintercept = 0, color = "red") +
  xlab("") +
  ylab(paste("Effect of number of stressors\non", resp_var)) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.title.y = element_text(
      margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0),
      size = 14
    ),
    axis.title.x = element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    ),
    axis.ticks.length = grid::unit(0.2, "cm")
  )

print(p_threshold)

############################################################
## 8B. Plot number of stressors above threshold vs response
############################################################

## RES already contains:
## MF        = response variable
## nMS       = number of stressors above each threshold
## threshold = threshold value

RES_plot <- RES %>%
  dplyr::filter(
    is.finite(MF),
    is.finite(nMS),
    is.finite(threshold)
  )

p_stressor_lines <- ggplot(
  data = RES_plot,
  aes(
    x = nMS,
    y = MF,
    group = threshold,
    color = threshold
  )
) +
  theme_classic() +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.35
  ) +
  scale_color_viridis_c(
    option = "viridis",
    direction = -1,
    name = "Threshold"
  ) +
  xlab("\nNumber of stressors\nabove threshold") +
  ylab(paste0(resp_var, "\n")) +
  ggtitle("Grassland subset") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_stressor_lines)

ggsave(
  filename = paste0("stressor_threshold_lines_", resp_var, ".tiff"),
  plot = p_stressor_lines,
  device = "tiff",
  units = "cm",
  width = 8,
  height = 6,
  dpi = 350
)

############################################################
## 9. Build MS25 / MS50 / MS75 for multi-model comparison
############################################################

dff$MS25 <- apply(MS, 1, function(rr) sum(rr >= 25, na.rm = TRUE))
dff$MS50 <- apply(MS, 1, function(rr) sum(rr >= 50, na.rm = TRUE))
dff$MS75 <- apply(MS, 1, function(rr) sum(rr >= 75, na.rm = TRUE))

############################################################
## 10. Build modelling dataframe
############################################################

dfi_mod <- dff[, c(resp_var, "MS25", "MS50", "MS75"), drop = FALSE]
dfi_mod <- cbind(dfi_mod, MS)

dfi_mod <- dfi_mod[complete.cases(dfi_mod), ]

cat("Rows used in multi-model analysis:", nrow(dfi_mod), "\n")

############################################################
## 11. Global linear model
############################################################

## Explicit predictor list -- avoid "~ ." for MuMIn::dredge()
predictor_vars <- setdiff(names(dfi_mod), resp_var)

form_global <- reformulate(
  termlabels = predictor_vars,
  response = resp_var
)

## Important:
## Fit the model using eval/substitute so the stored lm call contains
## the actual formula, not just the symbol "form_global".
mdl_global <- eval(substitute(
  lm(FORMULA, data = dfi_mod, na.action = na.fail),
  list(FORMULA = form_global)
))

############################################################
## 12. Diagnostics
############################################################

assumptions <- gvlma::gvlma(mdl_global)
vif_res <- car::vif(mdl_global)

cat("\n===== Global model summary =====\n")
print(summary(mdl_global))

cat("\n===== gvlma diagnostics =====\n")
print(summary(assumptions))

cat("\n===== VIF =====\n")
print(vif_res)

############################################################
## 13. Relative importance analysis
############################################################

RIMP <- relaimpo::calc.relimp(mdl_global, type = "lmg")

RIMP_df <- data.frame(
  variable = names(RIMP@lmg),
  lmg.imp = RIMP@lmg,
  perc.imp = RIMP@lmg / RIMP@R2
)

cat("\n===== Relative importance =====\n")
print(RIMP_df[order(-RIMP_df$perc.imp), ])

############################################################
## 14. MuMIn multi-model comparison
############################################################

## This is now safe because mdl_global has an explicit formula
dr1 <- MuMIn::dredge(
  mdl_global,
  rank = "BIC",
  extra = "R^2"
)

dr2 <- subset(dr1, delta < 4)

cat("\n===== Top dredge models (all models) =====\n")
print(head(dr1, 20))

cat("\n===== Best-supported models: delta < 4 =====\n")
print(dr2)

############################################################
## 15. Model averaging over best models
############################################################

if (nrow(dr2) > 1) {
  
  avg_mod <- MuMIn::model.avg(dr2)
  
  cat("\n===== Model-averaged summary =====\n")
  print(summary(avg_mod))
  
} else if (nrow(dr2) == 1) {
  
  cat("\n===== Only one best-supported model: no model averaging performed =====\n")
  print(dr2)
  
  best_mod <- MuMIn::get.models(dr1, subset = 1)[[1]]
  
  cat("\n===== Summary of the single best model =====\n")
  print(summary(best_mod))
  
} else {
  
  cat("\n===== No models found with delta < 4 =====\n")
}

############################################################
## 16. Optional stricter MuMIn analysis:
## keep all individual stressors fixed, only let MS25/MS50/MS75 vary
############################################################

fixed_stressors <- intersect(stressor_cols, names(dfi_mod))

dr1_fixed <- MuMIn::dredge(
  mdl_global,
  rank = "BIC",
  fixed = fixed_stressors,
  extra = "R^2"
)

dr2_fixed <- subset(dr1_fixed, delta < 4)

cat("\n===== Fixed-stressor dredge: delta < 4 =====\n")
print(dr2_fixed)

if (nrow(dr2_fixed) > 1) {
  
  avg_mod_fixed <- MuMIn::model.avg(dr2_fixed)
  
  cat("\n===== Fixed-stressor model-averaged summary =====\n")
  print(summary(avg_mod_fixed))
  
} else if (nrow(dr2_fixed) == 1) {
  
  cat("\n===== Only one fixed-stressor best-supported model: no model averaging performed =====\n")
  print(dr2_fixed)
  
  best_mod_fixed <- MuMIn::get.models(dr1_fixed, subset = 1)[[1]]
  
  cat("\n===== Summary of the single best fixed-stressor model =====\n")
  print(summary(best_mod_fixed))
  
} else {
  
  cat("\n===== No fixed-stressor models found with delta < 4 =====\n")
}

############################################################
## 17. Save outputs
############################################################

write.csv(mytable, "threshold_summary_NDVI.csv", row.names = FALSE)
write.csv(RES, "threshold_all_dataframes_NDVI.csv", row.names = FALSE)
write.csv(RES2, "threshold_bootstraps_NDVI.csv", row.names = FALSE)

write.csv(RIMP_df, "relative_importance_NDVI.csv", row.names = FALSE)
write.csv(as.data.frame(dr1), "MuMIn_all_models_NDVI.csv", row.names = FALSE)
write.csv(as.data.frame(dr2), "MuMIn_delta_lt_4_NDVI.csv", row.names = FALSE)
write.csv(as.data.frame(dr1_fixed), "MuMIn_fixed_all_models_NDVI.csv", row.names = FALSE)
write.csv(as.data.frame(dr2_fixed), "MuMIn_fixed_delta_lt_4_NDVI.csv", row.names = FALSE)

############################################################
## 18. Optional simple variable-importance plot
############################################################

p_imp <- ggplot(
  RIMP_df %>% arrange(desc(perc.imp)),
  aes(x = reorder(variable, perc.imp), y = perc.imp)
) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("Relative importance (% of model R2)") +
  ggtitle(paste("Relative importance for", resp_var))

print(p_imp)