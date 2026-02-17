# Declare global variables used in NSE (dplyr/ggplot2) to avoid R CMD check notes
utils::globalVariables(c(

  # Data columns
  "province", "province_std", "province_idx",
  "year", "year_idx", "year_factor",
  "age_group", "age_group_idx", "age_start", "age_end", "age_mid", "age_weight",
  "sex", "sex_std",
  "n", "row_id", "is_outlier", "is_valid", "study",

  # DMFT columns
  "mean_dmft", "mean_DMFT", "sd_dmft", "sd_DMFT", "se_dmft", "se_DMFT",
  "mean_d", "mean_m", "mean_f", "mean_D", "mean_M", "mean_F",
  "lci_dmft", "uci_dmft", "lci_DMFT", "uci_DMFT",
  "se_imputed", "sd_imputed", "se_source", "n_effective",
  "se_dmft_imputed", "se_DMFT_imputed", "precision",

  # Model columns
  "y", "region_std", "location",

  # AST columns
  "residual_AST", "residual_ast", "predicted_base",

  # Prediction columns
  "predicted", "lower", "upper",

  # Projection columns
  "scenario", "mean_proj",

  # Visualization
  "value", "PRENAME", "region",

  # Diagnostic
  "residual", "mean_resid"
))
