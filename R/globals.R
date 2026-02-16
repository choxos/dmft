# Declare global variables used in NSE (dplyr/ggplot2) to avoid R CMD check notes
utils::globalVariables(c(

  # Data columns
  "province", "province_std", "province_idx", "province_idx2",
  "year", "year_idx", "year_idx2",
  "age_group", "age_group_idx", "age_start", "age_end", "age_mid", "age_weight",
  "sex", "sex_std", "sex_idx",
  "n", "row_id", "is_outlier", "is_valid", "study",

  # DMFT columns
  "mean_dmft", "mean_DMFT", "sd_dmft", "sd_DMFT", "se_dmft", "se_DMFT",
  "mean_d", "mean_m", "mean_f", "mean_D", "mean_M", "mean_F",
  "lci_dmft", "uci_dmft", "lci_DMFT", "uci_DMFT",
  "se_imputed", "sd_imputed", "se_source", "n_effective",
  "se_dmft_imputed", "se_DMFT_imputed", "precision",

  # Model columns
  "y", "offset_log_n", "fold",
  "st_idx", "sa_idx", "ta_idx",

  # Prediction columns
  "eta", "eta_sd", "predicted_dmft", "predicted_DMFT",
  "dmft_lower", "dmft_upper", "DMFT_lower", "DMFT_upper",

  # Projection columns
  "scenario", "mean_proj", "lower", "upper",

  # Visualization
  "value", "PRENAME", "region", "predicted",

  # Diagnostic
  "residual", "mean_resid"
))
