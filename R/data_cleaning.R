#' Clean and standardize DMFT data
#'
#' Performs column standardization, region name mapping, sex standardization,
#' age group assignment, outlier flagging, and value validation. Returns
#' separate deciduous and permanent datasets ready for modeling.
#'
#' @param data A data frame from [dmft_load()].
#' @param config A [dmft_config] object.
#'
#' @returns A list with elements `deciduous` and `permanent`, each a tibble
#'   with standardized columns and imputed uncertainty.
#' @export
dmft_clean <- function(data, config) {
  stopifnot(inherits(config, "dmft_config"))

  # --- Harmonize DMFT column names BEFORE lowercasing ---
  # This preserves the deciduous (dmft) vs permanent (DMFT) distinction
  data <- harmonize_dmft_columns(data)

  # --- Normalize other column names ---
  # Protect DMFT-specific columns from case-folding
  dmft_protected <- c("mean_dmft", "mean_DMFT", "sd_dmft", "sd_DMFT",
                       "se_dmft", "se_DMFT", "lci_dmft", "lci_DMFT",
                       "uci_dmft", "uci_DMFT")
  is_protected <- names(data) %in% dmft_protected
  names(data)[!is_protected] <- tolower(names(data)[!is_protected])
  names(data)[!is_protected] <- gsub("[^a-z0-9_]", "_", names(data)[!is_protected])

  # --- Region standardization ---
  rcol <- config$region_col
  if (!rcol %in% names(data)) {
    cli_abort("Region column '{rcol}' not found in data.")
  }
  data$region_std <- standardize_regions(data[[rcol]], config)

  # --- Sex standardization ---
  if ("sex" %in% names(data)) {
    data$sex_std <- standardize_sex(data$sex)
  } else {
    data$sex_std <- "Both"
  }

  # --- Numeric coercion ---
  for (col in c("year", "age_start", "age_end", "n")) {
    if (col %in% names(data)) data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  }

  if (all(c("age_start", "age_end") %in% names(data))) {
    data$age_mid <- (data$age_start + data$age_end) / 2
  }

  # --- Split by dentition ---
  deciduous <- data[!is.na(data$mean_dmft) &
                      (is.na(data$age_end) | data$age_end <= 14), ]
  permanent <- data[!is.na(data$mean_DMFT) &
                      (is.na(data$age_start) | data$age_start >= 5), ]

  # --- Age groups ---
  if (nrow(deciduous) > 0) {
    deciduous <- add_age_groups(deciduous, "deciduous", config)
    deciduous <- flag_outliers(deciduous, "mean_dmft")
    deciduous <- validate_dmft(deciduous, "deciduous", config)
    deciduous <- impute_uncertainty(deciduous, "mean_dmft",
                                    sd_col = find_col(deciduous, "sd_dmft"),
                                    se_col = find_col(deciduous, "se_dmft"),
                                    cv_default = config$cv_default)
  }
  if (nrow(permanent) > 0) {
    permanent <- add_age_groups(permanent, "permanent", config)
    permanent <- flag_outliers(permanent, "mean_DMFT")
    permanent <- validate_dmft(permanent, "permanent", config)
    permanent <- impute_uncertainty(permanent, "mean_DMFT",
                                    sd_col = find_col(permanent, "sd_DMFT"),
                                    se_col = find_col(permanent, "se_DMFT"),
                                    cv_default = config$cv_default)
  }

  cli_alert_success("Cleaned: {nrow(deciduous)} deciduous, {nrow(permanent)} permanent records")

  list(deciduous = dplyr::as_tibble(deciduous),
       permanent = dplyr::as_tibble(permanent))
}


# -- internal helpers ----------------------------------------------------------

#' Standardize region names using config mapping
#' @keywords internal
standardize_regions <- function(region_col, config) {
  raw <- tolower(trimws(region_col))
  mapping <- config$region_mapping
  canonical_lower <- setNames(config$regions, tolower(config$regions))

  vapply(raw, function(r) {
    if (is.na(r) || r == "") return(NA_character_)
    # Exact canonical match
    if (r %in% names(canonical_lower)) return(canonical_lower[[r]])
    # User mapping
    if (!is.null(mapping) && r %in% tolower(names(mapping))) {
      key <- names(mapping)[tolower(names(mapping)) == r][1]
      return(mapping[[key]])
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}


#' Standardize sex column
#' @keywords internal
standardize_sex <- function(sex_col) {
  s <- tolower(trimws(sex_col))
  dplyr::case_when(
    s %in% c("male", "m", "1", "men", "boy", "boys") ~ "Male",
    s %in% c("female", "f", "0", "2", "women", "girl", "girls") ~ "Female",
    s %in% c("both", "all", "combined", "total", "mixed") ~ "Both",
    TRUE ~ NA_character_
  )
}


#' Harmonize DMFT column names
#' @keywords internal
harmonize_dmft_columns <- function(data) {
  # Mean
  if (!"mean_dmft" %in% names(data)) {
    col <- grep("^mean_dmft", names(data), value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(col)) data$mean_dmft <- suppressWarnings(as.numeric(data[[col]]))
    else data$mean_dmft <- NA_real_
  } else {
    data$mean_dmft <- suppressWarnings(as.numeric(data$mean_dmft))
  }

  if (!"mean_DMFT" %in% names(data)) {
    col <- grep("^mean_dmft[0-9]", names(data), value = TRUE, ignore.case = FALSE)[1]
    if (is.na(col)) col <- grep("^mean_DMFT", names(data), value = TRUE)[1]
    if (!is.na(col)) data$mean_DMFT <- suppressWarnings(as.numeric(data[[col]]))
    else data$mean_DMFT <- NA_real_
  } else {
    data$mean_DMFT <- suppressWarnings(as.numeric(data$mean_DMFT))
  }

  # SD / SE
  for (prefix in c("sd_dmft", "se_dmft", "sd_DMFT", "se_DMFT")) {
    if (!prefix %in% names(data)) {
      col <- grep(paste0("^", prefix), names(data), value = TRUE,
                  ignore.case = FALSE)[1]
      if (!is.na(col)) data[[prefix]] <- suppressWarnings(as.numeric(data[[col]]))
    }
  }

  data
}


#' Assign age group indices
#' @keywords internal
add_age_groups <- function(data, dentition, config) {
  if (!all(c("age_start", "age_end") %in% names(data))) {
    cli_alert_warning("age_start/age_end missing; setting age_group = NA")
    data$age_group <- NA_character_
    data$age_group_idx <- NA_integer_
    return(data)
  }

  if (dentition == "deciduous") {
    groups <- config$age_groups_deciduous
    lower  <- config$age_groups_deciduous_lower
    upper  <- config$age_groups_deciduous_upper
  } else {
    groups <- config$age_groups_permanent
    lower  <- config$age_groups_permanent_lower
    upper  <- config$age_groups_permanent_upper
  }

  mid <- (data$age_start + data$age_end) / 2
  ag <- rep(NA_character_, nrow(data))
  for (g in seq_along(groups)) {
    ag[!is.na(mid) & mid >= lower[g] & mid <= upper[g]] <- groups[g]
  }

  data$age_group <- ag
  data$age_group_idx <- match(ag, groups)
  data
}


#' Flag outliers via IQR
#' @keywords internal
flag_outliers <- function(data, value_col, multiplier = 3) {
  vals <- data[[value_col]]
  vals <- vals[!is.na(vals)]
  if (length(vals) < 4) { data$is_outlier <- FALSE; return(data) }
  q1 <- quantile(vals, 0.25)
  q3 <- quantile(vals, 0.75)
  iqr <- q3 - q1
  data$is_outlier <- !is.na(data[[value_col]]) &
    (data[[value_col]] < q1 - multiplier * iqr |
     data[[value_col]] > q3 + multiplier * iqr)
  data
}


#' Validate DMFT values
#' @keywords internal
validate_dmft <- function(data, dentition, config) {
  max_val <- if (dentition == "deciduous") config$dmft_max_deciduous
             else config$dmft_max_permanent
  vcol <- if (dentition == "deciduous") "mean_dmft" else "mean_DMFT"
  data$is_valid <- is.na(data[[vcol]]) |
    (data[[vcol]] >= 0 & data[[vcol]] <= max_val)
  data
}


#' Find a column name or return NULL
#' @keywords internal
find_col <- function(data, col) {
  if (col %in% names(data)) col else NULL
}
