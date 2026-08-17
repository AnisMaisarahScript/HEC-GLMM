# =============================================================================
# HUMAN–ELEPHANT CONFLICT (HEC) LAND-COVER TRANSITION WORKFLOW
# =============================================================================
#
# PURPOSE
# This script reproduces the two-stage analysis used to evaluate how specific
# land-cover (LC) transitions are associated with HEC intensity in Johor,
# Peninsular Malaysia.
#
# ANALYTICAL STRUCTURE
#
# Stage 1 — Gaussian GLMM
#   T1: LC transitions 2010–2015 aligned with HEC 2013–2017
#   T2: LC transitions 2015–2021 aligned with HEC 2018–2023
#
#   The T1 and T2 datasets are combined and fitted using:
#
#       HEC_kernel ~ LC transitions + Duration + (1 | Year)
#
#   Response:
#       KDE-derived HEC intensity (continuous)
#
#   Fixed effects:
#       Binary-coded LC transition variables
#       Duration (T1 vs T2)
#
#   Random effect:
#       Year
#
#   Error distribution:
#       Gaussian with identity link
#
#   Predictor screening:
#       Zero-variance predictors are removed.
#       Predictors with VIF > 5 are iteratively removed.
#
#   Significant LC transitions:
#       p < 0.05
#
#   Model performance:
#       Marginal R²  = variation explained by fixed effects
#       Conditional R² = variation explained by fixed + random effects
#
# Stage 2 — Cumulative impact + translation verification
#
#   For each statistically significant LC transition:
#
#       DeltaHEC = beta × n
#
#   where:
#       beta = GLMM coefficient for the transition
#       n    = observed HEC presence-cell count for that transition
#
#   The cumulative impact is translated into predicted HEC presence-cell counts
#   using:
#
#       ln(n_obs + 1) = alpha + lambda(DeltaHEC)
#
#   Predicted values are back-transformed to the original count scale:
#
#       n_pred = exp(predicted log count) - 1
#
#   Translation performance is evaluated using:
#       - slope p-value
#       - RMSE in presence-cell count units
#       - normalized RMSE (%) = RMSE / mean(observed counts) × 100
#
#   In this workflow, acceptable representation is defined as:
#       p < 0.05 AND normalized RMSE <= 25%
#
#
# SPATIAL RESOLUTION
# All raster-derived inputs use 250 m × 250 m cells.
# Each cell therefore represents 0.0625 km².
#
# =============================================================================
# REQUIRED PACKAGES
# =============================================================================
#
# Install once if needed:
# install.packages(c("glmmTMB", "car", "performance"))
#
# glmmTMB    : Gaussian generalized linear mixed-effects model
# car        : variance inflation factor (VIF)
# performance: marginal and conditional R²
#
# =============================================================================

library(glmmTMB)
library(car)
library(performance)

# =============================================================================
# 0) USER SETTINGS — EDIT THIS SECTION BEFORE RUNNING
# =============================================================================

# -----------------------------------------------------------------------------
# 0.1 WORKING DIRECTORY
# -----------------------------------------------------------------------------
# CHANGE this path to the folder containing your input CSV files.
# Use "/" or "\\" in Windows paths.

base_dir <- "C:/Users/HP-15FC/Desktop/Software/new R/Revised GLMM (Include ROC and PI)"

# -----------------------------------------------------------------------------
# 0.2 INPUT FILES
# -----------------------------------------------------------------------------

# DATA STRUCTURE
# Each row represents 250m x 250m raster cell.
# HEC_kernel contains the KDE-derived HEC intensity for that observation.
# LC transition variables are binary coded:
#     1 = the specified transition occurred
#     0 = the specified transition did not occur

# Required files:
#
# T1.csv
#   LC transitions 2010–2015 aligned with HEC 2013–2017.
#
# T2.csv
#   LC transitions 2015–2021 aligned with HEC 2018–2023.
#
# Presence_cells.csv
#   Must contain:
#       Predictor
#   and either:
#       T1 and T2
#   or:
#       Presence_Cells
#
# If T1 and T2 columns are supplied, the script calculates:
#       Presence_Cells = T1 + T2
#
# Predictor names must match the transition variables used in the GLMM.
# Presence_Cells represents the number of 250 × 250 m cells in which HEC
# presence was associated with each LC transition.

files <- list(
  T1 = file.path(base_dir, "T1.csv"),
  T2 = file.path(base_dir, "T2.csv")
)

presence_file <- file.path(base_dir, "Presence_cells.csv")

# -----------------------------------------------------------------------------
# 0.3 RESPONSE VARIABLE
# -----------------------------------------------------------------------------
# This must match the KDE-response column in T1.csv and T2.csv.

response_variable <- "HEC_kernel"

# -----------------------------------------------------------------------------
# 0.4 MODEL SETTINGS
# -----------------------------------------------------------------------------

vif_threshold <- 5
significance_threshold <- 0.05
rmse_threshold_percent <- 25

# -----------------------------------------------------------------------------
# 0.5 RASTER CELL AREA
# -----------------------------------------------------------------------------

cell_area_km2 <- 0.0625

# -----------------------------------------------------------------------------
# 0.6 LAND-COVER TRANSITION VARIABLES
# -----------------------------------------------------------------------------
#
# Codes:
#   F = Forest
#   P = Plantation
#   S = Settlement
#   I = Idle land
#   W = Water body
#
# "TO" indicates direction.
#
# Examples:
#   FTOP = Forest -> Plantation
#   PTOF = Plantation -> Forest
#   ITOP = Idle land -> Plantation
#
# Predictors must match the column names in T1.csv and T2.csv.

predictors <- c(
  "WTOF","WTOS","WTOI","WTOP",
  "STOF","STOW","STOI","STOP",
  "ITOF","ITOW","ITOS","ITOP",
  "PTOF","PTOW","PTOS","PTOI",
  "FTOP","FTOW","FTOS","FTOI"
)

predictors <- unique(predictors)

# =============================================================================
# 1) HELPER FUNCTIONS — NO EDITING NORMALLY REQUIRED
# =============================================================================

rename_response_if_needed <- function(df, response_variable) {

  possible_response_names <- c(
    response_variable,
    "Kernel13to17",
    "kernel13_17",
    "kernel18_23",
    "HEC"
  )

  found <- intersect(possible_response_names, names(df))

  if (!(response_variable %in% names(df)) && length(found) > 0) {
    names(df)[names(df) == found[1]] <- response_variable
  }

  df
}


read_timeline <- function(path, label, predictors, response_variable) {

  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }

  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- rename_response_if_needed(df, response_variable)

  # Duration is created from the file label.
  df$Duration <- label

  required_cols <- c(response_variable, "Year", "Duration")
  missing_required <- setdiff(required_cols, names(df))

  if (length(missing_required) > 0) {
    stop(
      "Missing required columns in ", basename(path), ": ",
      paste(missing_required, collapse = ", ")
    )
  }

  # Missing LC transition columns are treated as absence (0).
  missing_predictors <- setdiff(predictors, names(df))

  if (length(missing_predictors) > 0) {
    message(
      "Adding missing transition columns as 0 in ", label, ": ",
      paste(missing_predictors, collapse = ", ")
    )

    for (p in missing_predictors) {
      df[[p]] <- 0
    }
  }

  needed_cols <- c(response_variable, predictors, "Duration", "Year")
  df <- df[, needed_cols, drop = FALSE]

  # Ensure correct data types.
  for (p in predictors) {
    df[[p]] <- as.numeric(df[[p]])
  }

  df[[response_variable]] <- as.numeric(df[[response_variable]])
  df$Year <- as.factor(df$Year)
  df$Duration <- as.factor(df$Duration)

  df
}


remove_zero_variance_predictors <- function(df, predictors) {

  variances <- sapply(
    df[, predictors, drop = FALSE],
    var,
    na.rm = TRUE
  )

  zero_var <- names(
    variances[variances == 0 | is.na(variances)]
  )

  if (length(zero_var) > 0) {
    message(
      "Removed zero-variance predictors: ",
      paste(zero_var, collapse = ", ")
    )
  } else {
    message("No zero-variance predictors found.")
  }

  setdiff(predictors, zero_var)
}


iterative_vif_filter <- function(
    df,
    response_variable,
    predictors,
    threshold = 5,
    output_path = NULL
) {

  current_preds <- predictors
  vif_history <- data.frame()
  iteration <- 1

  if (length(current_preds) < 2) {
    return(
      list(
        predictors = current_preds,
        table = vif_history
      )
    )
  }

  repeat {

    # VIF is evaluated on the fixed-effect structure.
    vif_formula <- as.formula(
      paste(
        response_variable,
        "~",
        paste(c(current_preds, "Duration"), collapse = " + ")
      )
    )

    vif_fit <- lm(vif_formula, data = df)
    vif_values <- car::vif(vif_fit)

    if (is.matrix(vif_values)) {

      vif_table <- data.frame(
        Predictor = rownames(vif_values),
        GVIF = vif_values[, "GVIF"],
        Df = vif_values[, "Df"],
        VIF = vif_values[, "GVIF^(1/(2*Df))"],
        row.names = NULL
      )

    } else {

      vif_table <- data.frame(
        Predictor = names(vif_values),
        VIF = as.numeric(vif_values),
        row.names = NULL
      )
    }

    predictor_vifs <- vif_table[
      vif_table$Predictor %in% current_preds,
      ,
      drop = FALSE
    ]

    predictor_vifs <- predictor_vifs[
      order(-predictor_vifs$VIF),
      ,
      drop = FALSE
    ]

    vif_table$Iteration <- iteration
    vif_table$Decision <- "Kept"

    if (
      nrow(predictor_vifs) == 0 ||
      max(predictor_vifs$VIF, na.rm = TRUE) <= threshold
    ) {

      vif_history <- rbind(vif_history, vif_table)
      break
    }

    drop_pred <- predictor_vifs$Predictor[1]

    vif_table$Decision[
      vif_table$Predictor == drop_pred
    ] <- "Dropped"

    vif_history <- rbind(vif_history, vif_table)

    message(
      "Removed predictor because VIF > ",
      threshold,
      ": ",
      drop_pred
    )

    current_preds <- setdiff(current_preds, drop_pred)
    iteration <- iteration + 1

    if (length(current_preds) < 2) {
      break
    }
  }

  if (!is.null(output_path)) {
    write.csv(
      vif_history,
      output_path,
      row.names = FALSE
    )
  }

  list(
    predictors = current_preds,
    table = vif_history
  )
}


fit_glmm_workflow <- function(
    data_all,
    model_name,
    response_variable,
    predictors,
    vif_threshold,
    significance_threshold,
    base_dir
) {

  message("\n===================================================")
  message("STAGE 1: FITTING GAUSSIAN GLMM")
  message("Model: ", model_name)
  message("===================================================")

  # Set T1 as reference category.
  data_all$Duration <- factor(
    data_all$Duration,
    levels = c("T1", "T2")
  )

  predictors <- remove_zero_variance_predictors(
    data_all,
    predictors
  )

  needed_cols <- c(
    response_variable,
    predictors,
    "Duration",
    "Year"
  )

  data_cc <- data_all[
    complete.cases(
      data_all[, needed_cols, drop = FALSE]
    ),
    needed_cols,
    drop = FALSE
  ]

  if (nrow(data_cc) == 0) {
    stop(
      "No complete rows remain. Check missing values and column names."
    )
  }

  # ---------------------------------------------------------------------------
  # VIF SCREENING
  # ---------------------------------------------------------------------------

  vif_result <- iterative_vif_filter(
    df = data_cc,
    response_variable = response_variable,
    predictors = predictors,
    threshold = vif_threshold,
    output_path = file.path(
      base_dir,
      paste0(model_name, "_VIF_summary.csv")
    )
  )

  predictors <- vif_result$predictors

  if (length(predictors) == 0) {
    stop(
      "No LC transition predictors remain after VIF screening."
    )
  }

  # ---------------------------------------------------------------------------
  # GAUSSIAN GLMM
  # ---------------------------------------------------------------------------

  glmm_formula <- as.formula(
    paste(
      response_variable,
      "~",
      paste(predictors, collapse = " + "),
      "+ Duration + (1 | Year)"
    )
  )

  model <- glmmTMB(
    formula = glmm_formula,
    data = data_cc,
    family = gaussian(link = "identity")
  )

  model_summary <- summary(model)

  capture.output(
    model_summary,
    file = file.path(
      base_dir,
      paste0(model_name, "_GLMM_summary.txt")
    )
  )

  # ---------------------------------------------------------------------------
  # FIXED-EFFECT COEFFICIENT TABLE
  # ---------------------------------------------------------------------------

  coef_matrix <- model_summary$coefficients$cond

  coef_table <- data.frame(
    Term = rownames(coef_matrix),
    Estimate = coef_matrix[, "Estimate"],
    SE = coef_matrix[, "Std. Error"],
    z_value = coef_matrix[, "z value"],
    p_value = coef_matrix[, "Pr(>|z|)"],
    row.names = NULL
  )

  write.csv(
    coef_table,
    file.path(
      base_dir,
      paste0(model_name, "_all_fixed_effects.csv")
    ),
    row.names = FALSE
  )

  # Keep LC transition terms only.
  transition_results <- coef_table[
    coef_table$Term %in% predictors,
    ,
    drop = FALSE
  ]

  transition_results$Significant <- (
    transition_results$p_value < significance_threshold
  )

  write.csv(
    transition_results,
    file.path(
      base_dir,
      paste0(model_name, "_transition_results.csv")
    ),
    row.names = FALSE
  )

  significant_results <- transition_results[
    transition_results$Significant,
    ,
    drop = FALSE
  ]

  write.csv(
    significant_results,
    file.path(
      base_dir,
      paste0(model_name, "_significant_transitions.csv")
    ),
    row.names = FALSE
  )

  if (nrow(significant_results) == 0) {
    stop(
      "No LC transitions were significant at p < ",
      significance_threshold,
      "."
    )
  }

  # ---------------------------------------------------------------------------
  # MARGINAL AND CONDITIONAL R²
  # ---------------------------------------------------------------------------

  r2_result <- performance::r2_nakagawa(model)

  r2_table <- data.frame(
    Marginal_R2 = as.numeric(r2_result$R2_marginal),
    Conditional_R2 = as.numeric(r2_result$R2_conditional)
  )

  write.csv(
    r2_table,
    file.path(
      base_dir,
      paste0(model_name, "_R2.csv")
    ),
    row.names = FALSE
  )

  message(
    "Marginal R2 = ",
    round(r2_table$Marginal_R2, 3)
  )

  message(
    "Conditional R2 = ",
    round(r2_table$Conditional_R2, 3)
  )

  list(
    model = model,
    data = data_cc,
    predictors = predictors,
    transition_results = transition_results,
    significant_results = significant_results,
    r2 = r2_table
  )
}


prepare_presence_cells <- function(presence_file) {

  if (!file.exists(presence_file)) {
    stop("Presence-cell file not found: ", presence_file)
  }

  presence_cells <- read.csv(
    presence_file,
    stringsAsFactors = FALSE
  )

  if (!("Predictor" %in% names(presence_cells))) {
    stop(
      "Presence_cells.csv must contain a column named 'Predictor'."
    )
  }

  if (all(c("T1", "T2") %in% names(presence_cells))) {

    presence_cells$Presence_Cells <-
      as.numeric(presence_cells$T1) +
      as.numeric(presence_cells$T2)

  } else if ("Presence_Cells" %in% names(presence_cells)) {

    presence_cells$Presence_Cells <-
      as.numeric(presence_cells$Presence_Cells)

  } else {

    stop(
      "Presence_cells.csv must contain either T1 and T2 columns ",
      "or one Presence_Cells column."
    )
  }

  presence_cells[
    ,
    c("Predictor", "Presence_Cells"),
    drop = FALSE
  ]
}


calculate_cumulative_impact <- function(
    significant_results,
    presence_file,
    cell_area_km2,
    base_dir
) {

  message("\n===================================================")
  message("CALCULATING CUMULATIVE HEC IMPACT")
  message("===================================================")

  presence_cells <- prepare_presence_cells(
    presence_file
  )

  impact <- merge(
    significant_results,
    presence_cells,
    by.x = "Term",
    by.y = "Predictor",
    all.x = TRUE
  )

  if (any(is.na(impact$Presence_Cells))) {

    missing_names <- impact$Term[
      is.na(impact$Presence_Cells)
    ]

    stop(
      "Missing presence-cell counts for: ",
      paste(missing_names, collapse = ", ")
    )
  }

  # DeltaHEC = beta × n
  impact$DeltaHEC <- (
    impact$Estimate *
    impact$Presence_Cells
  )

  impact$Area_km2 <- (
    impact$Presence_Cells *
    cell_area_km2
  )

  impact <- impact[
    ,
    c(
      "Term",
      "Estimate",
      "SE",
      "z_value",
      "p_value",
      "Presence_Cells",
      "Area_km2",
      "DeltaHEC"
    )
  ]

  write.csv(
    impact,
    file.path(
      base_dir,
      "Cumulative_HEC_impact.csv"
    ),
    row.names = FALSE
  )

  impact
}


fit_translation_model <- function(
    impact_data,
    rmse_threshold_percent,
    significance_threshold,
    base_dir
) {

  message("\n===================================================")
  message("STAGE 2: TRANSLATION AND VERIFICATION")
  message("===================================================")

  if (nrow(impact_data) < 3) {
    stop(
      "At least three significant transitions are required ",
      "to fit the translation model."
    )
  }

  # ---------------------------------------------------------------------------
  # LOG-LINEAR TRANSLATION MODEL
  #
  # ln(n_obs + 1) = alpha + lambda(DeltaHEC)
  # ---------------------------------------------------------------------------

  translation_model <- lm(
    log(Presence_Cells + 1) ~ DeltaHEC,
    data = impact_data
  )

  capture.output(
    summary(translation_model),
    file = file.path(
      base_dir,
      "Translation_model_summary.txt"
    )
  )

  # Prediction on log scale.
  predicted_log <- predict(
    translation_model,
    newdata = impact_data
  )

  # Back-transform to presence-cell counts.
  impact_data$Predicted_Presence_Cells <-
    exp(predicted_log) - 1

  # Prevent tiny numerical negatives after back-transformation.
  impact_data$Predicted_Presence_Cells[
    impact_data$Predicted_Presence_Cells < 0
  ] <- 0

  # ---------------------------------------------------------------------------
  # VERIFICATION METRICS
  # ---------------------------------------------------------------------------

  residual_count <- (
    impact_data$Presence_Cells -
    impact_data$Predicted_Presence_Cells
  )

  rmse <- sqrt(
    mean(
      residual_count^2,
      na.rm = TRUE
    )
  )

  mean_observed <- mean(
    impact_data$Presence_Cells,
    na.rm = TRUE
  )

  nrmse_percent <- (
    rmse /
    mean_observed
  ) * 100

  model_coef <- summary(
    translation_model
  )$coefficients

  slope_p <- model_coef[
    "DeltaHEC",
    "Pr(>|t|)"
  ]

  model_r2 <- summary(
    translation_model
  )$r.squared

  accepted <- (
    slope_p < significance_threshold &&
    nrmse_percent <= rmse_threshold_percent
  )

  verification <- data.frame(
    Translation_Slope_p = slope_p,
    Translation_R2 = model_r2,
    RMSE_Count_Units = rmse,
    Mean_Observed_Count = mean_observed,
    NRMSE_Percent_of_Mean = nrmse_percent,
    Significance_Threshold = significance_threshold,
    RMSE_Threshold_Percent = rmse_threshold_percent,
    Representation = ifelse(
      accepted,
      "Accepted",
      "Needs refinement"
    )
  )

  write.csv(
    verification,
    file.path(
      base_dir,
      "Translation_verification_metrics.csv"
    ),
    row.names = FALSE
  )

  final_output <- impact_data[
    ,
    c(
      "Term",
      "Estimate",
      "SE",
      "p_value",
      "Presence_Cells",
      "Area_km2",
      "DeltaHEC",
      "Predicted_Presence_Cells"
    )
  ]

  names(final_output)[
    names(final_output) == "Presence_Cells"
  ] <- "Observed_Presence_Cells"

  write.csv(
    final_output,
    file.path(
      base_dir,
      "Observed_vs_Predicted_presence_cells.csv"
    ),
    row.names = FALSE
  )

  message(
    "Translation slope p = ",
    signif(slope_p, 4)
  )

  message(
    "RMSE = ",
    round(rmse, 2),
    " presence cells"
  )

  message(
    "Normalized RMSE = ",
    round(nrmse_percent, 2),
    "%"
  )

  message(
    "Representation: ",
    verification$Representation
  )

  list(
    model = translation_model,
    predictions = final_output,
    verification = verification
  )
}

# =============================================================================
# 2) READ T1 AND T2 DATA
# =============================================================================

data_t1 <- read_timeline(
  files$T1,
  "T1",
  predictors,
  response_variable
)

data_t2 <- read_timeline(
  files$T2,
  "T2",
  predictors,
  response_variable
)

# Combine the two temporally aligned datasets.
data_combined <- rbind(
  data_t1,
  data_t2
)

# =============================================================================
# 3) STAGE 1 — GLMM
# =============================================================================

# INTERPRETATION OF GLMM COEFFICIENTS
# Positive beta:
#   The LC transition is associated with higher HEC intensity.
#
# Negative beta:
#   The LC transition is associated with lower HEC intensity.
#
# Statistical significance is evaluated at p < 0.05.

glmm_result <- fit_glmm_workflow(
  data_all = data_combined,
  model_name = "T1_T2_HEC_GLMM",
  response_variable = response_variable,
  predictors = predictors,
  vif_threshold = vif_threshold,
  significance_threshold = significance_threshold,
  base_dir = base_dir
)

# =============================================================================
# 4) CUMULATIVE PREDICTED HEC IMPACT
# =============================================================================

impact_result <- calculate_cumulative_impact(
  significant_results = glmm_result$significant_results,
  presence_file = presence_file,
  cell_area_km2 = cell_area_km2,
  base_dir = base_dir
)

# =============================================================================
# 5) STAGE 2 — TRANSLATION AND VERIFICATION
# =============================================================================

# NOTE:
# The translation model provides an internal verification of agreement between
# cumulative HEC impact and observed presence-cell counts. It should not be
# interpreted as independent or external predictive validation.

translation_result <- fit_translation_model(
  impact_data = impact_result,
  rmse_threshold_percent = rmse_threshold_percent,
  significance_threshold = significance_threshold,
  base_dir = base_dir
)

# =============================================================================
# OUTPUT FILES
# =============================================================================
#
# T1_T2_HEC_GLMM_VIF_summary.csv
#   Full VIF screening history.
#
# T1_T2_HEC_GLMM_GLMM_summary.txt
#   Complete Gaussian GLMM summary.
#
# T1_T2_HEC_GLMM_all_fixed_effects.csv
#   All fixed-effect terms, including Duration.
#
# T1_T2_HEC_GLMM_transition_results.csv
#   LC transition coefficients, SE, z statistics, p-values and significance.
#
# T1_T2_HEC_GLMM_significant_transitions.csv
#   Only LC transitions with p < 0.05.
#
# T1_T2_HEC_GLMM_R2.csv
#   Marginal and conditional R².
#
# Cumulative_HEC_impact.csv
#   Significant transition coefficients, observed presence-cell counts,
#   spatial extent and DeltaHEC = beta × n.
#
# Translation_model_summary.txt
#   Log-linear translation model:
#       ln(n_obs + 1) = alpha + lambda(DeltaHEC)
#
# Translation_verification_metrics.csv
#   Translation p-value, R², RMSE, normalized RMSE (%) and acceptance status.
#
# Observed_vs_Predicted_presence_cells.csv
#   Final observed and predicted HEC presence-cell counts for Fig. 4 / reporting.
#
# =============================================================================

cat("\n===================================================\n")
cat("WORKFLOW COMPLETE\n")
cat("===================================================\n")
cat("Outputs saved to:\n", base_dir, "\n")
