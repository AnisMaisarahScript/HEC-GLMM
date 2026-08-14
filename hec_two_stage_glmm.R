# =========================================================
# TWO-STAGE HEC GLMM WORKFLOW
#   A) Concurrent model: T1 + T2 only
#   B) Delayed model: T3 only
#   Stage 2 calibration: predicted KDE-like impact -> HEC cases
# =========================================================

library(lme4)
library(car)
library(MASS)

# ---------------------------
# 0) USER SETTINGS
# ---------------------------

base_dir <- "C:/Users/HP-15FC/Desktop/Software/new R/Revised GLMM (Include ROC and PI)"

files <- list(
  T1 = file.path(base_dir, "T1.csv"),
  T2 = file.path(base_dir, "T2.csv"),
  T3 = file.path(base_dir, "T3.csv")
)

presence_file <- file.path(base_dir, "Presence_cells.csv")
calibration_file <- file.path(base_dir, "Calibration_HEC_cases.csv")

response_variable <- "HEC_kernel"
vif_threshold <- 5
cell_area_km2 <- 0.0625

predictors <- c(
  "WTOF","WTOS","WTOI","WTOP",
  "STOF","STOW","STOI","STOP",
  "ITOF","ITOW","ITOS","ITOP",
  "PTOF","PTOW","PTOS","PTOI",
  "FTOP","FTOW","FTOS","FTOI"
)

predictors <- unique(predictors)

# ---------------------------
# 1) HELPERS
# ---------------------------

rename_response_if_needed <- function(df, response_variable) {
  possible_response_names <- c(
    response_variable,
    "old_name",
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
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- rename_response_if_needed(df, response_variable)
  df$Duration <- label

  required_non_predictors <- c(response_variable, "Year", "Duration")
  missing_required <- setdiff(required_non_predictors, names(df))
  if (length(missing_required) > 0) {
    stop(
      "Missing required columns in ", basename(path), ": ",
      paste(missing_required, collapse = ", ")
    )
  }

  missing_predictors <- setdiff(predictors, names(df))
  if (length(missing_predictors) > 0) {
    message(
      "Adding missing predictor columns as 0 in ", label, ": ",
      paste(missing_predictors, collapse = ", ")
    )
    for (p in missing_predictors) {
      df[[p]] <- 0
    }
  }

  needed_cols <- c(response_variable, predictors, "Duration", "Year")
  df <- df[, needed_cols, drop = FALSE]

  for (p in predictors) {
    df[[p]] <- as.numeric(df[[p]])
  }

  df[[response_variable]] <- as.numeric(df[[response_variable]])
  df$Year <- as.factor(df$Year)
  df$Duration <- as.factor(df$Duration)

  df
}

remove_zero_variance_predictors <- function(df, predictors) {
  var_x <- sapply(df[, predictors, drop = FALSE], var, na.rm = TRUE)
  zero_var_preds <- names(var_x[var_x == 0 | is.na(var_x)])

  if (length(zero_var_preds) > 0) {
    message("Removed zero-variance predictors: ", paste(zero_var_preds, collapse = ", "))
  } else {
    message("No zero-variance predictors found.")
  }

  setdiff(predictors, zero_var_preds)
}

iterative_vif_filter <- function(df, response_variable, predictors, include_duration,
                                 threshold = 5, output_path = NULL) {
  current_preds <- predictors
  vif_history <- data.frame()

  if (length(current_preds) < 2) {
    return(list(predictors = current_preds, table = vif_history))
  }

  repeat {
    fixed_terms <- current_preds
    if (include_duration) {
      fixed_terms <- c(fixed_terms, "Duration")
    }

    vif_formula <- as.formula(
      paste(response_variable, "~", paste(fixed_terms, collapse = " + "))
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

    predictor_vifs <- vif_table[vif_table$Predictor %in% current_preds, , drop = FALSE]
    predictor_vifs <- predictor_vifs[order(-predictor_vifs$VIF), , drop = FALSE]

    vif_table$Iteration <- length(unique(vif_history$Iteration)) + 1
    vif_table$Decision <- "Kept"

    if (nrow(predictor_vifs) == 0 || max(predictor_vifs$VIF, na.rm = TRUE) <= threshold) {
      vif_history <- rbind(vif_history, vif_table)
      break
    }

    drop_pred <- predictor_vifs$Predictor[1]
    vif_table$Decision[vif_table$Predictor == drop_pred] <- "Dropped"
    vif_history <- rbind(vif_history, vif_table)

    message("Removed predictor due to high VIF: ", drop_pred)
    current_preds <- setdiff(current_preds, drop_pred)

    if (length(current_preds) < 2) {
      break
    }
  }

  if (!is.null(output_path)) {
    write.csv(vif_history, output_path, row.names = FALSE)
  }

  list(predictors = current_preds, table = vif_history)
}

fit_lmm_workflow <- function(data_all, model_name, response_variable, predictors,
                             include_duration, vif_threshold, base_dir) {
  message("\n==============================")
  message("Fitting ", model_name)
  message("==============================")

  if (include_duration) {
    data_all$Duration <- factor(data_all$Duration, levels = c("T1", "T2"))
  }

  predictors <- remove_zero_variance_predictors(data_all, predictors)

  needed_cols <- c(response_variable, predictors, "Year")
  if (include_duration) {
    needed_cols <- c(needed_cols, "Duration")
  }

  data_cc <- data_all[
    complete.cases(data_all[, needed_cols, drop = FALSE]),
    needed_cols,
    drop = FALSE
  ]

  if (nrow(data_cc) == 0) {
    stop("No complete rows remain for ", model_name, ". Check missing values and column names.")
  }

  vif_result <- iterative_vif_filter(
    df = data_cc,
    response_variable = response_variable,
    predictors = predictors,
    include_duration = include_duration,
    threshold = vif_threshold,
    output_path = file.path(base_dir, paste0(model_name, "_VIF_summary.csv"))
  )

  predictors <- vif_result$predictors

  if (length(predictors) == 0) {
    stop("No predictors remain after zero-variance and VIF filtering for ", model_name)
  }

  fixed_terms <- predictors
  if (include_duration) {
    fixed_terms <- c(fixed_terms, "Duration")
  }

  lmm_formula <- as.formula(
    paste(response_variable, "~", paste(fixed_terms, collapse = " + "), "+ (1 | Year)")
  )

  model <- lmer(lmm_formula, data = data_cc, REML = FALSE)
  model_summary <- summary(model)

  capture.output(
    model_summary,
    file = file.path(base_dir, paste0(model_name, "_lmm_summary.txt"))
  )

  if (isSingular(model, tol = 1e-4)) {
    warning(
      model_name,
      " has a singular random-effect fit. Consider Year as a fixed effect, ",
      "or compare with a non-random Year model."
    )
  }

  fixed_effects <- fixef(model)
  coef_out <- data.frame(
    Predictor = predictors,
    Coefficient = as.numeric(fixed_effects[predictors]),
    row.names = NULL
  )

  write.csv(
    coef_out,
    file.path(base_dir, paste0(model_name, "_coefficients.csv")),
    row.names = FALSE
  )

  write.csv(
    data.frame(Final_Predictors = predictors),
    file.path(base_dir, paste0(model_name, "_final_predictors.csv")),
    row.names = FALSE
  )

  list(
    model = model,
    data = data_cc,
    predictors = predictors,
    coefficients = coef_out
  )
}

prepare_presence_cells <- function(presence_file, model_type) {
  presence_cells <- read.csv(presence_file, stringsAsFactors = FALSE)

  if (!("Predictor" %in% names(presence_cells))) {
    stop("Presence file must contain a Predictor column.")
  }

  if (model_type == "Concurrent") {
    if (all(c("T1", "T2") %in% names(presence_cells))) {
      presence_cells$Presence_Cells <- presence_cells$T1 + presence_cells$T2
    } else if (!("Presence_Cells" %in% names(presence_cells))) {
      stop("Concurrent presence file must contain T1 and T2 columns, or Presence_Cells.")
    }
  }

  if (model_type == "Delayed") {
    if ("T3" %in% names(presence_cells)) {
      presence_cells$Presence_Cells <- presence_cells$T3
    } else if (!("Presence_Cells" %in% names(presence_cells))) {
      stop("Delayed presence file must contain T3 or Presence_Cells.")
    }
  }

  presence_cells[, c("Predictor", "Presence_Cells"), drop = FALSE]
}

calculate_predicted_impact <- function(coef_out, presence_file, model_type,
                                       cell_area_km2, base_dir) {
  presence_cells <- prepare_presence_cells(presence_file, model_type)

  impact_out <- merge(coef_out, presence_cells, by = "Predictor", all.x = TRUE)

  if (any(is.na(impact_out$Presence_Cells))) {
    missing_presence <- impact_out$Predictor[is.na(impact_out$Presence_Cells)]
    stop(
      "Missing presence-cell counts for: ",
      paste(missing_presence, collapse = ", ")
    )
  }

  impact_out$Predicted_Impact <- impact_out$Coefficient * impact_out$Presence_Cells
  impact_out$Area_km2 <- impact_out$Presence_Cells * cell_area_km2

  write.csv(
    impact_out,
    file.path(base_dir, paste0(model_type, "_Predicted_impact.csv")),
    row.names = FALSE
  )

  impact_out
}

fit_calibration_model <- function(calibration_file, base_dir) {
  calibration_data <- read.csv(calibration_file, stringsAsFactors = FALSE)

  required_cols <- c("Observed_KDE", "HEC_Cases")
  missing_cols <- setdiff(required_cols, names(calibration_data))
  if (length(missing_cols) > 0) {
    stop("Calibration file is missing: ", paste(missing_cols, collapse = ", "))
  }

  case_model_poisson <- glm(
    HEC_Cases ~ Observed_KDE,
    family = poisson(link = "log"),
    data = calibration_data
  )

  dispersion_ratio <- sum(residuals(case_model_poisson, type = "pearson")^2) /
    case_model_poisson$df.residual

  if (dispersion_ratio > 1.5) {
    case_model <- glm.nb(HEC_Cases ~ Observed_KDE, data = calibration_data)
    model_type <- "Negative Binomial"
  } else {
    case_model <- case_model_poisson
    model_type <- "Poisson"
  }

  capture.output(
    paste("Calibration model used:", model_type),
    paste("Poisson dispersion ratio:", round(dispersion_ratio, 3)),
    summary(case_model),
    file = file.path(base_dir, "Calibration_model_summary.txt")
  )

  list(
    model = case_model,
    model_type = model_type,
    dispersion_ratio = dispersion_ratio
  )
}

estimate_cases <- function(impact_out, calibration_result, model_type, base_dir) {
  case_model <- calibration_result$model

  model_frame <- model.frame(case_model)
  observed_kde_range <- range(model_frame$Observed_KDE, na.rm = TRUE)
  predicted_impact_range <- range(impact_out$Predicted_Impact, na.rm = TRUE)

  if (
    predicted_impact_range[1] < observed_kde_range[1] ||
      predicted_impact_range[2] > observed_kde_range[2]
  ) {
    warning(
      model_type,
      " predicted impacts extend outside the calibration Observed_KDE range. ",
      "Treat case estimates as extrapolated."
    )
  }

  prediction_data <- data.frame(
    Observed_KDE = impact_out$Predicted_Impact
  )

  baseline_data <- data.frame(Observed_KDE = 0)
  baseline_cases <- as.numeric(predict(case_model, newdata = baseline_data, type = "response"))

  impact_out$Estimated_HEC_Cases <- as.numeric(
    predict(case_model, newdata = prediction_data, type = "response")
  )

  impact_out$Estimated_Additional_HEC_Cases <- impact_out$Estimated_HEC_Cases - baseline_cases
  impact_out$Estimated_HEC_Cases_Rounded <- round(impact_out$Estimated_HEC_Cases, 0)
  impact_out$Estimated_Additional_HEC_Cases_Rounded <- round(
    impact_out$Estimated_Additional_HEC_Cases,
    0
  )

  final_cols <- c(
    "Predictor",
    "Coefficient",
    "Presence_Cells",
    "Area_km2",
    "Predicted_Impact",
    "Estimated_HEC_Cases",
    "Estimated_HEC_Cases_Rounded",
    "Estimated_Additional_HEC_Cases",
    "Estimated_Additional_HEC_Cases_Rounded"
  )

  impact_out <- impact_out[, final_cols, drop = FALSE]

  write.csv(
    impact_out,
    file.path(base_dir, paste0(model_type, "_Final_predicted_impact_HEC_cases.csv")),
    row.names = FALSE
  )

  impact_out
}

# ---------------------------
# 2) READ DATA
# ---------------------------

data_t1 <- read_timeline(files$T1, "T1", predictors, response_variable)
data_t2 <- read_timeline(files$T2, "T2", predictors, response_variable)
data_t3 <- read_timeline(files$T3, "T3", predictors, response_variable)

# ---------------------------
# 3) STAGE 1: GLMM MODELS
# ---------------------------

data_concurrent <- rbind(data_t1, data_t2)
data_delayed <- data_t3

concurrent_result <- fit_lmm_workflow(
  data_all = data_concurrent,
  model_name = "Concurrent_T1_T2",
  response_variable = response_variable,
  predictors = predictors,
  include_duration = TRUE,
  vif_threshold = vif_threshold,
  base_dir = base_dir
)

delayed_result <- fit_lmm_workflow(
  data_all = data_delayed,
  model_name = "Delayed_T3",
  response_variable = response_variable,
  predictors = predictors,
  include_duration = FALSE,
  vif_threshold = vif_threshold,
  base_dir = base_dir
)

# ---------------------------
# 4) STAGE 1 OUTPUT: PREDICTED IMPACT
# ---------------------------

concurrent_impact <- calculate_predicted_impact(
  coef_out = concurrent_result$coefficients,
  presence_file = presence_file,
  model_type = "Concurrent",
  cell_area_km2 = cell_area_km2,
  base_dir = base_dir
)

delayed_impact <- calculate_predicted_impact(
  coef_out = delayed_result$coefficients,
  presence_file = presence_file,
  model_type = "Delayed",
  cell_area_km2 = cell_area_km2,
  base_dir = base_dir
)

# ---------------------------
# 5) STAGE 2: TRANSLATION
# ---------------------------

calibration_result <- fit_calibration_model(calibration_file, base_dir)

concurrent_final <- estimate_cases(
  impact_out = concurrent_impact,
  calibration_result = calibration_result,
  model_type = "Concurrent",
  base_dir = base_dir
)

delayed_final <- estimate_cases(
  impact_out = delayed_impact,
  calibration_result = calibration_result,
  model_type = "Delayed",
  base_dir = base_dir
)

cat("\nWorkflow complete.\n")
cat("Outputs saved to:", base_dir, "\n")
