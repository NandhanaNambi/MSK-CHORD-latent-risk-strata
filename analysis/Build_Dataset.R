
# 01_build_dataset.R
# Purpose:
#   Validate and save the already-created patient-level
#   dataset (data_model) for the project pipeline.
#
# Notes:
#   This script does NOT rebuild from raw mutation files.
#   It assumes you already created a correct `data_model`
#   object in R, OR you have data_model.csv available.


rm(list = ls())

library(dplyr)
library(readr)

# Project paths

project_dir <- "/Users/WalnutDragon/Downloads/student-template"
data_dir <- file.path(project_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(data_dir, "data_model.rds")
csv_path <- file.path(data_dir, "data_model.csv")


# Load dataset
# Priority:
#   1. existing object in memory
#   2. existing CSV in project folder

if (exists("data_model")) {
  message("Using existing `data_model` object from memory.")
} else if (file.exists(csv_path)) {
  message("Loading data_model from CSV: ", csv_path)
  data_model <- read_csv(csv_path, show_col_types = FALSE)
} else {
  stop(
    "No `data_model` object in memory and no data_model.csv found at:\n",
    csv_path,
    "\n\nCreate/load your working dataset first, then rerun this script."
  )
}

# Required columns

required_cols <- c(
  "PATIENT_ID", "OS_TIME", "OS_EVENT", "AGE", "SEX", "STAGE",
  "TP53", "KRAS", "APC", "EGFR", "BRAF", "PTEN", "RBM10", "PIK3CA"
)

missing_cols <- setdiff(required_cols, names(data_model))
if (length(missing_cols) > 0) {
  stop("Missing required columns in data_model: ",
       paste(missing_cols, collapse = ", "))
}

# Standardize types

data_model <- data_model %>%
  mutate(
    PATIENT_ID = as.character(PATIENT_ID),
    SEX = as.character(SEX),
    STAGE = as.character(STAGE),
    OS_TIME = as.numeric(OS_TIME),
    OS_EVENT = as.integer(OS_EVENT),
    AGE = as.numeric(AGE)
  )

gene_cols <- c("TP53", "KRAS", "APC", "EGFR", "BRAF", "PTEN", "RBM10", "PIK3CA")

for (g in gene_cols) {
  data_model[[g]] <- as.numeric(data_model[[g]])
  data_model[[g]][is.na(data_model[[g]])] <- 0
}

# Derived variable

if (!("STAGE4" %in% names(data_model))) {
  data_model$STAGE4 <- ifelse(data_model$STAGE == "Stage 4", 1, 0)
}

# Remove exact duplicate patients if any
# Keep first occurrence

data_model <- data_model %>%
  distinct(PATIENT_ID, .keep_all = TRUE)

# Validation summary

cat("\nFinal data_model dimensions:", dim(data_model), "\n")
cat("Unique patients:", dplyr::n_distinct(data_model$PATIENT_ID), "\n")

cat("\nSummary of OS_TIME:\n")
print(summary(data_model$OS_TIME))

cat("\nOS_EVENT counts:\n")
print(table(data_model$OS_EVENT, useNA = "ifany"))

cat("\nStage counts:\n")
print(table(data_model$STAGE, useNA = "ifany"))

cat("\nMutation frequencies:\n")
print(round(colMeans(data_model[, gene_cols], na.rm = TRUE), 3))

# Save outputs

saveRDS(data_model, rds_path)
write.csv(data_model, csv_path, row.names = FALSE)

cat("\nSaved files:\n")
cat(rds_path, "\n")
cat(csv_path, "\n")

cat("\nDone.\n")

