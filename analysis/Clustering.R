rm(list = ls())

library(dplyr)
library(survival)

# PATHS
project_dir <- "/Users/WalnutDragon/Downloads/student-template"
data_path <- file.path(project_dir, "data/data_model.rds")
results_dir <- file.path(project_dir, "results")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# LOAD DATA
if (!file.exists(data_path)) {
  stop("data_model.rds not found at: ", data_path)
}

data_model <- readRDS(data_path)

cat("Loaded:", dim(data_model), "\n")

# FEATURE PREP

data_model$STAGE4 <- ifelse(data_model$STAGE == "Stage 4", 1, 0)

cluster_features <- data_model %>%
  select(AGE, STAGE4, TP53, KRAS, APC, EGFR, BRAF, PTEN, RBM10, PIK3CA)

cluster_features_scaled <- scale(cluster_features)

# PCA
pca <- prcomp(cluster_features_scaled)

data_model$PC1 <- pca$x[,1]
data_model$PC2 <- pca$x[,2]

# KMEANS

set.seed(123)
km <- kmeans(cluster_features_scaled, centers = 3, nstart = 25)

data_model$cluster_raw <- as.factor(km$cluster)

print(table(data_model$cluster_raw))

# RISK LABELING

cluster_medians <- data_model %>%
  group_by(cluster_raw) %>%
  summarise(median_os = median(OS_TIME), .groups = "drop") %>%
  arrange(desc(median_os))

labels <- c("Low-Risk", "Intermediate-Risk", "High-Risk")

map <- setNames(labels, cluster_medians$cluster_raw)

data_model$class_bio <- map[as.character(data_model$cluster_raw)]

print(table(data_model$class_bio))

# SURVIVAL

fit <- survfit(Surv(OS_TIME, OS_EVENT) ~ class_bio, data = data_model)

print(survdiff(Surv(OS_TIME, OS_EVENT) ~ class_bio, data = data_model))

cox <- coxph(Surv(OS_TIME, OS_EVENT) ~ class_bio, data = data_model)
print(summary(cox))

# SAVE

saveRDS(data_model, file.path(results_dir, "clustered_data_model.rds"))

cat("Saved clustered data\n")

