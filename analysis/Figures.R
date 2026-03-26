
# PACKAGES

library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)


# OUTPUT FOLDER

outdir <- "/Users/WalnutDragon/Downloads/student-template/docs/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# BASIC CHECKS

required_cols <- c("OS_TIME", "OS_EVENT", "AGE", "STAGE", "class_bio",
                   "TP53", "KRAS", "APC", "EGFR", "BRAF", "PTEN", "RBM10")

missing_cols <- setdiff(required_cols, names(data_model))
if (length(missing_cols) > 0) {
  stop("Missing required columns in data_model: ", paste(missing_cols, collapse = ", "))
}

data_model$class_bio <- as.factor(data_model$class_bio)


if (length(levels(data_model$class_bio)) == 3) {
  levels(data_model$class_bio) <- c(
    "Cluster 1 (Low-Risk)",
    "Cluster 2 (Intermediate-Risk)",
    "Cluster 3 (High-Risk)"
  )
}

# PREPARE FEATURES FOR PCA / CLUSTER VIS

cluster_features <- data_model %>%
  transmute(
    AGE = AGE,
    STAGE4 = ifelse(STAGE == "Stage 4", 1, 0),
    TP53 = TP53,
    KRAS = KRAS,
    APC = APC,
    EGFR = EGFR,
    BRAF = BRAF,
    PTEN = PTEN,
    RBM10 = RBM10,
    PIK3CA = if ("PIK3CA" %in% names(data_model)) PIK3CA else 0
  )

cluster_features_scaled <- scale(cluster_features)

pca_obj <- prcomp(cluster_features_scaled, center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  Cluster = data_model$class_bio
)


# FIGURE 1: OVERALL SURVIVAL KM

fit_overall <- survfit(Surv(OS_TIME, OS_EVENT) ~ 1, data = data_model)

km_overall <- ggsurvplot(
  fit_overall,
  data = data_model,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "Overall Survival in the MSK-CHORD Baseline Cohort",
  xlab = "Time (months)",
  ylab = "Survival probability"
)

ggsave(
  filename = file.path(outdir, "figure1.png"),
  plot = km_overall$plot,
  width = 7,
  height = 5,
  dpi = 144
)


# FIGURE 2: PCA STRUCTURE

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 1.8) +
  labs(
    title = "PCA Projection of Clinicogenomic Features",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "figure2.png"),
  plot = p2,
  width = 7,
  height = 5,
  dpi = 144
)

# FIGURE 3: ELBOW PLOT
set.seed(123)

wss <- sapply(2:8, function(k) {
  kmeans(cluster_features_scaled, centers = k, nstart = 25)$tot.withinss
})

elbow_df <- data.frame(
  k = 2:8,
  wss = wss
)

p3 <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Elbow Plot for Cluster Selection",
    x = "Number of clusters",
    y = "Within-cluster sum of squares"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "figure3.png"),
  plot = p3,
  width = 7,
  height = 5,
  dpi = 144
)

# FIGURE 4: KM BY CLUSTER

fit_cluster <- survfit(Surv(OS_TIME, OS_EVENT) ~ class_bio, data = data_model)

km_cluster <- ggsurvplot(
  fit_cluster,
  data = data_model,
  risk.table = TRUE,
  pval = TRUE,
  title = "Overall Survival Across Clinicogenomic Clusters",
  xlab = "Time (months)",
  ylab = "Survival probability",
  legend.title = "Cluster"
)

ggsave(
  filename = file.path(outdir, "figure4.png"),
  plot = km_cluster$plot,
  width = 7,
  height = 5,
  dpi = 144
)


# FIGURE 5: AGE ACROSS CLUSTERS

p5 <- ggplot(data_model, aes(x = class_bio, y = AGE, fill = class_bio)) +
  geom_boxplot(alpha = 0.85) +
  labs(
    title = "Age Distribution Across Clusters",
    x = "Cluster",
    y = "Age"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(outdir, "figure5.png"),
  plot = p5,
  width = 7,
  height = 5,
  dpi = 144
)

# FIGURE 6: TUMOR STAGE ACROSS CLUSTERS

stage_df <- data_model %>%
  mutate(StageGroup = ifelse(STAGE == "Stage 4", "Stage 4", "Stage 1-3")) %>%
  count(class_bio, StageGroup) %>%
  group_by(class_bio) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p6 <- ggplot(stage_df, aes(x = class_bio, y = prop, fill = StageGroup)) +
  geom_col(position = "fill") +
  labs(
    title = "Tumor Stage Distribution Across Clusters",
    x = "Cluster",
    y = "Proportion",
    fill = "Stage"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "figure6.png"),
  plot = p6,
  width = 7,
  height = 5,
  dpi = 144
)

# FIGURE 7: MUTATION PREVALENCE ACROSS CLUSTERS

genes_to_plot <- c("TP53", "KRAS", "APC", "EGFR", "BRAF", "PTEN", "RBM10")
genes_to_plot <- genes_to_plot[genes_to_plot %in% names(data_model)]

mut_cluster_df <- data_model %>%
  group_by(class_bio) %>%
  summarise(across(all_of(genes_to_plot), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "Gene",
    values_to = "MutationFrequency"
  )

p7 <- ggplot(mut_cluster_df, aes(x = Gene, y = MutationFrequency, fill = class_bio)) +
  geom_col(position = "dodge") +
  labs(
    title = "Mutation Prevalence Across Clusters",
    x = "Gene",
    y = "Mutation frequency",
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "figure7.png"),
  plot = p7,
  width = 8,
  height = 5,
  dpi = 144
)

# FIGURE 8: OVERALL MUTATION FREQUENCY

mut_overall_df <- data.frame(
  Gene = genes_to_plot,
  MutationFrequency = sapply(data_model[genes_to_plot], mean, na.rm = TRUE)
)

p8 <- ggplot(mut_overall_df, aes(x = reorder(Gene, MutationFrequency), y = MutationFrequency)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Overall Mutation Frequency of Key Driver Genes",
    x = "Gene",
    y = "Mutation frequency"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(outdir, "figure8.png"),
  plot = p8,
  width = 7,
  height = 5,
  dpi = 144
)


# FIGURE 9: CLUSTER SIZES
cluster_size_df <- data_model %>%
  count(class_bio)

p9 <- ggplot(cluster_size_df, aes(x = class_bio, y = n, fill = class_bio)) +
  geom_col() +
  labs(
    title = "Cluster Sizes",
    x = "Cluster",
    y = "Number of patients"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(outdir, "figure9.png"),
  plot = p9,
  width = 7,
  height = 5,
  dpi = 144
)

# ggsave(file.path(outdir, "figure1.svg"), plot = km_overall$plot, width = 7, height = 5)
# ggsave(file.path(outdir, "figure2.svg"), plot = p2, width = 7, height = 5)
# ggsave(file.path(outdir, "figure3.svg"), plot = p3, width = 7, height = 5)
# ggsave(file.path(outdir, "figure4.svg"), plot = km_cluster$plot, width = 7, height = 5)
# ggsave(file.path(outdir, "figure5.svg"), plot = p5, width = 7, height = 5)
# ggsave(file.path(outdir, "figure6.svg"), plot = p6, width = 7, height = 5)
# ggsave(file.path(outdir, "figure7.svg"), plot = p7, width = 8, height = 5)
# ggsave(file.path(outdir, "figure8.svg"), plot = p8, width = 7, height = 5)
# ggsave(file.path(outdir, "figure9.svg"), plot = p9, width = 7, height = 5)


# CONFIRM FILES SAVED

print(list.files(outdir))

