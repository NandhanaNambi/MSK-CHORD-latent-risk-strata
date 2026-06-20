
# PHASE 2: Survival Modeling of Genomic Risk

# 0. Packages

packages <- c(
  "tidyverse",
  "survival",
  "broom",
  "forcats",
  "scales",
  "glue"
)

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(packages, library, character.only = TRUE))

set.seed(123)

# 1. Output folders

fig_dir <- "docs/figures"
results_dir <- "results"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# 2. Load data

data_candidates <- c(
  "data_model.rds",
  "data/data_model.rds",
  "C:/Users/nandh/Downloads/data_model.rds"
)

data_path <- data_candidates[file.exists(data_candidates)][1]

if (is.na(data_path)) {
  stop("data_model.rds not found. Check whether the file exists in Downloads or project folder.")
}

data_model <- readRDS(data_path)

cat("Loaded data from:", data_path, "\n")
cat("Rows:", nrow(data_model), "\n")
cat("Columns:", ncol(data_model), "\n")


# 3. Required variables

required_vars <- c("OS_TIME", "OS_EVENT", "AGE", "CANCER_TYPE")
missing_required <- setdiff(required_vars, names(data_model))

if (length(missing_required) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_required, collapse = ", ")
  )
}

stage_candidates <- c(
  "STAGE",
  "stage",
  "AJCC_STAGE",
  "AJCC_PATHOLOGIC_STAGE",
  "TUMOR_STAGE",
  "CLINICAL_STAGE"
)

stage_col <- intersect(stage_candidates, names(data_model))[1]

if (is.na(stage_col)) {
  stop(
    "No stage column found. Expected one of: ",
    paste(stage_candidates, collapse = ", ")
  )
}

cat("Using stage column:", stage_col, "\n")

# 4. Helper functions

clean_event <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  
  if (is.numeric(x)) {
    return(as.integer(x > 0))
  }
  
  z <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  
  dplyr::case_when(
    z %in% c("1", "yes", "y", "true", "dead", "deceased", "event") ~ 1L,
    z %in% c("0", "no", "n", "false", "alive", "censored") ~ 0L,
    TRUE ~ NA_integer_
  )
}

clean_mutation <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  
  if (is.numeric(x)) {
    return(as.integer(dplyr::coalesce(x, 0) > 0))
  }
  
  z <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  
  dplyr::case_when(
    z %in% c("1", "yes", "y", "true", "mut", "mutated", "mutation", "altered") ~ 1L,
    z %in% c("0", "no", "n", "false", "wt", "wildtype", "wild-type", "none") ~ 0L,
    is.na(z) | z == "" ~ 0L,
    TRUE ~ 0L
  )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

# 5. Define genes

target_genes <- c(
  "TP53",
  "KRAS",
  "APC",
  "EGFR",
  "BRAF",
  "PTEN",
  "RBM10",
  "PIK3CA"
)

gene_vars <- intersect(target_genes, names(data_model))
missing_genes <- setdiff(target_genes, gene_vars)

if (length(gene_vars) == 0) {
  stop("None of the target mutation genes were found in data_model.")
}

if (length(missing_genes) > 0) {
  warning(
    "These genes were not found and will be skipped: ",
    paste(missing_genes, collapse = ", ")
  )
}

cat("Genes included in Phase 2:", paste(gene_vars, collapse = ", "), "\n")

# 6. Clean analysis dataset

df2 <- data_model %>%
  mutate(
    OS_TIME = as.numeric(OS_TIME),
    OS_EVENT = clean_event(OS_EVENT),
    AGE = as.numeric(AGE),
    STAGE = stringr::str_trim(as.character(.data[[stage_col]])),
    STAGE = na_if(STAGE, ""),
    STAGE = na_if(STAGE, "NA"),
    STAGE = as.factor(STAGE),
    CANCER_TYPE = as.factor(CANCER_TYPE),
    across(all_of(gene_vars), clean_mutation)
  ) %>%
  filter(
    !is.na(OS_TIME),
    !is.na(OS_EVENT),
    !is.na(AGE),
    !is.na(STAGE),
    !is.na(CANCER_TYPE),
    OS_TIME > 0
  )

cat("Final Phase 2 dataset rows:", nrow(df2), "\n")
cat("Number of deaths/events:", sum(df2$OS_EVENT == 1, na.rm = TRUE), "\n")


# 7. Table 1: analytic cohort summary

phase2_summary <- df2 %>%
  summarise(
    N = n(),
    deaths = sum(OS_EVENT == 1, na.rm = TRUE),
    event_rate_percent = round(mean(OS_EVENT == 1, na.rm = TRUE) * 100, 1),
    age_mean = round(mean(AGE, na.rm = TRUE), 1),
    age_median = round(median(AGE, na.rm = TRUE), 1),
    cancer_types = n_distinct(CANCER_TYPE),
    stages = n_distinct(STAGE)
  )

print(phase2_summary)

readr::write_csv(
  phase2_summary,
  file.path(results_dir, "phase2_cohort_summary.csv")
)

# RESULT 1: Mutation prevalence

gene_prevalence <- purrr::map_dfr(gene_vars, function(g) {
  tibble(
    gene = g,
    mutated_n = sum(df2[[g]] == 1, na.rm = TRUE),
    total_n = sum(!is.na(df2[[g]])),
    prevalence = mutated_n / total_n
  )
})

print(gene_prevalence)

readr::write_csv(
  gene_prevalence,
  file.path(results_dir, "gene_prevalence.csv")
)

p_gene_prev <- gene_prevalence %>%
  mutate(
    gene = forcats::fct_reorder(gene, prevalence)
  ) %>%
  ggplot(aes(x = prevalence, y = gene)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Mutation Prevalence Across Selected Driver Genes",
    x = "Prevalence",
    y = NULL
  ) +
  theme_bw(base_size = 12)

print(p_gene_prev)

ggsave(
  filename = file.path(fig_dir, "phase2_gene_prevalence.png"),
  plot = p_gene_prev,
  width = 7,
  height = 5,
  dpi = 300
)

# RESULT 2: Gene-level adjusted Cox models

fit_binary_cox <- function(data, feature, model_name, feature_label = feature,
                           min_mutated = 10, min_events = 20,
                           stratify_cancer = TRUE) {
  
  d <- data %>%
    filter(
      !is.na(.data[[feature]]),
      !is.na(OS_TIME),
      !is.na(OS_EVENT),
      !is.na(AGE),
      !is.na(STAGE),
      !is.na(CANCER_TYPE)
    )
  
  mutated_n <- sum(d[[feature]] == 1, na.rm = TRUE)
  unmutated_n <- sum(d[[feature]] == 0, na.rm = TRUE)
  event_n <- sum(d$OS_EVENT == 1, na.rm = TRUE)
  
  if (mutated_n < min_mutated | unmutated_n < min_mutated | event_n < min_events) {
    return(
      tibble(
        model_name = model_name,
        feature = feature,
        feature_label = feature_label,
        term = feature,
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        ph_p_value = NA_real_,
        mutated_n = mutated_n,
        unmutated_n = unmutated_n,
        event_n = event_n,
        model_ok = FALSE,
        note = "Skipped: insufficient mutated/unmutated count or event count"
      )
    )
  }
  
  rhs_terms <- c(feature, "AGE")
  
  if (n_distinct(d$STAGE) > 1) {
    rhs_terms <- c(rhs_terms, "STAGE")
  }
  
  if (stratify_cancer && n_distinct(d$CANCER_TYPE) > 1) {
    rhs_terms <- c(rhs_terms, "strata(CANCER_TYPE)")
  }
  
  fml <- as.formula(
    paste0(
      "Surv(OS_TIME, OS_EVENT) ~ ",
      paste(rhs_terms, collapse = " + ")
    )
  )
  
  fit <- tryCatch(
    coxph(fml, data = d, x = TRUE),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        model_name = model_name,
        feature = feature,
        feature_label = feature_label,
        term = feature,
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_,
        ph_p_value = NA_real_,
        mutated_n = mutated_n,
        unmutated_n = unmutated_n,
        event_n = event_n,
        model_ok = FALSE,
        note = paste("Model failed:", fit$message)
      )
    )
  }
  
  tidy_fit <- broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(term == feature)
  
  ph_test <- tryCatch(
    cox.zph(fit),
    error = function(e) NULL
  )
  
  ph_p <- NA_real_
  
  if (!is.null(ph_test)) {
    ph_table <- as.data.frame(ph_test$table)
    ph_table$term <- rownames(ph_table)
    
    if (feature %in% ph_table$term) {
      ph_p <- ph_table %>%
        filter(term == feature) %>%
        pull(p)
    }
  }
  
  tidy_fit %>%
    mutate(
      model_name = model_name,
      feature = feature,
      feature_label = feature_label,
      ph_p_value = ph_p,
      mutated_n = mutated_n,
      unmutated_n = unmutated_n,
      event_n = event_n,
      model_ok = TRUE,
      note = NA_character_
    ) %>%
    select(
      model_name,
      feature,
      feature_label,
      term,
      estimate,
      conf.low,
      conf.high,
      p.value,
      ph_p_value,
      mutated_n,
      unmutated_n,
      event_n,
      model_ok,
      note
    )
}

gene_cox_results <- purrr::map_dfr(gene_vars, function(g) {
  fit_binary_cox(
    data = df2,
    feature = g,
    model_name = "Gene-level adjusted Cox model",
    feature_label = g,
    stratify_cancer = TRUE
  )
})

gene_cox_results_clean <- gene_cox_results %>%
  mutate(
    HR = round(estimate, 3),
    CI = paste0(round(conf.low, 3), " - ", round(conf.high, 3)),
    p = format_p(p.value),
    ph_p = format_p(ph_p_value)
  )

print(gene_cox_results_clean)

readr::write_csv(
  gene_cox_results_clean,
  file.path(results_dir, "gene_cox_results.csv")
)

p_gene_cox <- gene_cox_results %>%
  filter(model_ok) %>%
  mutate(
    feature_label = forcats::fct_reorder(feature_label, estimate)
  ) %>%
  ggplot(aes(y = feature_label, x = estimate)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(
    aes(x = conf.low, xend = conf.high, y = feature_label, yend = feature_label),
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(
    title = "Adjusted Gene-Level Cox Models",
    subtitle = "Adjusted for age and stage; stratified by cancer type",
    x = "Hazard Ratio, log scale",
    y = NULL
  ) +
  theme_bw(base_size = 12)

print(p_gene_cox)

ggsave(
  filename = file.path(fig_dir, "phase2_gene_cox_results.png"),
  plot = p_gene_cox,
  width = 7,
  height = 5,
  dpi = 300
)

# RESULT 3: Within-cancer-type mutation effects

top_cancers <- df2 %>%
  count(CANCER_TYPE, sort = TRUE) %>%
  filter(n >= 100) %>%
  slice_head(n = 5)

print(top_cancers)

fit_gene_within_cancer <- function(gene, this_cancer) {
  
  d <- df2 %>%
    filter(CANCER_TYPE == this_cancer)
  
  fit_binary_cox(
    data = d,
    feature = gene,
    model_name = paste("Within-cancer Cox:", this_cancer),
    feature_label = gene,
    min_mutated = 15,
    min_events = 20,
    stratify_cancer = FALSE
  ) %>%
    mutate(cancer_type = as.character(this_cancer))
}

within_cancer_results <- purrr::map_dfr(top_cancers$CANCER_TYPE, function(ct) {
  purrr::map_dfr(gene_vars, function(g) {
    fit_gene_within_cancer(gene = g, this_cancer = ct)
  })
})

within_cancer_results_clean <- within_cancer_results %>%
  mutate(
    HR = round(estimate, 3),
    CI = paste0(round(conf.low, 3), " - ", round(conf.high, 3)),
    p = format_p(p.value),
    ph_p = format_p(ph_p_value)
  )

print(within_cancer_results_clean)

readr::write_csv(
  within_cancer_results_clean,
  file.path(results_dir, "within_cancer_gene_cox_results.csv")
)

heat_df <- within_cancer_results %>%
  filter(model_ok) %>%
  mutate(
    log2_HR = log2(estimate),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_within_cancer_heatmap <- ggplot(
  heat_df,
  aes(x = feature_label, y = cancer_type, fill = log2_HR)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 4) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    title = "Cancer-Type-Specific Mutation Effects on Survival",
    subtitle = "Cells show log2(HR); stars mark nominal significance",
    x = NULL,
    y = NULL,
    fill = "log2(HR)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_within_cancer_heatmap)

ggsave(
  filename = file.path(fig_dir, "phase2_within_cancer_heatmap.png"),
  plot = p_within_cancer_heatmap,
  width = 8,
  height = 5,
  dpi = 300
)

# RESULT 4: Pathway-level mutation groups

add_any_mut <- function(data, genes, new_name) {
  present_genes <- intersect(genes, names(data))
  
  if (length(present_genes) == 0) {
    data[[new_name]] <- NA_integer_
  } else {
    data[[new_name]] <- as.integer(
      rowSums(data[, present_genes, drop = FALSE], na.rm = TRUE) > 0
    )
  }
  
  data
}

df2 <- df2 %>%
  mutate(
    pathway_p53 = if ("TP53" %in% names(.)) TP53 else NA_integer_,
    pathway_wnt = if ("APC" %in% names(.)) APC else NA_integer_,
    pathway_splicing = if ("RBM10" %in% names(.)) RBM10 else NA_integer_
  )

df2 <- add_any_mut(df2, c("KRAS", "BRAF", "EGFR"), "pathway_mapk")
df2 <- add_any_mut(df2, c("PTEN", "PIK3CA"), "pathway_pi3k")

pathway_info <- tibble(
  pathway_var = c(
    "pathway_p53",
    "pathway_mapk",
    "pathway_pi3k",
    "pathway_wnt",
    "pathway_splicing"
  ),
  pathway_label = c(
    "Selected p53 / genomic instability alteration",
    "Selected MAPK-related alteration",
    "Selected PI3K-AKT-related alteration",
    "Selected Wnt-related alteration",
    "Selected RNA-splicing-related alteration"
  ),
  genes_represented = c(
    "TP53",
    "KRAS, BRAF, EGFR",
    "PTEN, PIK3CA",
    "APC",
    "RBM10"
  )
)

pathway_vars <- pathway_info$pathway_var

pathway_cox_results <- purrr::map_dfr(seq_len(nrow(pathway_info)), function(i) {
  fit_binary_cox(
    data = df2,
    feature = pathway_info$pathway_var[i],
    model_name = "Pathway-level adjusted Cox model",
    feature_label = pathway_info$pathway_label[i],
    stratify_cancer = TRUE
  ) %>%
    mutate(
      genes_represented = pathway_info$genes_represented[i]
    )
})

pathway_cox_results_clean <- pathway_cox_results %>%
  mutate(
    HR = round(estimate, 3),
    CI = paste0(round(conf.low, 3), " - ", round(conf.high, 3)),
    p = format_p(p.value),
    ph_p = format_p(ph_p_value)
  )

print(pathway_cox_results_clean)

readr::write_csv(
  pathway_cox_results_clean,
  file.path(results_dir, "pathway_cox_results.csv")
)

p_pathway_cox <- pathway_cox_results %>%
  filter(model_ok) %>%
  mutate(
    feature_label = forcats::fct_reorder(feature_label, estimate)
  ) %>%
  ggplot(aes(y = feature_label, x = estimate)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(
    aes(x = conf.low, xend = conf.high, y = feature_label, yend = feature_label),
    linewidth = 0.7
  ) +
  geom_point(size = 3, color = "darkgreen") +
  scale_x_log10() +
  labs(
    title = "Adjusted Pathway-Level Cox Models",
    subtitle = "Pathway groups are based on selected representative genes",
    x = "Hazard Ratio, log scale",
    y = NULL
  ) +
  theme_bw(base_size = 12)

print(p_pathway_cox)

ggsave(
  filename = file.path(fig_dir, "phase2_pathway_cox_results.png"),
  plot = p_pathway_cox,
  width = 8,
  height = 5,
  dpi = 300
)


# RESULT 5: Driver mutation burden


df2 <- df2 %>%
  mutate(
    mut_count = rowSums(across(all_of(gene_vars)), na.rm = TRUE),
    mut_group = case_when(
      mut_count == 0 ~ "0 mutations",
      mut_count == 1 ~ "1 mutation",
      mut_count >= 2 ~ "2+ mutations",
      TRUE ~ NA_character_
    ),
    mut_group = factor(
      mut_group,
      levels = c("0 mutations", "1 mutation", "2+ mutations")
    )
  ) %>%
  filter(!is.na(mut_group))

mutation_burden_table <- df2 %>%
  count(mut_group) %>%
  mutate(percent = round(n / sum(n) * 100, 1))

print(mutation_burden_table)

readr::write_csv(
  mutation_burden_table,
  file.path(results_dir, "mutation_burden_distribution.csv")
)

# 5A. Kaplan-Meier curve

fit_km_burden <- survfit(
  Surv(OS_TIME, OS_EVENT) ~ mut_group,
  data = df2
)

km_summary <- summary(fit_km_burden)

km_df <- tibble(
  time = km_summary$time,
  survival = km_summary$surv,
  n_risk = km_summary$n.risk,
  n_event = km_summary$n.event,
  strata = km_summary$strata
) %>%
  mutate(
    mut_group = stringr::str_remove(strata, "mut_group=")
  )

logrank_burden <- survdiff(
  Surv(OS_TIME, OS_EVENT) ~ mut_group,
  data = df2
)

logrank_p <- 1 - pchisq(
  logrank_burden$chisq,
  df = length(logrank_burden$n) - 1
)

p_km_burden <- ggplot(
  km_df,
  aes(x = time, y = survival, color = mut_group)
) +
  geom_step(linewidth = 1) +
  annotate(
    "text",
    x = max(km_df$time, na.rm = TRUE) * 0.65,
    y = 0.15,
    label = paste0("Log-rank p = ", format_p(logrank_p)),
    hjust = 0
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Overall Survival by Selected Driver Mutation Burden",
    x = "Overall survival time",
    y = "Survival probability",
    color = "Mutation burden"
  ) +
  theme_bw(base_size = 12)

print(p_km_burden)

ggsave(
  filename = file.path(fig_dir, "phase2_mutation_burden_km.png"),
  plot = p_km_burden,
  width = 7,
  height = 5,
  dpi = 300
)


# 5B. Adjusted Cox model for mutation burden

fit_burden_cox <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ mut_group + AGE + STAGE + strata(CANCER_TYPE),
  data = df2,
  x = TRUE
)

burden_cox_results <- broom::tidy(
  fit_burden_cox,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(stringr::str_detect(term, "^mut_group")) %>%
  mutate(
    comparison = case_when(
      term == "mut_group1 mutation" ~ "1 mutation vs 0 mutations",
      term == "mut_group2+ mutations" ~ "2+ mutations vs 0 mutations",
      TRUE ~ term
    ),
    HR = round(estimate, 3),
    CI = paste0(round(conf.low, 3), " - ", round(conf.high, 3)),
    p = format_p(p.value)
  )

print(burden_cox_results)

readr::write_csv(
  burden_cox_results,
  file.path(results_dir, "mutation_burden_adjusted_cox_results.csv")
)

p_burden_cox <- burden_cox_results %>%
  mutate(
    comparison = forcats::fct_reorder(comparison, estimate)
  ) %>%
  ggplot(aes(y = comparison, x = estimate)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(
    aes(x = conf.low, xend = conf.high, y = comparison, yend = comparison),
    linewidth = 0.7
  ) +
  geom_point(size = 3, color = "purple") +
  scale_x_log10() +
  labs(
    title = "Adjusted Cox Model for Driver Mutation Burden",
    subtitle = "Adjusted for age and stage; stratified by cancer type",
    x = "Hazard Ratio, log scale",
    y = NULL
  ) +
  theme_bw(base_size = 12)

print(p_burden_cox)

ggsave(
  filename = file.path(fig_dir, "phase2_mutation_burden_cox.png"),
  plot = p_burden_cox,
  width = 7,
  height = 4,
  dpi = 300
)


# RESULT 6: Proportional hazards assumption checks

burden_ph <- cox.zph(fit_burden_cox)

burden_ph_table <- as.data.frame(burden_ph$table) %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  mutate(
    model = "Mutation burden adjusted Cox model",
    p_formatted = format_p(p)
  )

print(burden_ph_table)

readr::write_csv(
  burden_ph_table,
  file.path(results_dir, "mutation_burden_ph_assumption_check.csv")
)

ph_summary <- bind_rows(
  gene_cox_results %>%
    filter(model_ok) %>%
    transmute(
      model = "Gene-level adjusted Cox model",
      feature = feature_label,
      ph_p_value = ph_p_value,
      ph_p = format_p(ph_p_value)
    ),
  pathway_cox_results %>%
    filter(model_ok) %>%
    transmute(
      model = "Pathway-level adjusted Cox model",
      feature = feature_label,
      ph_p_value = ph_p_value,
      ph_p = format_p(ph_p_value)
    )
)

print(ph_summary)

readr::write_csv(
  ph_summary,
  file.path(results_dir, "gene_pathway_ph_assumption_summary.csv")
)



# RESULT 7: Website-ready interpretation table

website_key_results <- bind_rows(
  gene_cox_results %>%
    filter(model_ok) %>%
    transmute(
      analysis = "Gene-level Cox model",
      feature = feature_label,
      HR = round(estimate, 2),
      CI = paste0(round(conf.low, 2), "–", round(conf.high, 2)),
      p_value = format_p(p.value),
      interpretation = case_when(
        estimate > 1 ~ "Higher hazard / worse survival association",
        estimate < 1 ~ "Lower hazard / better survival association",
        TRUE ~ "No difference"
      )
    ),
  pathway_cox_results %>%
    filter(model_ok) %>%
    transmute(
      analysis = "Pathway-level Cox model",
      feature = feature_label,
      HR = round(estimate, 2),
      CI = paste0(round(conf.low, 2), "–", round(conf.high, 2)),
      p_value = format_p(p.value),
      interpretation = case_when(
        estimate > 1 ~ "Higher hazard / worse survival association",
        estimate < 1 ~ "Lower hazard / better survival association",
        TRUE ~ "No difference"
      )
    ),
  burden_cox_results %>%
    transmute(
      analysis = "Mutation burden Cox model",
      feature = comparison,
      HR = round(estimate, 2),
      CI = paste0(round(conf.low, 2), "–", round(conf.high, 2)),
      p_value = format_p(p.value),
      interpretation = case_when(
        estimate > 1 ~ "Higher hazard compared with 0 mutations",
        estimate < 1 ~ "Lower hazard compared with 0 mutations",
        TRUE ~ "No difference"
      )
    )
)

print(website_key_results)

readr::write_csv(
  website_key_results,
  file.path(results_dir, "phase2_website_key_results.csv")
)


# ============================================================
# 8. Final Phase 2 completion message
# ============================================================

cat("\nPhase 2 complete.\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Result tables saved to:", results_dir, "\n")

cat("\nWebsite figures created:\n")
cat("- phase2_gene_prevalence.png\n")
cat("- phase2_gene_cox_results.png\n")
cat("- phase2_within_cancer_heatmap.png\n")
cat("- phase2_pathway_cox_results.png\n")
cat("- phase2_mutation_burden_km.png\n")
cat("- phase2_mutation_burden_cox.png\n")


