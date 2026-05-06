library(tidyverse)

dir.create("results/statistics", recursive = TRUE, showWarnings = FALSE)
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

gene_matrix <- read_csv("data/processed/gene_matrix_with_markers.csv", show_col_types = FALSE)

exclude_cols <- c(
  "genome_id",
  "phenotype",
  "total_AMR_genes",
  "carbapenemase_present",
  "blaKPC_any",
  "blaNDM_any",
  "blaOXA48_any",
  "blaVIM_any",
  "blaIMP_any"
)

gene_cols <- setdiff(colnames(gene_matrix), exclude_cols)

run_fisher <- function(gene_name) {
  
  tab <- table(
    Gene = factor(gene_matrix[[gene_name]], levels = c(1, 0)),
    Phenotype = factor(gene_matrix$phenotype, levels = c("Resistant", "Susceptible"))
  )
  
  test <- fisher.test(tab)
  
  tibble(
    gene = gene_name,
    resistant_present = tab["1", "Resistant"],
    susceptible_present = tab["1", "Susceptible"],
    resistant_absent = tab["0", "Resistant"],
    susceptible_absent = tab["0", "Susceptible"],
    odds_ratio = as.numeric(test$estimate),
    p_value = test$p.value
  )
}

association_results <- map_dfr(gene_cols, run_fisher) |>
  mutate(
    adjusted_p_value = p.adjust(p_value, method = "BH"),
    log2_odds_ratio = log2(odds_ratio),
    minus_log10_adj_p = -log10(adjusted_p_value),
    log2_odds_ratio_plot = pmax(pmin(log2_odds_ratio, 10), -10),
    association_category = case_when(
      adjusted_p_value < 0.05 & odds_ratio > 1 ~ "Resistance-associated",
      adjusted_p_value < 0.05 & odds_ratio < 1 ~ "Susceptibility-associated",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(adjusted_p_value)

write_csv(association_results, "results/statistics/association_results.csv")

plot_data <- association_results |>
  filter(is.finite(minus_log10_adj_p))

p <- ggplot(plot_data, aes(x = log2_odds_ratio_plot, y = minus_log10_adj_p)) +
  geom_point(aes(color = association_category), alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Single-Gene Association Analysis",
    x = "log2(Odds Ratio), capped at ±10",
    y = "-log10(Adjusted p-value)"
  )

ggsave("plots/association_volcano_plot.png", p, width = 7, height = 5, dpi = 300)

cat("\nSingle-gene association analysis completed.\n")
cat("Top 15 genes:\n")
print(head(association_results, 15))