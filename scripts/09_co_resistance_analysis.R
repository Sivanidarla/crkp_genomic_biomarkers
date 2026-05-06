library(tidyverse)
library(pheatmap)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/statistics", recursive = TRUE, showWarnings = FALSE)
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

gene_matrix <- read_csv("data/processed/gene_matrix_with_markers.csv", show_col_types = FALSE)

# -----------------------------
# 1. AMR burden comparison
# -----------------------------

burden_summary <- gene_matrix |>
  group_by(phenotype) |>
  summarise(
    mean = mean(total_AMR_genes),
    median = median(total_AMR_genes),
    sd = sd(total_AMR_genes),
    .groups = "drop"
  )

burden_test <- wilcox.test(total_AMR_genes ~ phenotype, data = gene_matrix)

write_csv(burden_summary, "results/tables/amr_burden_summary.csv")

write_csv(
  tibble(
    test = "Wilcoxon rank-sum test",
    p_value = burden_test$p.value
  ),
  "results/statistics/amr_burden_test.csv"
)

p1 <- ggplot(gene_matrix, aes(x = phenotype, y = total_AMR_genes)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "AMR Gene Burden by Phenotype",
    x = "Phenotype",
    y = "Total AMR Genes"
  )

ggsave("plots/amr_burden_by_phenotype.png", p1, width = 6, height = 4, dpi = 300)

# -----------------------------
# 2. Biologically filtered co-occurrence
# -----------------------------

marker_cols <- c(
  "genome_id", "phenotype",
  "total_AMR_genes",
  "carbapenemase_present",
  "blaKPC_any", "blaNDM_any",
  "blaOXA48_any", "blaVIM_any", "blaIMP_any"
)

gene_cols <- setdiff(colnames(gene_matrix), marker_cols)

amr_like_genes <- gene_cols[
  str_detect(
    gene_cols,
    regex(
      "bla|aac|ant|aph|oqx|tet|cat|dfr|sul|erm|qnr|rmt|fos|arr|mph",
      ignore_case = TRUE
    )
  )
]

gene_freq <- gene_matrix |>
  summarise(across(all_of(amr_like_genes), mean)) |>
  pivot_longer(everything(), names_to = "gene", values_to = "freq")

informative_genes <- gene_freq |>
  filter(freq > 0.05, freq < 0.95) |>
  arrange(desc(freq)) |>
  slice_head(n = 30) |>
  pull(gene)

mat <- gene_matrix |>
  select(all_of(informative_genes)) |>
  as.matrix()

cooccurrence <- t(mat) %*% mat

cooccurrence_df <- as.data.frame(as.table(cooccurrence)) |>
  rename(gene1 = Var1, gene2 = Var2, count = Freq) |>
  mutate(
    gene1 = as.character(gene1),
    gene2 = as.character(gene2),
    count = as.numeric(count)
  ) |>
  filter(gene1 != gene2) |>
  filter(gene1 < gene2) |>
  arrange(desc(count))

write_csv(cooccurrence_df, "results/tables/gene_cooccurrence_table.csv")

# -----------------------------
# 3. Clustered heatmap
# -----------------------------

heatmap_data <- gene_matrix |>
  select(genome_id, phenotype, all_of(informative_genes))

annotation <- heatmap_data |>
  select(phenotype) |>
  as.data.frame()

rownames(annotation) <- heatmap_data$genome_id

heatmap_matrix <- heatmap_data |>
  select(all_of(informative_genes)) |>
  as.matrix()

rownames(heatmap_matrix) <- heatmap_data$genome_id

png("plots/clustered_amr_heatmap.png", width = 1400, height = 1000, res = 150)

pheatmap(
  heatmap_matrix,
  annotation_row = annotation,
  show_rownames = FALSE,
  fontsize_col = 8,
  main = "Clustered AMR Gene Presence/Absence"
)

dev.off()

cat("\nCo-resistance analysis completed.\n")
cat("\nAMR burden summary:\n")
print(burden_summary)

cat("\nAMR burden Wilcoxon p-value:\n")
print(burden_test$p.value)

cat("\nTop 10 co-occurring AMR-like gene pairs:\n")
print(head(cooccurrence_df, 10))