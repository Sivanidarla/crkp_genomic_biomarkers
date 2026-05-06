library(tidyverse)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

phenotypes <- read_csv("data/processed/final_analysis_cohort.csv", show_col_types = FALSE)
genes <- read_csv("data/processed/amr_genes_final.csv", show_col_types = FALSE)

genes_filtered <- genes |>
  filter(genome_id %in% phenotypes$genome_id)

gene_matrix <- genes_filtered |>
  mutate(present = 1) |>
  select(genome_id, gene, present) |>
  pivot_wider(
    names_from = gene,
    values_from = present,
    values_fill = 0
  )

final_matrix <- phenotypes |>
  select(genome_id, phenotype) |>
  left_join(gene_matrix, by = "genome_id") |>
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

write_csv(final_matrix, "data/processed/gene_matrix.csv")

cat("\nGene matrix saved.\n")
cat("Dimensions:\n")
print(dim(final_matrix))
cat("\nPhenotype distribution:\n")
print(table(final_matrix$phenotype))