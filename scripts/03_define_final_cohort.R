library(tidyverse)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

phenotypes <- read_csv("data/processed/clean_phenotype_cohort.csv", show_col_types = FALSE)
genes <- read_csv("data/processed/amr_genes_final.csv", show_col_types = FALSE)

final_cohort <- phenotypes |>
  semi_join(genes, by = "genome_id")

write_csv(final_cohort, "data/processed/final_analysis_cohort.csv")

cat("\nFinal analysis cohort saved.\n")
cat("Final cohort size:", nrow(final_cohort), "\n")
print(table(final_cohort$phenotype))