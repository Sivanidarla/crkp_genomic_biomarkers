library(tidyverse)
library(janitor)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

genes <- read_csv("data/raw/amr_genes.csv") |> clean_names()

amr_genes_clean <- genes |>
  filter(!is.na(gene), gene != "") |>
  mutate(
    gene = str_trim(gene),
    gene = str_replace_all(gene, "[^A-Za-z0-9_]", "_")
  ) |>
  distinct()

amr_genes_final <- amr_genes_clean |>
  distinct(genome_id, gene, .keep_all = TRUE)

write_csv(amr_genes_final, "data/processed/amr_genes_final.csv")

cat("\nAMR genes cleaned and saved.\n")
cat("Final AMR gene rows:", nrow(amr_genes_final), "\n")