library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

gene_matrix <- read_csv("data/processed/gene_matrix_with_markers.csv", show_col_types = FALSE)

marker_cols <- c(
  "genome_id", "phenotype",
  "total_AMR_genes", "carbapenemase_present",
  "blaKPC_any", "blaNDM_any", "blaOXA48_any", "blaVIM_any", "blaIMP_any"
)

gene_cols <- setdiff(colnames(gene_matrix), marker_cols)

sample_summary <- gene_matrix |>
  count(phenotype, name = "number_of_genomes")

gene_frequency <- gene_matrix |>
  summarise(across(all_of(gene_cols), sum)) |>
  pivot_longer(everything(), names_to = "gene", values_to = "frequency") |>
  arrange(desc(frequency))

write_csv(sample_summary, "results/tables/sample_summary.csv")
write_csv(gene_frequency, "results/tables/gene_frequency_summary.csv")

p1 <- ggplot(sample_summary, aes(x = phenotype, y = number_of_genomes)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Meropenem Phenotype Distribution",
    x = "Phenotype",
    y = "Number of Genomes"
  )

ggsave("plots/phenotype_distribution.png", p1, width = 6, height = 4, dpi = 300)

top_genes <- gene_frequency |>
  slice_head(n = 20)

p2 <- ggplot(top_genes, aes(x = reorder(gene, frequency), y = frequency)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top AMR Gene Frequencies",
    x = "AMR Gene",
    y = "Number of Genomes"
  )

ggsave("plots/top_amr_gene_frequencies.png", p2, width = 8, height = 6, dpi = 300)

cat("\nDataset summary completed.\n")
print(sample_summary)
print(head(gene_frequency, 10))