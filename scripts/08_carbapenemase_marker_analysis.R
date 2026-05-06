library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

gene_matrix <- read_csv("data/processed/gene_matrix_with_markers.csv", show_col_types = FALSE)

markers <- c(
  "blaKPC_any",
  "blaNDM_any",
  "blaOXA48_any",
  "blaVIM_any",
  "blaIMP_any"
)

carbapenemase_frequency <- gene_matrix |>
  select(phenotype, all_of(markers)) |>
  pivot_longer(
    cols = all_of(markers),
    names_to = "marker",
    values_to = "present"
  ) |>
  group_by(phenotype, marker) |>
  summarise(
    count_present = sum(present),
    total = n(),
    frequency = count_present / total,
    .groups = "drop"
  )

carbapenemase_summary <- gene_matrix |>
  group_by(phenotype) |>
  summarise(
    total_genomes = n(),
    with_carbapenemase = sum(carbapenemase_present),
    without_carbapenemase = sum(!carbapenemase_present),
    percentage_with = with_carbapenemase / total_genomes * 100,
    .groups = "drop"
  )

write_csv(carbapenemase_frequency, "results/tables/carbapenemase_frequency_by_phenotype.csv")
write_csv(carbapenemase_summary, "results/tables/carbapenemase_summary.csv")

p <- ggplot(carbapenemase_frequency, aes(x = marker, y = frequency, fill = phenotype)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Carbapenemase Marker Frequency by Phenotype",
    x = "Marker",
    y = "Frequency"
  )

ggsave("plots/carbapenemase_frequency_by_phenotype.png", p, width = 8, height = 5, dpi = 300)

cat("\nCarbapenemase marker analysis completed.\n")
print(carbapenemase_summary)

cat("\nTop marker frequencies:\n")
print(
  carbapenemase_frequency |>
    arrange(desc(frequency)) |>
    head(10)
)