library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

gene_matrix <- read_csv("data/processed/gene_matrix_with_markers.csv", show_col_types = FALSE)

resistant_without_carbapenemase <- gene_matrix |>
  filter(phenotype == "Resistant", carbapenemase_present == FALSE)

susceptible_with_carbapenemase <- gene_matrix |>
  filter(phenotype == "Susceptible", carbapenemase_present == TRUE)

write_csv(resistant_without_carbapenemase, "results/tables/resistant_without_carbapenemase.csv")
write_csv(susceptible_with_carbapenemase, "results/tables/susceptible_with_carbapenemase.csv")

exception_summary <- tibble(
  exception_type = c(
    "Resistant without carbapenemase",
    "Susceptible with carbapenemase"
  ),
  count = c(
    nrow(resistant_without_carbapenemase),
    nrow(susceptible_with_carbapenemase)
  )
)

write_csv(exception_summary, "results/tables/exception_summary.csv")

p <- ggplot(exception_summary, aes(x = exception_type, y = count)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Genotype–Phenotype Exceptions",
    x = "",
    y = "Number of Genomes"
  )

ggsave("plots/exception_summary_plot.png", p, width = 7, height = 4, dpi = 300)

cat("\nException analysis completed.\n")
print(exception_summary)