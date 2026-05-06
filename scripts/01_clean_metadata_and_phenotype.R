library(tidyverse)
library(janitor)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

genome <- read_csv("genome.csv") |> clean_names()
phenotype <- read_csv("phenotype.csv") |> clean_names()

genome_kp_clean <- genome |>
  filter(
    str_detect(str_to_lower(genome_name), "^klebsiella pneumoniae"),
    genome_status == "Complete",
    genome_quality == "Good"
  ) |>
  distinct(genome_id, .keep_all = TRUE)

phenotype_meropenem_rs <- phenotype |>
  filter(
    str_to_lower(antibiotic) == "meropenem",
    resistant_phenotype %in% c("Resistant", "Susceptible")
  )

clean_phenotype_cohort <- phenotype_meropenem_rs |>
  inner_join(genome_kp_clean, by = "genome_id") |>
  group_by(genome_id) |>
  filter(n_distinct(resistant_phenotype) == 1) |>
  slice(1) |>
  ungroup() |>
  transmute(
    genome_id,
    genome_name = genome_name.y,
    antibiotic,
    phenotype = resistant_phenotype
  )

write_csv(genome_kp_clean, "data/processed/clean_genome_metadata.csv")
write_csv(clean_phenotype_cohort, "data/processed/clean_phenotype_cohort.csv")

cat("\nClean metadata and phenotype files saved.\n")
cat("Clean phenotype cohort size:", nrow(clean_phenotype_cohort), "\n")
print(table(clean_phenotype_cohort$phenotype))