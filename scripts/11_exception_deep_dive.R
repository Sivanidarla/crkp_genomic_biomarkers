library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

resistant_no_carb <- read_csv(
  "results/tables/resistant_without_carbapenemase.csv",
  show_col_types = FALSE
)

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

gene_cols <- setdiff(colnames(resistant_no_carb), exclude_cols)

amr_like_genes <- gene_cols[
  str_detect(
    gene_cols,
    regex(
      "bla|aac|ant|aph|oqx|tet|cat|dfr|sul|erm|qnr|rmt|fos|arr|mph|omp|emr|acr|arn",
      ignore_case = TRUE
    )
  )
]

resistant_no_carb_gene_summary <- resistant_no_carb |>
  summarise(across(all_of(amr_like_genes), sum)) |>
  pivot_longer(everything(), names_to = "gene", values_to = "count") |>
  arrange(desc(count))

write_csv(
  resistant_no_carb_gene_summary,
  "results/tables/resistant_without_carbapenemase_amr_gene_summary.csv"
)

cat("\nException deep-dive completed.\n")
cat("Top AMR-like genes in resistant isolates without carbapenemase:\n")
print(head(resistant_no_carb_gene_summary, 30))