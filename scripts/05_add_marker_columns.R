library(tidyverse)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

gene_matrix <- read_csv("data/processed/gene_matrix.csv", show_col_types = FALSE)

gene_cols <- setdiff(colnames(gene_matrix), c("genome_id", "phenotype"))

mk_any <- function(pattern) {
  cols <- colnames(gene_matrix)[
    str_detect(colnames(gene_matrix), regex(pattern, ignore_case = TRUE))
  ]
  
  if (length(cols) == 0) {
    return(rep(FALSE, nrow(gene_matrix)))
  }
  
  apply(gene_matrix[, cols, drop = FALSE] == 1, 1, any)
}

gene_matrix_updated <- gene_matrix |>
  mutate(
    blaKPC_any = mk_any("blaKPC|KPC"),
    blaNDM_any = mk_any("blaNDM|NDM"),
    blaOXA48_any = mk_any("OXA[_-]?48|OXA[_-]?181|OXA[_-]?232|OXA[_-]?244"),
    blaVIM_any = mk_any("blaVIM|VIM"),
    blaIMP_any = mk_any("blaIMP|IMP"),
    
    carbapenemase_present = blaKPC_any | blaNDM_any | blaOXA48_any | blaVIM_any | blaIMP_any,
    
    total_AMR_genes = rowSums(across(all_of(gene_cols)))
  ) |>
  relocate(
    total_AMR_genes,
    carbapenemase_present,
    blaKPC_any,
    blaNDM_any,
    blaOXA48_any,
    blaVIM_any,
    blaIMP_any,
    .after = phenotype
  )

write_csv(gene_matrix_updated, "data/processed/gene_matrix_with_markers.csv")

cat("\nMarker column summary:\n")

print(
  gene_matrix_updated |>
    summarise(
      blaKPC = sum(blaKPC_any),
      blaNDM = sum(blaNDM_any),
      blaOXA48 = sum(blaOXA48_any),
      blaVIM = sum(blaVIM_any),
      blaIMP = sum(blaIMP_any),
      any_carbapenemase = sum(carbapenemase_present)
    )
)

cat("\nCarbapenemase presence by phenotype:\n")
print(
  gene_matrix_updated |>
    count(phenotype, carbapenemase_present)
)

cat("\nMarker-enhanced gene matrix saved.\n")