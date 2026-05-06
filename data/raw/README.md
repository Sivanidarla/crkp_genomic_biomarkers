Raw data were obtained from the BV-BRC database.

A hybrid acquisition approach was used due to dataset size and platform constraints:

- Genome metadata were downloaded manually after filtering for:
  - Organism: *Klebsiella pneumoniae*
  - Genome status: Complete genomes
  - Genome quality: High-quality assemblies

- Phenotype data were downloaded manually after filtering for:
  - Antibiotic: meropenem
  - Phenotype: Resistant or Susceptible (intermediate and missing values excluded)

- AMR gene annotations were downloaded separately and processed in R due to file size and formatting constraints.