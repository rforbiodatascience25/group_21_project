CPM_normalization <- function(raw_counts){
  CPM_normalized_counts <- raw_counts |>
    mutate(across(!GeneFeature, ~(.*1000000)/sum(.)))

  return(CPM_normalized_counts)
}

quantile_normalization <- function(raw_counts){
  raw_counts_ranks <- raw_counts |>
    mutate(across(!GeneFeature, ~rank(., ties.method = "min")))

  raw_counts_sorted_cols <- raw_counts |>
    mutate(across(!GeneFeature, sort))

  raw_counts_means <- raw_counts_sorted_cols |>
    rowwise() |>
    summarise(m = mean(c_across(!GeneFeature))) |>
    ungroup() |>
    rowid_to_column("ID")

  lookup_vec <- setNames(raw_counts_means$m, raw_counts_means$ID)
  quantile_normalized_counts <- raw_counts_ranks |>
    mutate(across(!GeneFeature, ~ lookup_vec[.] |> as.numeric()))

  return(quantile_normalized_counts)
}

#Method from: https://pubmed.ncbi.nlm.nih.gov/28650338/ (including quantile normalization and log10 transformation)
#Gene weights from: https://www.science.org/doi/10.1126/science.aar3593 (Supplementary table B)
#Returns sample metadata with augmented GEP scores
GEP_score_calculation <- function(quantile_normalized_counts, metadata){
  signature_genes <- tibble(
    gene = c("CCL5", "CD27", "CD274", "CD276","CD8A","CMKLR1","CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1",
               "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
  )

  augmented_metadata <- quantile_normalized_counts |>
    mutate(across(!GeneFeature, ~log10(.+1))) |>
    right_join(signature_genes,
               join_by(GeneFeature == gene)) |>
    summarise(across(!GeneFeature, mean)) |>
    pivot_longer(everything(),
                 names_to = "Sample_names",
                 values_to = "GEP_score") |>
    right_join(metadata,
               join_by(Sample_names))

  return(augmented_metadata)
}

run_fisher_test <- function(data, var1, var2) {
  # var1, var2: column names as strings

  tab <- table(data[[var1]], data[[var2]])
  fisher.test(tab)
}

