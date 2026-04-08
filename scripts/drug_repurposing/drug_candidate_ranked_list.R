library(tidyverse)

# 1. Load PPI network - FIXED VERSION
load_ppi <- function(ppi_path) {
  # Read all columns but only keep needed ones
  ppi <- read_tsv(ppi_path, comment = "#", show_col_types = FALSE,
                  col_names = FALSE)
  
  # Check how many columns your file has
  cat("PPI columns detected:", ncol(ppi), "\n")
  
  # Based on your original script, column names are:
  # gene1, gene2, string_id1, string_id2, neighborhood, gene_fusion, 
  # phylogenetic_cooccurrence, homology, coexpression, experimental, 
  # database, textmining, combined_score
  
  # Extract just gene names and combined_score (should be column 13)
  if (ncol(ppi) >= 13) {
    ppi <- ppi %>%
      select(gene1 = X1, gene2 = X2, combined_score = X13) %>%
      filter(!is.na(gene1), !is.na(gene2), !is.na(combined_score))
  } else if (ncol(ppi) == 3) {
    # If file only has 3 columns, use as-is
    ppi <- ppi %>%
      rename(gene1 = X1, gene2 = X2, combined_score = X3) %>%
      filter(!is.na(gene1), !is.na(gene2), !is.na(combined_score))
  } else {
    stop("PPI file has unexpected number of columns: ", ncol(ppi))
  }
  
  # Convert combined_score to numeric if needed
  ppi <- ppi %>%
    mutate(combined_score = as.numeric(combined_score))
  
  return(ppi)
}

# 2. Load drug-gene interactions - IMPROVED VERSION
load_drugs <- function(drug_path) {
  # Try to read with different column names
  drug_data <- read_tsv(drug_path, show_col_types = FALSE)
  
  # Check what columns we have
  cat("Drug file columns:", paste(colnames(drug_data), collapse = ", "), "\n")
  
  # Find drug name and gene name columns (case-insensitive)
  col_names_lower <- tolower(colnames(drug_data))
  
  drug_col <- ifelse("drug_name" %in% colnames(drug_data), "drug_name",
                    ifelse("drug_name" %in% col_names_lower, 
                           colnames(drug_data)[which(col_names_lower == "drug_name")[1]],
                    ifelse("drug" %in% col_names_lower,
                           colnames(drug_data)[which(col_names_lower == "drug")[1]],
                           NA_character_)))
  
  gene_col <- ifelse("gene_name" %in% colnames(drug_data), "gene_name",
                    ifelse("gene_name" %in% col_names_lower,
                           colnames(drug_data)[which(col_names_lower == "gene_name")[1]],
                    ifelse("gene" %in% col_names_lower,
                           colnames(drug_data)[which(col_names_lower == "gene")[1]],
                           NA_character_)))
  
  if (is.na(drug_col) || is.na(gene_col)) {
    stop("Could not find drug_name and gene_name columns in drug file")
  }
  
  cat("Using columns:", drug_col, "and", gene_col, "\n")
  
  # Clean and select
  drug_data %>%
    rename(drug_name = all_of(drug_col),
           gene_name = all_of(gene_col)) %>%
    select(drug_name, gene_name) %>%
    mutate(
      drug_name = as.character(str_trim(drug_name)),
      gene_name = as.character(str_trim(gene_name))
    ) %>%
    filter(
      !is.na(drug_name), drug_name != "",
      !is.na(gene_name), gene_name != "",
      !str_detect(tolower(drug_name), "^null$|^na$|^$")
    ) %>%
    distinct(drug_name, gene_name, .keep_all = TRUE)
}

# 3. Calculate network centrality for genes - FIXED
get_gene_centrality <- function(ppi) {
  # Degree centrality
  degree_counts <- c(ppi$gene1, ppi$gene2) %>%
    table() %>%
    as.data.frame() %>%
    rename(gene = 1, degree = 2)
  
  # Average connection strength
  gene_strength <- ppi %>%
    pivot_longer(cols = c(gene1, gene2),
                 names_to = "position",
                 values_to = "gene") %>%
    group_by(gene) %>%
    summarise(avg_score = mean(combined_score, na.rm = TRUE),
              .groups = "drop")
  
  # Combine metrics
  degree_counts %>%
    inner_join(gene_strength, by = "gene") %>%
    mutate(degree_norm = degree / max(degree, na.rm = TRUE),
           score_norm = avg_score / max(avg_score, na.rm = TRUE))
}

# 4. Score drugs based on target network properties
score_drugs <- function(drugs, gene_scores, min_targets = 2) {
  drugs %>%
    inner_join(gene_scores, by = c("gene_name" = "gene")) %>%
    group_by(drug_name) %>%
    summarise(
      n_targets = n(),
      avg_degree = mean(degree_norm, na.rm = TRUE),
      avg_strength = mean(score_norm, na.rm = TRUE),
      target_genes = paste(sort(unique(gene_name)), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(n_targets >= min_targets) %>%
    mutate(
      # Handle cases where scores might be NaN
      avg_degree = ifelse(is.na(avg_degree), 0, avg_degree),
      avg_strength = ifelse(is.na(avg_strength), 0, avg_strength),
      
      # Simple weighted score
      network_score = (avg_degree * 0.6 + avg_strength * 0.4) * log10(n_targets + 1)
    ) %>%
    arrange(desc(network_score)) %>%
    select(drug_name, network_score, n_targets, avg_degree, avg_strength, target_genes)
}

# 5. Main function
main <- function(ppi_path, drug_path, output_prefix = "als_drugs") {
  cat("Loading PPI network...\n")
  ppi <- load_ppi(ppi_path)
  cat("PPI edges:", nrow(ppi), "\n")
  cat("Unique genes in PPI:", n_distinct(c(ppi$gene1, ppi$gene2)), "\n")
  
  cat("\nLoading drug-gene interactions...\n")
  drugs <- load_drugs(drug_path)
  cat("Drug-gene pairs:", nrow(drugs), "\n")
  cat("Unique drugs:", n_distinct(drugs$drug_name), "\n")
  cat("Unique genes in drug data:", n_distinct(drugs$gene_name), "\n")
  
  cat("\nCalculating gene centrality...\n")
  gene_scores <- get_gene_centrality(ppi)
  cat("Genes with centrality scores:", nrow(gene_scores), "\n")
  
  cat("\nScoring drugs...\n")
  ranked_drugs <- score_drugs(drugs, gene_scores)
  
  cat("\n=== TOP 20 DRUG CANDIDATES ===\n")
  print(head(ranked_drugs, 20))
  
  # Save results
  output_file <- paste0(output_prefix, "_ranked.csv")
  write_csv(ranked_drugs, output_file)
  write_csv(gene_scores, paste0(output_prefix, "_gene_scores.csv"))
  
  cat("\nResults saved to:", output_file, "\n")
  cat("Total drugs ranked:", nrow(ranked_drugs), "\n")
  
  # Summary stats
  cat("\n=== SUMMARY ===\n")
  cat("Network covers", nrow(gene_scores), "genes\n")
  cat("Drugs with at least 2 targets in network:", nrow(ranked_drugs), "\n")
  cat("Average targets per drug:", mean(ranked_drugs$n_targets), "\n")
  cat("Max network score:", max(ranked_drugs$network_score, na.rm = TRUE), "\n")
  cat("Min network score:", min(ranked_drugs$network_score, na.rm = TRUE), "\n")
  
  return(list(ranked_drugs = ranked_drugs, gene_scores = gene_scores))
}

# 6. Run analysis
results <- main(
  ppi_path = "C:/Users/noahm/PycharmProjects/MarbleProject/ppi_results/results_moderate/filtered_ppi.tsv",
  drug_path = "C:/Users/noahm/PycharmProjects/MarbleProject/datasets/drug_targets/interactions.tsv",
  output_prefix = "moderate"
)