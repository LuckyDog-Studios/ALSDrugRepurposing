# Load required libraries
library(tidyverse)
library(igraph)
library(caret)
library(randomForest)
library(xgboost)
library(doParallel)

# Enable parallel processing
registerDoParallel(cores = parallel::detectCores() - 1)

# ==================== CONFIGURATION ====================
ppi_file <- "C:/Users/noahm/PycharmProjects/MarbleProject/ppi_results/results_lenient/filtered_ppi.tsv"
drug_files <- c(
  "C:/Users/noahm/PycharmProjects/MarbleProject/datasets/drug_targets/interactions.tsv"
)
output_dir <- "drug_candidate_results/results_lenient/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==================== DATA LOADING ====================

load_ppi_network <- function(file_path) {
  cat("Loading PPI network from:", file_path, "\n")
  
  # First, check what columns exist
  first_line <- readLines(file_path, n = 1)
  cat("First line of PPI file:", first_line, "\n")
  
  # Read PPI data with proper column specification
  ppi_data <- tryCatch({
    read_tsv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    cat("Error reading as TSV, trying with delim\n")
    read_delim(file_path, delim = "\t", show_col_types = FALSE)
  })
  
  # Show column names
  cat("Columns in PPI data:\n")
  print(colnames(ppi_data))
  cat("\nFirst few rows:\n")
  print(head(ppi_data))
  
  # Check which column names contain "node1" (case-insensitive)
  node1_col <- grep("node1", colnames(ppi_data), ignore.case = TRUE, value = TRUE)[1]
  node2_col <- grep("node2", colnames(ppi_data), ignore.case = TRUE, value = TRUE)[1]
  
  if (is.na(node1_col) || is.na(node2_col)) {
    # Try to find gene/protein columns
    possible_cols <- colnames(ppi_data)
    node1_col <- possible_cols[1]
    node2_col <- possible_cols[2]
    cat("Using columns", node1_col, "and", node2_col, "as node identifiers\n")
  }
  
  # Filter out NA values
  ppi_data <- ppi_data %>%
    filter(!is.na(!!sym(node1_col)) & !is.na(!!sym(node2_col)))
  
  # Create igraph object
  g <- graph_from_data_frame(ppi_data %>% select(all_of(c(node1_col, node2_col))), directed = FALSE)
  
  # Store edge attributes if they exist
  if ("combined_score" %in% colnames(ppi_data)) {
    edge_attr(g, "combined_score") <- ppi_data$combined_score
  }
  
  if ("experimentally_determined_interaction" %in% colnames(ppi_data)) {
    edge_attr(g, "experimental_score") <- ppi_data$experimentally_determined_interaction
  }
  
  cat("Network loaded:", vcount(g), "nodes,", ecount(g), "edges\n")
  return(g)
}

load_drug_interactions <- function(file_paths) {
  cat("Loading drug interactions...\n")
  
  drug_data <- list()
  for (file_path in file_paths) {
    cat("  Loading:", file_path, "\n")
    
    # Check file structure first
    test_data <- readLines(file_path, n = 5)
    cat("First few lines of drug file:\n")
    cat(paste(test_data, collapse = "\n"), "\n\n")
    
    data <- read_tsv(file_path, show_col_types = FALSE)
    
    # Show column names
    cat("Columns in drug data:\n")
    print(colnames(data))
    cat("\n")
    
    # Clean column names
    data <- data %>%
      rename_with(~ tolower(gsub(" ", "_", .x)))
    
    # Find the relevant columns
    gene_col <- NULL
    drug_col <- NULL
    
    # Try to find gene column
    possible_gene_cols <- c("gene_name", "gene_claim_name", "gene_concept_id", "gene")
    for (col in possible_gene_cols) {
      if (col %in% colnames(data)) {
        gene_col <- col
        break
      }
    }
    
    # Try to find drug column
    possible_drug_cols <- c("drug_name", "drug_claim_name", "drug_concept_id", "drug")
    for (col in possible_drug_cols) {
      if (col %in% colnames(data)) {
        drug_col <- col
        break
      }
    }
    
    if (is.null(gene_col) || is.null(drug_col)) {
      cat("Warning: Could not find standard columns in", file_path, "\n")
      cat("Available columns:", paste(colnames(data), collapse = ", "), "\n")
      next
    }
    
    # Select relevant columns
    data <- data %>%
      select(
        gene_name = !!sym(gene_col),
        drug_name = !!sym(drug_col),
        any_of(c("interaction_type", "interaction_score", "approved", 
                "anti_neoplastic", "immunotherapy"))
      ) %>%
      filter(!is.na(gene_name) & !is.na(drug_name))
    
    drug_data[[file_path]] <- data
  }
  
  if (length(drug_data) == 0) {
    stop("No drug data could be loaded")
  }
  
  # Combine all drug data
  combined_drugs <- bind_rows(drug_data) %>%
    distinct() %>%
    mutate(
      gene_name = as.character(gene_name),
      drug_name = as.character(drug_name)
    )
  
  # Convert numeric columns if they exist
  if ("interaction_score" %in% colnames(combined_drugs)) {
    combined_drugs <- combined_drugs %>%
      mutate(interaction_score = as.numeric(interaction_score))
  }
  
  # Convert logical columns if they exist
  logical_cols <- c("approved", "anti_neoplastic", "immunotherapy")
  for (col in logical_cols) {
    if (col %in% colnames(combined_drugs)) {
      combined_drugs <- combined_drugs %>%
        mutate(!!col := as.logical(tolower(!!sym(col)) %in% 
                                   c("true", "t", "yes", "y", "1")))
    }
  }
  
  cat("Drug interactions loaded:", nrow(combined_drugs), "interactions\n")
  cat("Unique genes:", n_distinct(combined_drugs$gene_name), "\n")
  cat("Unique drugs:", n_distinct(combined_drugs$drug_name), "\n")
  
  return(combined_drugs)
}

# ==================== NETWORK ANALYSIS ====================

calculate_network_properties <- function(g) {
  cat("Calculating network properties...\n")
  
  # Calculate basic properties first
  node_names <- V(g)$name
  
  if (length(node_names) == 0) {
    stop("Graph has no vertices")
  }
  
  # Initialize data frame
  node_properties <- data.frame(
    gene = node_names,
    stringsAsFactors = FALSE
  )
  
  # Calculate degree (always works)
  node_properties$degree <- degree(g)
  
  # Calculate betweenness (might fail for disconnected graphs)
  tryCatch({
    node_properties$betweenness <- betweenness(g, normalized = TRUE)
  }, error = function(e) {
    cat("Warning: Could not calculate betweenness:", e$message, "\n")
    node_properties$betweenness <- NA
  })
  
  # Calculate closeness (might fail for disconnected graphs)
  tryCatch({
    node_properties$closeness <- closeness(g, normalized = TRUE)
  }, error = function(e) {
    cat("Warning: Could not calculate closeness:", e$message, "\n")
    node_properties$closeness <- NA
  })
  
  # Calculate eigenvector centrality
  tryCatch({
    node_properties$eigenvector <- eigen_centrality(g)$vector
  }, error = function(e) {
    cat("Warning: Could not calculate eigenvector centrality:", e$message, "\n")
    node_properties$eigenvector <- NA
  })
  
  # Calculate PageRank
  tryCatch({
    node_properties$page_rank <- page_rank(g)$vector
  }, error = function(e) {
    cat("Warning: Could not calculate PageRank:", e$message, "\n")
    node_properties$page_rank <- NA
  })
  
  # Calculate clustering coefficient
  tryCatch({
    node_properties$clustering_coef <- transitivity(g, type = "local", isolates = "zero")
  }, error = function(e) {
    cat("Warning: Could not calculate clustering coefficient:", e$message, "\n")
    node_properties$clustering_coef <- NA
  })
  
  # Calculate community structure (simplified)
  tryCatch({
    communities <- cluster_louvain(g)
    node_properties$community <- membership(communities)
  }, error = function(e) {
    cat("Warning: Could not calculate communities:", e$message, "\n")
    node_properties$community <- 1  # Assign all to same community
  })
  
  # Calculate edge-based metrics
  node_properties$avg_edge_score <- sapply(node_names, function(node) {
    edges <- incident(g, node)
    if (length(edges) > 0 && "combined_score" %in% list.edge.attributes(g)) {
      scores <- edge_attr(g, "combined_score", edges)
      mean(scores, na.rm = TRUE)
    } else {
      0
    }
  })
  
  # Calculate neighborhood connectivity
  node_properties$neighbor_connectivity <- sapply(node_names, function(node) {
    neighbors <- neighbors(g, node)$name
    if (length(neighbors) > 0) {
      mean(degree(g, neighbors), na.rm = TRUE)
    } else {
      0
    }
  })
  
  # Identify hub genes (top 10% by degree)
  if (all(!is.na(node_properties$degree))) {
    degree_threshold <- quantile(node_properties$degree, 0.9, na.rm = TRUE)
    node_properties$is_hub <- node_properties$degree >= degree_threshold
  } else {
    node_properties$is_hub <- FALSE
  }
  
  # Identify bottleneck genes (top 10% by betweenness)
  if (all(!is.na(node_properties$betweenness))) {
    betweenness_threshold <- quantile(node_properties$betweenness, 0.9, na.rm = TRUE)
    node_properties$is_bottleneck <- node_properties$betweenness >= betweenness_threshold
  } else {
    node_properties$is_bottleneck <- FALSE
  }
  
  # Replace NA values with 0 or appropriate defaults
  node_properties <- node_properties %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  cat("Network properties calculated for", nrow(node_properties), "genes\n")
  cat("Hub genes:", sum(node_properties$is_hub), "\n")
  cat("Bottleneck genes:", sum(node_properties$is_bottleneck), "\n")
  
  return(node_properties)
}

# ==================== FEATURE ENGINEERING ====================

create_training_features <- function(node_properties, drug_interactions) {
  cat("Creating training features...\n")
  
  # Check overlap between network genes and drug target genes
  network_genes <- unique(node_properties$gene)
  drug_genes <- unique(drug_interactions$gene_name)
  
  cat("Genes in network:", length(network_genes), "\n")
  cat("Genes targeted by drugs:", length(drug_genes), "\n")
  cat("Overlap:", length(intersect(network_genes, drug_genes)), "\n")
  
  # Merge network properties with drug interactions
  features <- drug_interactions %>%
    left_join(node_properties, by = c("gene_name" = "gene")) %>%
    # Create interaction type features if available
    mutate(
      interaction_type_numeric = if ("interaction_type" %in% colnames(.)) {
        case_when(
          grepl("inhibitor|antagonist", tolower(interaction_type)) ~ -1,
          grepl("agonist|activator", tolower(interaction_type)) ~ 1,
          TRUE ~ 0
        )
      } else {
        0
      }
    ) %>%
    # Create derived features
    group_by(drug_name) %>%
    mutate(
      num_targets = n_distinct(gene_name),
      avg_target_degree = mean(degree, na.rm = TRUE),
      avg_target_betweenness = mean(betweenness, na.rm = TRUE),
      hub_target_ratio = sum(is_hub, na.rm = TRUE) / max(num_targets, 1),
      bottleneck_target_ratio = sum(is_bottleneck, na.rm = TRUE) / max(num_targets, 1)
    ) %>%
    ungroup() %>%
    # Aggregate features at drug level
    group_by(drug_name) %>%
    summarise(
      # Drug properties if available
      is_approved = if ("approved" %in% colnames(.)) any(approved, na.rm = TRUE) else FALSE,
      is_anti_neoplastic = if ("anti_neoplastic" %in% colnames(.)) any(anti_neoplastic, na.rm = TRUE) else FALSE,
      is_immunotherapy = if ("immunotherapy" %in% colnames(.)) any(immunotherapy, na.rm = TRUE) else FALSE,
      
      # Network topology features
      avg_degree = mean(degree, na.rm = TRUE),
      max_degree = max(degree, na.rm = TRUE),
      avg_betweenness = mean(betweenness, na.rm = TRUE),
      max_betweenness = max(betweenness, na.rm = TRUE),
      avg_closeness = mean(closeness, na.rm = TRUE),
      avg_eigenvector = mean(eigenvector, na.rm = TRUE),
      avg_page_rank = mean(page_rank, na.rm = TRUE),
      avg_clustering_coef = mean(clustering_coef, na.rm = TRUE),
      
      # Community features
      num_communities = n_distinct(community),
      
      # Edge score features
      avg_edge_score = mean(avg_edge_score, na.rm = TRUE),
      
      # Neighborhood features
      avg_neighbor_connectivity = mean(neighbor_connectivity, na.rm = TRUE),
      
      # Target composition features
      num_targets = first(num_targets),
      hub_target_ratio = first(hub_target_ratio),
      bottleneck_target_ratio = first(bottleneck_target_ratio),
      
      # Interaction features if available
      avg_interaction_score = if ("interaction_score" %in% colnames(.)) 
        mean(interaction_score, na.rm = TRUE) else 0,
      avg_interaction_type = mean(interaction_type_numeric, na.rm = TRUE)
    ) %>%
    # Handle NA values
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
    # Create drug ID
    mutate(drug_id = row_number())
  
  cat("Features created for", nrow(features), "drugs\n")
  
  if (nrow(features) == 0) {
    stop("No features could be created. Check if network genes and drug targets overlap.")
  }
  
  return(features)
}

# ==================== SIMPLIFIED DRUG RANKING ====================
# Since we don't have ground truth labels, we'll use unsupervised ranking

rank_drugs_simple <- function(features) {
  cat("Ranking drugs using unsupervised approach...\n")
  
  # Create a composite score based on network properties
  features <- features %>%
    mutate(
      # Normalize key features
      norm_degree = scale(avg_degree)[,1],
      norm_betweenness = scale(avg_betweenness)[,1],
      norm_hub_ratio = scale(hub_target_ratio)[,1],
      norm_num_targets = scale(num_targets)[,1],
      
      # Create composite score (weighted combination)
      network_score = (
        0.3 * ifelse(is.finite(norm_degree), norm_degree, 0) +
        0.3 * ifelse(is.finite(norm_betweenness), norm_betweenness, 0) +
        0.2 * ifelse(is.finite(norm_hub_ratio), norm_hub_ratio, 0) +
        0.2 * ifelse(is.finite(norm_num_targets), norm_num_targets, 0)
      ),
      
      # Apply bonus for approved/anti-neoplastic drugs if available
      bonus = 0.5 * as.numeric(is_approved) + 0.3 * as.numeric(is_anti_neoplastic)
    ) %>%
    mutate(
      final_score = network_score + bonus,
      repurposing_score = (final_score - min(final_score, na.rm = TRUE)) / 
                         (max(final_score, na.rm = TRUE) - min(final_score, na.rm = TRUE))
    )
  
  # Rank drugs
  drug_ranking <- features %>%
    select(drug_name, repurposing_score, 
           avg_degree, avg_betweenness, hub_target_ratio, num_targets,
           is_approved, is_anti_neoplastic) %>%
    mutate(
      rank = rank(-repurposing_score, ties.method = "min"),
      confidence = case_when(
        repurposing_score >= 0.7 ~ "high",
        repurposing_score >= 0.4 ~ "medium",
        TRUE ~ "low"
      )
    ) %>%
    arrange(rank)
  
  cat("Drugs ranked:", nrow(drug_ranking), "\n")
  return(drug_ranking)
}

# ==================== VISUALIZATION ====================

generate_visualizations <- function(g, node_properties, drug_ranking) {
  cat("Generating visualizations...\n")
  
  # 1. Network visualization (simplified)
  tryCatch({
    # Take top 30 nodes by degree for visualization
    if (nrow(node_properties) > 30) {
      top_genes <- node_properties %>%
        arrange(desc(degree)) %>%
        head(30) %>%
        pull(gene)
    } else {
      top_genes <- node_properties$gene
    }
    
    if (length(top_genes) > 1) {
      subgraph <- induced_subgraph(g, top_genes)
      
      png(file.path(output_dir, "network_plot.png"), width = 1000, height = 800)
      plot(subgraph,
           vertex.size = sqrt(degree(subgraph)) * 3,
           vertex.color = ifelse(V(subgraph)$name %in% 
                                node_properties$gene[node_properties$is_hub],
                              "red", "lightblue"),
           vertex.label.cex = 0.7,
           vertex.label.dist = 1,
           edge.width = 1,
           main = "PPI Network (Top 30 Genes)")
      dev.off()
      cat("Network plot saved\n")
    }
  }, error = function(e) {
    cat("Could not generate network plot:", e$message, "\n")
  })
  
  # 2. Top drugs bar plot
  tryCatch({
    top_drugs <- drug_ranking %>%
      head(20) %>%
      mutate(drug_name_short = substr(drug_name, 1, 30))
    
    p <- ggplot(top_drugs, aes(x = reorder(drug_name_short, repurposing_score), 
                              y = repurposing_score, fill = confidence)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = "Drug (truncated)", y = "Repurposing Score", 
           title = "Top 20 Drug Repurposing Candidates") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8)) +
      scale_fill_manual(values = c("high" = "#2E8B57", "medium" = "#FFD700", "low" = "#DC143C"))
    
    ggsave(file.path(output_dir, "top_drugs_plot.png"), p, width = 10, height = 8)
    cat("Top drugs plot saved\n")
  }, error = function(e) {
    cat("Could not generate top drugs plot:", e$message, "\n")
  })
  
  # 3. Feature distribution plot
  tryCatch({
    features_plot <- node_properties %>%
      select(degree, betweenness, closeness) %>%
      pivot_longer(everything(), names_to = "metric", values_to = "value")
    
    p2 <- ggplot(features_plot, aes(x = value, fill = metric)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~metric, scales = "free", ncol = 3) +
      theme_minimal() +
      labs(title = "Distribution of Network Metrics") +
      theme(legend.position = "none")
    
    ggsave(file.path(output_dir, "feature_distribution.png"), p2, width = 12, height = 4)
    cat("Feature distribution plot saved\n")
  }, error = function(e) {
    cat("Could not generate feature distribution plot:", e$message, "\n")
  })
}

# ==================== MAIN PIPELINE ====================

main <- function() {
  cat("Starting drug repurposing pipeline...\n")
  
  tryCatch({
    # Step 1: Load data
    ppi_network <- load_ppi_network(ppi_file)
    drug_interactions <- load_drug_interactions(drug_files)
    
    # Step 2: Analyze network
    node_properties <- calculate_network_properties(ppi_network)
    
    # Save node properties
    write_csv(node_properties, file.path(output_dir, "node_properties.csv"))
    cat("Node properties saved to:", file.path(output_dir, "node_properties.csv"), "\n")
    
    # Step 3: Create features
    features <- create_training_features(node_properties, drug_interactions)
    
    # Save features
    write_csv(features, file.path(output_dir, "drug_features.csv"))
    cat("Drug features saved to:", file.path(output_dir, "drug_features.csv"), "\n")
    
    # Step 4: Rank drugs (using simplified unsupervised approach)
    drug_ranking <- rank_drugs_simple(features)
    
    # Save ranking
    write_csv(drug_ranking, file.path(output_dir, "drug_ranking.csv"))
    cat("Drug ranking saved to:", file.path(output_dir, "drug_ranking.csv"), "\n")
    
    # Step 5: Generate visualizations
    generate_visualizations(ppi_network, node_properties, drug_ranking)
    
    # Step 6: Generate summary report
    generate_summary_report(drug_ranking, node_properties)
    
    cat("\n=== PIPELINE COMPLETED SUCCESSFULLY ===\n")
    cat("Results saved in:", output_dir, "\n")
    
    # Print top 10 candidates
    cat("\n=== TOP 10 DRUG REPURPOSING CANDIDATES ===\n")
    print(drug_ranking %>% 
          select(rank, drug_name, repurposing_score, confidence, 
                 avg_degree, avg_betweenness) %>%
          head(10))
    
  }, error = function(e) {
    cat("\n=== ERROR IN PIPELINE ===\n")
    cat("Error message:", e$message, "\n")
    cat("Traceback:\n")
    print(e)
  })
}

# ==================== REPORT GENERATION ====================

generate_summary_report <- function(drug_ranking, node_properties) {
  cat("Generating summary report...\n")
  
  report <- c(
    "DRUG REPURPOSING ANALYSIS REPORT",
    paste("Generated:", Sys.time()),
    paste("=", paste(rep("=", 48), collapse = ""), sep = ""),
    "",
    paste("Total drugs analyzed:", nrow(drug_ranking)),
    paste("High confidence candidates:", sum(drug_ranking$confidence == "high")),
    paste("Medium confidence candidates:", sum(drug_ranking$confidence == "medium")),
    paste("Low confidence candidates:", sum(drug_ranking$confidence == "low")),
    "",
    "NETWORK STATISTICS:",
    paste("  Total genes:", nrow(node_properties)),
    paste("  Hub genes:", sum(node_properties$is_hub)),
    paste("  Bottleneck genes:", sum(node_properties$is_bottleneck)),
    paste("  Average degree:", round(mean(node_properties$degree, na.rm = TRUE), 2)),
    paste("  Maximum degree:", max(node_properties$degree, na.rm = TRUE)),
    "",
    "TOP 5 CANDIDATES:"
  )
  
  if (nrow(drug_ranking) > 0) {
    top5 <- drug_ranking %>% head(5)
    for (i in 1:min(5, nrow(top5))) {
      report <- c(report,
                 paste(i, ". ", top5$drug_name[i], 
                       " (Score: ", round(top5$repurposing_score[i], 3),
                       ", Confidence: ", top5$confidence[i], ")", sep = ""))
    }
  } else {
    report <- c(report, "No candidates found")
  }
  
  report <- c(report,
              "",
              "KEY METRICS:",
              paste("  Average repurposing score:", round(mean(drug_ranking$repurposing_score), 3)),
              paste("  Median repurposing score:", round(median(drug_ranking$repurposing_score), 3)),
              paste("  Standard deviation:", round(sd(drug_ranking$repurposing_score), 3)),
              "",
              "FILES GENERATED:",
              "  1. node_properties.csv - Network metrics for each gene",
              "  2. drug_features.csv - Features for each drug",
              "  3. drug_ranking.csv - Ranked list of drug candidates",
              "  4. network_plot.png - Visualization of PPI network",
              "  5. top_drugs_plot.png - Bar chart of top candidates",
              "  6. feature_distribution.png - Distribution of network metrics"
  )
  
  # Write report to file
  writeLines(report, file.path(output_dir, "analysis_report.txt"))
  cat("Summary report saved to:", file.path(output_dir, "analysis_report.txt"), "\n")
}

# ==================== RUN PIPELINE ====================

# Run the pipeline
cat("Drug Repurposing Analysis Script\n")
cat("===============================\n")
cat("PPI File:", ppi_file, "\n")
cat("Drug Files:", paste(drug_files, collapse = ", "), "\n")
cat("Output Directory:", output_dir, "\n")
cat("\n")

main()