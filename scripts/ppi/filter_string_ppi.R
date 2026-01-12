#!/usr/bin/env Rscript

# SIMPLE PPI FILTERING SCRIPT WITH NETWORK SUMMARY
# Usage: Rscript filter_ppi.R

library(data.table)
library(igraph)

# ========== USER INPUT ==========
cat("SIMPLE PPI FILTERING TOOL\n")
cat("===========================\n\n")

# File paths (you can change these if needed)
de_file <- "C:/Users/noahm/PycharmProjects/MarbleProject/de_results/combined/combined_DE_results_annotated_cleaned.csv"
ppi_file <- "C:/Users/noahm/PycharmProjects/MarbleProject/datasets/ppi/string_interactions.tsv"

# Ask for thresholds
cat("Enter filtering thresholds:\n")
adj_p <- as.numeric(readline("adj.P.Val threshold (e.g., 0.05): "))
logfc <- as.numeric(readline("Minimum |logFC| (e.g., 0.5): "))

use_datasets <- readline("Filter by minimum datasets? (y/n): ")
if(tolower(use_datasets) == "y") {
  min_datasets <- as.integer(readline("Minimum datasets (1-5): "))
} else {
  min_datasets <- NULL
}

output_name <- readline("Output directory name: ")
output_dir <- paste0("./", output_name)

# ========== LOAD DATA ==========
cat("\nLoading data...\n")
de_data <- fread(de_file)
ppi_data <- fread(ppi_file)

# ========== FILTER DE GENES ==========
cat("Filtering DE genes...\n")
filtered_de <- de_data[adj.P.Val < adj_p & abs(logFC) > logfc, ]

if(!is.null(min_datasets) && "n_datasets" %in% colnames(filtered_de)) {
  filtered_de <- filtered_de[n_datasets >= min_datasets, ]
}

sig_genes <- unique(filtered_de$SYMBOL)
sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]

cat(paste0("Found ", length(sig_genes), " significant genes\n"))

# ========== FILTER PPI NETWORK ==========
cat("Filtering PPI network...\n")

# Find gene columns in PPI data
if(all(c("node1", "node2") %in% colnames(ppi_data))) {
  gene1_col <- "node1"
  gene2_col <- "node2"
} else if(all(c("gene1", "gene2") %in% colnames(ppi_data))) {
  gene1_col <- "gene1"
  gene2_col <- "gene2"
} else {
  # Try to guess
  possible_cols <- colnames(ppi_data)[sapply(ppi_data, function(x) 
    all(grepl("^[A-Z][A-Z0-9-]+$", head(x, 20))))]
  if(length(possible_cols) >= 2) {
    gene1_col <- possible_cols[1]
    gene2_col <- possible_cols[2]
  } else {
    stop("Cannot find gene columns in PPI file")
  }
}

# Filter network
filtered_ppi <- ppi_data[get(gene1_col) %in% sig_genes & get(gene2_col) %in% sig_genes, ]

# ========== CREATE NETWORK SUMMARY ==========
create_network_summary <- function(filtered_ppi, gene1_col, gene2_col, sig_genes, output_dir) {
  cat("\nCreating network summary...\n")
  
  # Create igraph object
  edges <- data.frame(
    from = filtered_ppi[[gene1_col]],
    to = filtered_ppi[[gene2_col]],
    stringsAsFactors = FALSE
  )
  
  # Get all genes that should be in the network
  all_network_genes <- unique(c(filtered_ppi[[gene1_col]], filtered_ppi[[gene2_col]]))
  
  # Create graph (include isolated nodes)
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = all_network_genes))
  
  # ========== BASIC NETWORK STATISTICS ==========
  basic_stats <- list(
    total_nodes = vcount(g),
    total_edges = ecount(g),
    network_density = edge_density(g),
    average_degree = mean(degree(g)),
    max_degree = max(degree(g)),
    diameter = diameter(g),
    average_path_length = mean_distance(g),
    clustering_coefficient = transitivity(g, type = "average")
  )
  
  # ========== NODE-LEVEL STATISTICS ==========
  # Calculate centrality measures
  degree_vals <- degree(g)
  betweenness_vals <- betweenness(g)
  closeness_vals <- closeness(g)
  
  # Create node statistics table
  node_stats <- data.frame(
    gene = V(g)$name,
    degree = degree_vals,
    betweenness = betweenness_vals,
    closeness = closeness_vals,
    stringsAsFactors = FALSE
  )
  
  # Sort by degree (most connected first)
  node_stats <- node_stats[order(-node_stats$degree), ]
  
  # ========== TOP HUBS ==========
  top_hubs <- head(node_stats, 10)
  
  # ========== CONNECTED COMPONENTS ==========
  components <- components(g)
  component_summary <- data.frame(
    component_id = 1:components$no,
    size = components$csize,
    stringsAsFactors = FALSE
  )
  component_summary <- component_summary[order(-component_summary$size), ]
  
  # ========== ALS GENE ANALYSIS ==========
  als_genes <- c("SOD1", "TARDBP", "C9ORF72", "FUS", "VCP", "OPTN", "UBQLN2")
  als_in_network <- intersect(als_genes, V(g)$name)
  
  als_stats <- NULL
  if(length(als_in_network) > 0) {
    als_stats <- node_stats[node_stats$gene %in% als_in_network, ]
  }
  
  # ========== DEGREE DISTRIBUTION ==========
  degree_dist <- table(degree_vals)
  degree_dist_df <- data.frame(
    degree = as.numeric(names(degree_dist)),
    frequency = as.numeric(degree_dist),
    stringsAsFactors = FALSE
  )
  
  # ========== SAVE SUMMARY FILES ==========
  # 1. Basic network statistics
  basic_stats_df <- data.frame(
    metric = names(basic_stats),
    value = unlist(basic_stats),
    stringsAsFactors = FALSE
  )
  fwrite(basic_stats_df, file.path(output_dir, "network_basic_stats.tsv"), sep = "\t")
  
  # 2. Node statistics
  fwrite(node_stats, file.path(output_dir, "node_statistics.tsv"), sep = "\t")
  
  # 3. Top hubs
  fwrite(top_hubs, file.path(output_dir, "top_hubs.tsv"), sep = "\t")
  
  # 4. Connected components
  fwrite(component_summary, file.path(output_dir, "connected_components.tsv"), sep = "\t")
  
  # 5. Degree distribution
  fwrite(degree_dist_df, file.path(output_dir, "degree_distribution.tsv"), sep = "\t")
  
  # 6. ALS gene statistics (if any)
  if(!is.null(als_stats) && nrow(als_stats) > 0) {
    fwrite(als_stats, file.path(output_dir, "als_gene_statistics.tsv"), sep = "\t")
  }
  
  # ========== PRINT SUMMARY TO CONSOLE ==========
  cat("\n" , strrep("=", 60), "\n", sep = "")
  cat("NETWORK SUMMARY\n")
  cat(strrep("=", 60), "\n")
  
  cat("\nBASIC STATISTICS:\n")
  cat(paste0("  Total nodes: ", basic_stats$total_nodes, "\n"))
  cat(paste0("  Total edges: ", basic_stats$total_edges, "\n"))
  cat(paste0("  Network density: ", round(basic_stats$network_density, 4), "\n"))
  cat(paste0("  Average degree: ", round(basic_stats$average_degree, 2), "\n"))
  cat(paste0("  Maximum degree: ", basic_stats$max_degree, "\n"))
  cat(paste0("  Diameter: ", basic_stats$diameter, "\n"))
  cat(paste0("  Average path length: ", round(basic_stats$average_path_length, 3), "\n"))
  cat(paste0("  Clustering coefficient: ", round(basic_stats$clustering_coefficient, 4), "\n"))
  
  cat("\nCONNECTED COMPONENTS:\n")
  cat(paste0("  Number of components: ", components$no, "\n"))
  cat(paste0("  Largest component: ", max(components$csize), " genes\n"))
  if(components$no > 1) {
    cat(paste0("  Second largest: ", component_summary$size[2], " genes\n"))
  }
  
  cat("\nTOP 10 HUB GENES:\n")
  for(i in 1:nrow(top_hubs)) {
    cat(sprintf("  %2d. %-10s (degree: %3d, betweenness: %8.1f)\n", 
                i, top_hubs$gene[i], top_hubs$degree[i], top_hubs$betweenness[i]))
  }
  
  cat("\nALS GENE ANALYSIS:\n")
  cat(paste0("  ALS genes in network: ", length(als_in_network), "/7\n"))
  if(length(als_in_network) > 0) {
    cat("  Present: ", paste(als_in_network, collapse = ", "), "\n")
    cat("\n  ALS gene statistics:\n")
    for(gene in als_in_network) {
      stats <- node_stats[node_stats$gene == gene, ]
      cat(sprintf("    %-8s: degree=%2d, betweenness=%8.1f\n", 
                  gene, stats$degree, stats$betweenness))
    }
  }
  
  cat("\nDEGREE DISTRIBUTION:\n")
  cat(paste0("  Isolated nodes (degree 0): ", sum(degree_vals == 0), "\n"))
  cat(paste0("  Low connectivity (degree 1-2): ", sum(degree_vals >= 1 & degree_vals <= 2), "\n"))
  cat(paste0("  Medium connectivity (degree 3-10): ", sum(degree_vals >= 3 & degree_vals <= 10), "\n"))
  cat(paste0("  High connectivity (degree >10): ", sum(degree_vals > 10), "\n"))
  
  return(list(
    graph = g,
    basic_stats = basic_stats,
    node_stats = node_stats,
    top_hubs = top_hubs,
    components = components,
    als_in_network = als_in_network
  ))
}

# ========== SAVE RESULTS ==========
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save filtered network
output_file <- file.path(output_dir, "filtered_ppi.tsv")
fwrite(filtered_ppi, output_file, sep = "\t")
cat(paste0("Filtered network saved to: ", output_file, "\n"))

# Save gene list
genes_file <- file.path(output_dir, "significant_genes.txt")
writeLines(sort(sig_genes), genes_file)

# Save simple edge list
edges_file <- file.path(output_dir, "edge_list.tsv")
fwrite(filtered_ppi[, .(get(gene1_col), get(gene2_col))], 
       edges_file, sep = "\t", col.names = FALSE)

# ========== CREATE NETWORK SUMMARY ==========
if(nrow(filtered_ppi) > 0) {
  network_summary <- create_network_summary(filtered_ppi, gene1_col, gene2_col, sig_genes, output_dir)
} else {
  cat("\nWARNING: No interactions in filtered network. Skipping network analysis.\n")
}

# ========== FINAL SUMMARY ==========
cat("\n" , strrep("=", 60), "\n", sep = "")
cat("FILTERING COMPLETE\n")
cat(strrep("=", 60), "\n")

cat(paste0("\nTHRESHOLDS USED:\n"))
cat(paste0("  adj.P.Val < ", adj_p, "\n"))
cat(paste0("  |logFC| > ", logfc, "\n"))
if(!is.null(min_datasets)) {
  cat(paste0("  datasets >= ", min_datasets, "\n"))
}

cat(paste0("\nRESULTS SUMMARY:\n"))
cat(paste0("  Significant DE genes: ", length(sig_genes), "\n"))
if(exists("filtered_ppi") && nrow(filtered_ppi) > 0) {
  network_genes <- unique(c(filtered_ppi[[gene1_col]], filtered_ppi[[gene2_col]]))
  cat(paste0("  Genes in PPI network: ", length(network_genes), "\n"))
  cat(paste0("  Interactions in network: ", nrow(filtered_ppi), "\n"))
  
  # Check ALS genes
  als_genes <- c("SOD1", "TARDBP", "C9ORF72", "FUS", "VCP", "OPTN", "UBQLN2")
  als_in_network <- intersect(als_genes, network_genes)
  cat(paste0("  ALS genes in network: ", length(als_in_network), "/7\n"))
  if(length(als_in_network) > 0) {
    cat(paste0("    Present: ", paste(als_in_network, collapse = ", "), "\n"))
  }
}

cat(paste0("\nOUTPUT FILES CREATED:\n"))
cat(paste0("  ", output_dir, "/filtered_ppi.tsv          - Filtered PPI network\n"))
cat(paste0("  ", output_dir, "/significant_genes.txt     - List of significant genes\n"))
cat(paste0("  ", output_dir, "/edge_list.tsv            - Simple edge list\n"))
if(exists("network_summary")) {
  cat(paste0("  ", output_dir, "/network_basic_stats.tsv  - Basic network statistics\n"))
  cat(paste0("  ", output_dir, "/node_statistics.tsv      - Node-level statistics\n"))
  cat(paste0("  ", output_dir, "/top_hubs.tsv            - Top 10 hub genes\n"))
  cat(paste0("  ", output_dir, "/connected_components.tsv - Connected components\n"))
  cat(paste0("  ", output_dir, "/degree_distribution.tsv  - Degree distribution\n"))
  if(file.exists(file.path(output_dir, "als_gene_statistics.tsv"))) {
    cat(paste0("  ", output_dir, "/als_gene_statistics.tsv - ALS gene statistics\n"))
  }
}

cat(paste0("\nAll analysis complete! Check the ", output_dir, " directory for results.\n"))