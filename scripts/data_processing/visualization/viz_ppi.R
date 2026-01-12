#!/usr/bin/env Rscript

# Script to visualize filtered PPI network
# Usage: Rscript visualize_ppi.R <network_file.tsv> <output_dir>

# Load required libraries
library(igraph)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(viridis)
library(dplyr)
library(data.table)
library(scales)

# Function to create directory if it doesn't exist
create_dir <- function(dir_path) {
  if(!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat(paste0("Created directory: ", dir_path, "\n"))
  }
}

# Function to load and prepare network data
load_network <- function(network_file, gene_symbols_file = NULL) {
  cat(paste0("Loading network from: ", network_file, "\n"))
  
  # Check if file exists
  if(!file.exists(network_file)) {
    stop(paste0("Network file not found: ", network_file))
  }
  
  # Read network file
  network <- fread(network_file)
  cat(paste0("Loaded network with ", nrow(network), " interactions\n"))
  
  # Standardize column names (remove special characters, convert to lowercase)
  colnames(network) <- tolower(gsub("[^a-zA-Z0-9]", "", colnames(network)))
  
  # Try to identify edge columns (gene/protein columns)
  edge_cols <- c()
  for(col in colnames(network)) {
    if(grepl("gene|protein|node|symbol|name", col) && 
       !grepl("id|string|ensembl|entrez", col)) {
      edge_cols <- c(edge_cols, col)
    }
  }
  
  # If we found at least 2 edge columns, use them
  if(length(edge_cols) >= 2) {
    gene1_col <- edge_cols[1]
    gene2_col <- edge_cols[2]
    cat(paste0("Using columns: '", gene1_col, "' and '", gene2_col, "' as edge list\n"))
  } else {
    # Otherwise use first two columns
    gene1_col <- colnames(network)[1]
    gene2_col <- colnames(network)[2]
    cat(paste0("Using first two columns: '", gene1_col, "' and '", gene2_col, "' as edge list\n"))
  }
  
  # Check for interaction score column
  score_cols <- colnames(network)[grepl("score|combined|confidence|weight", colnames(network))]
  if(length(score_cols) > 0) {
    score_col <- score_cols[1]
    cat(paste0("Using '", score_col, "' as interaction score\n"))
  } else {
    score_col <- NULL
    cat("No interaction score column found\n")
  }
  
  # Create edge list
  edges <- data.frame(
    from = network[[gene1_col]],
    to = network[[gene2_col]],
    stringsAsFactors = FALSE
  )
  
  # Add score if available
  if(!is.null(score_col)) {
    edges$score <- as.numeric(network[[score_col]])
  }
  
  # Remove any self-loops
  edges <- edges[edges$from != edges$to, ]
  
  # Remove any duplicates
  edges <- unique(edges)
  
  cat(paste0("Edge list contains ", nrow(edges), " unique interactions\n"))
  
  # Load gene symbols if provided
  gene_symbols <- NULL
  if(!is.null(gene_symbols_file) && file.exists(gene_symbols_file)) {
    gene_symbols <- readLines(gene_symbols_file)
    cat(paste0("Loaded ", length(gene_symbols), " gene symbols\n"))
  }
  
  return(list(
    edges = edges,
    score_col = score_col,
    gene_symbols = gene_symbols
  ))
}

# Function to create and analyze network graph
create_graph <- function(edges, min_score = NULL) {
  cat("Creating network graph...\n")
  
  # Filter by score if min_score is provided
  if(!is.null(min_score) && "score" %in% colnames(edges)) {
    edges <- edges[edges$score >= min_score, ]
    cat(paste0("Filtered to ", nrow(edges), " interactions with score >= ", min_score, "\n"))
  }
  
  # Create igraph object
  if(nrow(edges) == 0) {
    stop("No edges remaining after filtering")
  }
  
  # Create graph from edge list
  if("score" %in% colnames(edges)) {
    g <- graph_from_data_frame(edges, directed = FALSE)
    E(g)$weight <- edges$score
  } else {
    g <- graph_from_data_frame(edges, directed = FALSE)
  }
  
  # Remove isolated nodes (if any were added but have no edges)
  g <- delete_vertices(g, degree(g) == 0)
  
  # Calculate network statistics
  cat("\nNetwork Statistics:\n")
  cat(paste0("Number of nodes: ", vcount(g), "\n"))
  cat(paste0("Number of edges: ", ecount(g), "\n"))
  cat(paste0("Network density: ", round(graph.density(g), 4), "\n"))
  
  # Calculate degree distribution
  degrees <- degree(g)
  cat(paste0("Average degree: ", round(mean(degrees), 2), "\n"))
  cat(paste0("Maximum degree: ", max(degrees), "\n"))
  
  # Calculate betweenness centrality
  betweenness_vals <- betweenness(g, directed = FALSE, normalized = TRUE)
  
  # Calculate clustering coefficient
  clustering_coef <- transitivity(g, type = "local")
  cat(paste0("Average clustering coefficient: ", round(mean(clustering_coef, na.rm = TRUE), 4), "\n"))
  
  # Add node attributes
  V(g)$degree <- degrees
  V(g)$betweenness <- betweenness_vals
  V(g)$clustering <- clustering_coef
  
  # Identify hub genes (top 10% by degree)
  degree_threshold <- quantile(degrees, 0.9)
  hub_genes <- names(degrees[degrees >= degree_threshold])
  V(g)$is_hub <- V(g)$name %in% hub_genes
  
  cat(paste0("Hub genes (top 10% by degree): ", length(hub_genes), "\n"))
  if(length(hub_genes) > 0) {
    cat("Top hub genes:\n")
    top_hubs <- head(sort(degrees, decreasing = TRUE), min(10, length(hub_genes)))
    for(i in seq_along(top_hubs)) {
      cat(paste0("  ", names(top_hubs)[i], ": degree = ", top_hubs[i], "\n"))
    }
  }
  
  # Detect communities (if network is large enough)
  if(vcount(g) > 10) {
    cat("\nCommunity Detection:\n")
    
    # Try different community detection algorithms
    communities_louvain <- cluster_louvain(g)
    V(g)$community_louvain <- communities_louvain$membership
    
    cat(paste0("Louvain communities: ", length(unique(communities_louvain$membership)), "\n"))
    
    # Calculate modularity
    modularity_louvain <- modularity(communities_louvain)
    cat(paste0("Louvain modularity: ", round(modularity_louvain, 4), "\n"))
  }
  
  return(g)
}

# Function to visualize network with different layouts
visualize_network <- function(g, output_dir, plot_title = "PPI Network") {
  cat(paste0("\nCreating visualizations in: ", output_dir, "\n"))
  
  # Create output directory
  create_dir(output_dir)
  
  # Get basic network info for titles
  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  
  # 1. Basic network visualization
  cat("Creating basic network visualization...\n")
  
  # Try different layouts
  layouts <- list(
    "Fruchterman-Reingold" = layout_with_fr(g),
    "Kamada-Kawai" = layout_with_kk(g),
    "Large Graph Layout" = layout_with_lgl(g),
    "Circle" = layout_in_circle(g)
  )
  
  for(layout_name in names(layouts)) {
    cat(paste0("  Creating plot with ", layout_name, " layout...\n"))
    
    png_file <- file.path(output_dir, paste0("network_", gsub("[^a-zA-Z0-9]", "_", layout_name), ".png"))
    png(png_file, width = 12, height = 10, units = "in", res = 300)
    
    # Set up plot parameters based on network size
    if(n_nodes > 100) {
      vertex_size <- 3
      vertex_label_size <- 0.8
      show_labels <- FALSE  # Don't show labels for large networks
    } else if(n_nodes > 50) {
      vertex_size <- 4
      vertex_label_size <- 1
      show_labels <- TRUE
    } else {
      vertex_size <- 6
      vertex_label_size <- 1.2
      show_labels <- TRUE
    }
    
    # Color nodes by degree
    degree_colors <- colorRampPalette(c("lightblue", "darkblue"))(100)
    node_colors <- degree_colors[cut(V(g)$degree, breaks = 100)]
    
    # Create plot
    par(mar = c(0, 0, 2, 0))
    plot(g,
         layout = layouts[[layout_name]],
         vertex.size = vertex_size,
         vertex.color = node_colors,
         vertex.frame.color = "white",
         vertex.label = if(show_labels) V(g)$name else NA,
         vertex.label.color = "black",
         vertex.label.cex = vertex_label_size,
         vertex.label.dist = 0.5,
         edge.color = adjustcolor("gray", alpha.f = 0.5),
         edge.width = if("weight" %in% edge_attr_names(g)) E(g)$weight * 2 else 1,
         main = paste0(plot_title, "\n", layout_name, " Layout\n",
                       n_nodes, " nodes, ", n_edges, " edges"))
    
    # Add legend for degree
    legend("bottomleft", 
           legend = c("Low degree", "High degree"),
           fill = c("lightblue", "darkblue"),
           bty = "n",
           cex = 0.8)
    
    dev.off()
    cat(paste0("    Saved to: ", png_file, "\n"))
  }
  
  # 2. Enhanced visualization with ggraph (if network is not too large)
  if(n_nodes <= 200) {
    cat("Creating enhanced visualization with ggraph...\n")
    
    # Convert to tidygraph for ggraph
    tg <- as_tbl_graph(g)
    
    # Create ggraph plot
    p <- ggraph(tg, layout = 'fr') + 
      geom_edge_link(aes(alpha = if("weight" %in% edge_attr_names(g)) weight else 0.5), 
                     color = "gray") +
      geom_node_point(aes(size = degree, color = degree)) +
      geom_node_text(aes(label = ifelse(degree > quantile(degree, 0.8), name, "")), 
                     repel = TRUE, size = 3) +
      scale_color_viridis(name = "Degree", option = "plasma") +
      scale_size_continuous(range = c(2, 8), name = "Degree") +
      theme_void() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      labs(title = paste0(plot_title, "\n", n_nodes, " nodes, ", n_edges, " edges"))
    
    ggsave(file.path(output_dir, "network_enhanced.png"), p, 
           width = 12, height = 10, dpi = 300)
    cat(paste0("Enhanced visualization saved to: ", file.path(output_dir, "network_enhanced.png"), "\n"))
  }
  
  # 3. Community visualization (if communities were detected)
  if("community_louvain" %in% vertex_attr_names(g)) {
    cat("Creating community visualization...\n")
    
    png_file <- file.path(output_dir, "network_communities.png")
    png(png_file, width = 12, height = 10, units = "in", res = 300)
    
    # Color by community
    n_communities <- length(unique(V(g)$community_louvain))
    community_colors <- rainbow(n_communities, alpha = 0.7)
    
    par(mar = c(0, 0, 2, 0))
    plot(g,
         layout = layout_with_fr(g),
         vertex.size = 5,
         vertex.color = community_colors[V(g)$community_louvain],
         vertex.frame.color = "white",
         vertex.label = if(n_nodes <= 100) V(g)$name else NA,
         vertex.label.color = "black",
         vertex.label.cex = 0.8,
         edge.color = adjustcolor("gray", alpha.f = 0.3),
         edge.width = 1,
         main = paste0(plot_title, "\nCommunities (Louvain algorithm)\n",
                       n_communities, " communities detected"))
    
    dev.off()
    cat(paste0("Community visualization saved to: ", png_file, "\n"))
  }
  
  # 4. Degree distribution plot
  cat("Creating degree distribution plot...\n")
  
  degree_df <- data.frame(
    degree = V(g)$degree,
    gene = V(g)$name
  )
  
  p_degree <- ggplot(degree_df, aes(x = degree)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_density(aes(y = after_stat(count)), color = "red", size = 1) +
    theme_minimal() +
    labs(
      title = "Degree Distribution",
      subtitle = paste0("Mean degree: ", round(mean(degree_df$degree), 2)),
      x = "Degree (number of connections)",
      y = "Number of genes"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  ggsave(file.path(output_dir, "degree_distribution.png"), p_degree, 
         width = 8, height = 6, dpi = 300)
  cat(paste0("Degree distribution saved to: ", file.path(output_dir, "degree_distribution.png"), "\n"))
  
  # 5. Hub gene network (subnetwork of top hubs)
  if(length(V(g)$is_hub[V(g)$is_hub]) > 1) {
    cat("Creating hub gene subnetwork...\n")
    
    # Get hub genes and their immediate neighbors
    hub_nodes <- V(g)$name[V(g)$is_hub]
    hub_neighbors <- unique(unlist(neighborhood(g, order = 1, nodes = hub_nodes)))
    hub_subgraph <- induced_subgraph(g, hub_neighbors)
    
    # Highlight hubs in red
    V(hub_subgraph)$color <- ifelse(V(hub_subgraph)$name %in% hub_nodes, "red", "lightblue")
    V(hub_subgraph)$size <- ifelse(V(hub_subgraph)$name %in% hub_nodes, 8, 4)
    
    png_file <- file.path(output_dir, "hub_network.png")
    png(png_file, width = 10, height = 8, units = "in", res = 300)
    
    par(mar = c(0, 0, 2, 0))
    plot(hub_subgraph,
         layout = layout_with_fr(hub_subgraph),
         vertex.size = V(hub_subgraph)$size,
         vertex.color = V(hub_subgraph)$color,
         vertex.frame.color = "white",
         vertex.label = V(hub_subgraph)$name,
         vertex.label.color = "black",
         vertex.label.cex = 0.8,
         edge.color = adjustcolor("gray", alpha.f = 0.5),
         edge.width = 1,
         main = paste0("Hub Gene Subnetwork\n",
                       length(hub_nodes), " hub genes + their neighbors\n",
                       vcount(hub_subgraph), " total nodes"))
    
    legend("bottomleft", 
           legend = c("Hub genes", "Neighbors"),
           fill = c("red", "lightblue"),
           bty = "n",
           cex = 0.8)
    
    dev.off()
    cat(paste0("Hub subnetwork saved to: ", png_file, "\n"))
  }
  
  # 6. Save network statistics
  cat("Saving network statistics...\n")
  
  # Create node statistics table
  node_stats <- data.frame(
    Gene = V(g)$name,
    Degree = V(g)$degree,
    Betweenness = round(V(g)$betweenness, 4),
    Clustering = round(V(g)$clustering, 4),
    Is_Hub = V(g)$is_hub,
    stringsAsFactors = FALSE
  )
  
  if("community_louvain" %in% vertex_attr_names(g)) {
    node_stats$Community <- V(g)$community_louvain
  }
  
  # Sort by degree (descending)
  node_stats <- node_stats[order(-node_stats$Degree), ]
  
  # Save node statistics
  stats_file <- file.path(output_dir, "node_statistics.tsv")
  write.table(node_stats, stats_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste0("Node statistics saved to: ", stats_file, "\n"))
  
  # Save edge list
  edge_file <- file.path(output_dir, "edge_list.tsv")
  write.table(get.data.frame(g, what = "edges"), edge_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste0("Edge list saved to: ", edge_file, "\n"))
  
  # Create summary statistics
  summary_stats <- data.frame(
    Metric = c("Number of nodes", "Number of edges", "Network density", 
               "Average degree", "Maximum degree", "Average clustering coefficient"),
    Value = c(
      vcount(g),
      ecount(g),
      round(graph.density(g), 4),
      round(mean(V(g)$degree), 2),
      max(V(g)$degree),
      round(mean(V(g)$clustering, na.rm = TRUE), 4)
    )
  )
  
  if("community_louvain" %in% vertex_attr_names(g)) {
    summary_stats <- rbind(summary_stats,
                           data.frame(
                             Metric = c("Number of communities", "Modularity"),
                             Value = c(
                               length(unique(V(g)$community_louvain)),
                               round(modularity(g, V(g)$community_louvain), 4)
                             )
                           ))
  }
  
  summary_file <- file.path(output_dir, "network_summary.tsv")
  write.table(summary_stats, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste0("Network summary saved to: ", summary_file, "\n"))
  
  cat("\nAll visualizations and statistics saved successfully!\n")
}

# Main function
main <- function() {
  cat("PPI Network Visualization Script\n")
  cat("================================\n")
  
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if(length(args) < 1) {
    cat("Usage: Rscript visualize_ppi.R <network_file.tsv> [output_dir] [gene_symbols_file]\n")
    cat("\nArguments:\n")
    cat("  network_file.tsv    : Path to your filtered PPI network file\n")
    cat("  output_dir          : (Optional) Output directory for visualizations\n")
    cat("  gene_symbols_file   : (Optional) File with list of gene symbols\n")
    cat("\nExample:\n")
    cat('  Rscript visualize_ppi.R "./ppi_results/filtered_string_network.tsv" "./network_viz"\n')
    return()
  }
  
  network_file <- args[1]
  output_dir <- if(length(args) >= 2) args[2] else "./network_visualizations"
  gene_symbols_file <- if(length(args) >= 3) args[3] else NULL
  
  # Clean paths
  network_file <- normalizePath(network_file, mustWork = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  if(!is.null(gene_symbols_file)) {
    gene_symbols_file <- normalizePath(gene_symbols_file, mustWork = FALSE)
  }
  
  # Create output directory
  create_dir(output_dir)
  
  # Load network data
  network_data <- load_network(network_file, gene_symbols_file)
  
  # Create graph
  tryCatch({
    g <- create_graph(network_data$edges)
    
    # Create plot title from filename
    plot_title <- gsub("_", " ", tools::file_path_sans_ext(basename(network_file)))
    
    # Visualize network
    visualize_network(g, output_dir, plot_title)
    
    cat("\n================================\n")
    cat("Visualization complete!\n")
    cat(paste0("Output directory: ", output_dir, "\n"))
    cat("================================\n")
    
  }, error = function(e) {
    cat(paste0("Error: ", e$message, "\n"))
  })
}

# Run the script
if(!interactive()) {
  main()
} else {
  cat("Running in interactive mode...\n")
  
  # Interactive example - modify these paths as needed
  network_file <- "./datasets/ppi/string_interactions.tsv"
  output_dir <- "./network_visualizations_original"
  gene_symbols_file <- "./ppi_results/all_gene_symbols.txt"
  
  # Create output directory
  create_dir(output_dir)
  
  # Load network data
  network_data <- load_network(network_file, gene_symbols_file)
  
  # Create graph
  g <- create_graph(network_data$edges)
  
  # Create plot title from filename
  plot_title <- gsub("_", " ", tools::file_path_sans_ext(basename(network_file)))
  
  # Visualize network
  visualize_network(g, output_dir, plot_title)
  
  cat("\n================================\n")
  cat("Visualization complete!\n")
  cat(paste0("Output directory: ", output_dir, "\n"))
  cat("================================\n")
}