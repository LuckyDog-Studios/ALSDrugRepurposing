# Enhanced visualization function with dataset-specific thresholds
create_comprehensive_visualizations <- function(significant_lists, consensus_results, all_data = NULL, output_dir = "results", 
                                               dataset_thresholds = NULL) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Default thresholds if not provided
  if (is.null(dataset_thresholds)) {
    dataset_thresholds <- list()
    for (ds_name in names(significant_lists)) {
      dataset_thresholds[[ds_name]] <- list(
        adj_pval_threshold = 0.05,
        logFC_threshold = 0.5
      )
    }
  }
  
  # Helper function to get thresholds for a dataset
  get_thresholds <- function(dataset_name) {
    if (dataset_name %in% names(dataset_thresholds)) {
      return(dataset_thresholds[[dataset_name]])
    } else {
      # Return defaults if not specified
      return(list(
        adj_pval_threshold = 0.05,
        logFC_threshold = 0.5
      ))
    }
  }
  
  # 1. VOLCANO PLOTS for each dataset
  create_volcano_plots <- function(all_data, output_dir) {
    if (!is.null(all_data)) {
      cat("Creating volcano plots...\n")
      for (ds_name in names(all_data)) {
        tryCatch({
          data <- all_data[[ds_name]]
          if (nrow(data) > 0) {
            # Get dataset-specific thresholds
            thresholds <- get_thresholds(ds_name)
            adj_pval_thresh <- thresholds$adj_pval_threshold
            logFC_thresh <- thresholds$logFC_threshold
            
            # Create volcano plot data
            volcano_data <- data %>%
              mutate(
                significance = case_when(
                  adj.P.Val < adj_pval_thresh & abs(logFC) > logFC_thresh ~ "Significant",
                  TRUE ~ "Not significant"
                ),
                log10_pval = -log10(P.Value)
              )
            
            p <- ggplot(volcano_data, aes(x = logFC, y = log10_pval, color = significance)) +
              geom_point(alpha = 0.6, size = 1) +
              scale_color_manual(values = c("Not significant" = "gray", "Significant" = "red")) +
              labs(
                title = paste("Volcano Plot -", ds_name),
                subtitle = paste("Thresholds: adj.P.Val <", adj_pval_thresh, ", |logFC| >", logFC_thresh),
                x = "log2 Fold Change", 
                y = "-log10 P-value"
              ) +
              theme_minimal() +
              geom_hline(yintercept = -log10(adj_pval_thresh), linetype = "dashed", alpha = 0.5) +
              geom_vline(xintercept = c(-logFC_thresh, logFC_thresh), linetype = "dashed", alpha = 0.5)
            
            ggsave(file.path(output_dir, paste0("volcano_", ds_name, ".pdf")), p, width = 8, height = 6)
            ggsave(file.path(output_dir, paste0("volcano_", ds_name, ".png")), p, width = 8, height = 6)
          }
        }, error = function(e) {
          cat("Error creating volcano plot for", ds_name, ":", e$message, "\n")
        })
      }
    }
  }
  
  # 2. MA PLOTS for each dataset
  create_ma_plots <- function(all_data, output_dir) {
    if (!is.null(all_data)) {
      cat("Creating MA plots...\n")
      for (ds_name in names(all_data)) {
        tryCatch({
          data <- all_data[[ds_name]]
          if (nrow(data) > 0) {
            # Get dataset-specific thresholds
            thresholds <- get_thresholds(ds_name)
            adj_pval_thresh <- thresholds$adj_pval_threshold
            logFC_thresh <- thresholds$logFC_threshold
            
            # For MA plot we need average expression - if not available, use rank
            if ("AveExpr" %in% colnames(data)) {
              ma_data <- data
            } else {
              # Create pseudo average expression based on rank
              ma_data <- data %>%
                mutate(AveExpr = rank(logFC) / n())
            }
            
            ma_data <- ma_data %>%
              mutate(significance = ifelse(adj.P.Val < adj_pval_thresh & abs(logFC) > logFC_thresh, "Significant", "Not significant"))
            
            p <- ggplot(ma_data, aes(x = AveExpr, y = logFC, color = significance)) +
              geom_point(alpha = 0.6, size = 1) +
              scale_color_manual(values = c("Not significant" = "gray", "Significant" = "red")) +
              labs(
                title = paste("MA Plot -", ds_name),
                subtitle = paste("Thresholds: adj.P.Val <", adj_pval_thresh, ", |logFC| >", logFC_thresh),
                x = "Average Expression", 
                y = "log2 Fold Change"
              ) +
              theme_minimal() +
              geom_hline(yintercept = 0, color = "black") +
              geom_hline(yintercept = c(-logFC_thresh, logFC_thresh), linetype = "dashed", alpha = 0.5) +
              geom_smooth(method = "loess", color = "blue", se = FALSE)
            
            ggsave(file.path(output_dir, paste0("ma_plot_", ds_name, ".pdf")), p, width = 8, height = 6)
            ggsave(file.path(output_dir, paste0("ma_plot_", ds_name, ".png")), p, width = 8, height = 6)
          }
        }, error = function(e) {
          cat("Error creating MA plot for", ds_name, ":", e$message, "\n")
        })
      }
    }
  }
  
  # 3. P-VALUE DISTRIBUTION plots
  create_pvalue_distributions <- function(all_data, output_dir) {
    if (!is.null(all_data)) {
      cat("Creating p-value distribution plots...\n")
      
      # Individual p-value distributions
      for (ds_name in names(all_data)) {
        tryCatch({
          data <- all_data[[ds_name]]
          if (nrow(data) > 0) {
            # Get dataset-specific thresholds
            thresholds <- get_thresholds(ds_name)
            adj_pval_thresh <- thresholds$adj_pval_threshold
            
            p <- ggplot(data, aes(x = P.Value)) +
              geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
              labs(
                title = paste("P-value Distribution -", ds_name),
                subtitle = paste("adj.P.Val threshold:", adj_pval_thresh),
                x = "P-value", y = "Frequency"
              ) +
              theme_minimal() +
              geom_vline(xintercept = adj_pval_thresh, linetype = "dashed", color = "red", linewidth = 0.5)
            
            ggsave(file.path(output_dir, paste0("pvalue_dist_", ds_name, ".pdf")), p, width = 8, height = 6)
            ggsave(file.path(output_dir, paste0("pvalue_dist_", ds_name, ".png")), p, width = 8, height = 6)
          }
        }, error = function(e) {
          cat("Error creating p-value distribution for", ds_name, ":", e$message, "\n")
        })
      }
      
      # Combined p-value distribution
      pval_data <- data.frame()
      threshold_data <- data.frame()
      
      for (ds_name in names(all_data)) {
        temp_data <- all_data[[ds_name]] %>%
          select(P.Value) %>%
          mutate(Dataset = ds_name)
        pval_data <- rbind(pval_data, temp_data)
        
        # Add threshold info
        thresholds <- get_thresholds(ds_name)
        threshold_data <- rbind(threshold_data, data.frame(
          Dataset = ds_name,
          Threshold = thresholds$adj_pval_threshold
        ))
      }
      
      if (nrow(pval_data) > 0) {
        p <- ggplot(pval_data, aes(x = P.Value, fill = Dataset)) +
          geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
          labs(
            title = "P-value Distributions Across Datasets",
            x = "P-value", y = "Frequency"
          ) +
          theme_minimal() +
          facet_wrap(~Dataset, scales = "free_y") +
          # Add threshold lines
          geom_vline(data = threshold_data, aes(xintercept = Threshold), 
                     linetype = "dashed", color = "red", linewidth = 0.5)
        
        ggsave(file.path(output_dir, "pvalue_distributions_combined.pdf"), p, width = 10, height = 8)
        ggsave(file.path(output_dir, "pvalue_distributions_combined.png"), p, width = 10, height = 8)
      }
    }
  }
  
  # 4. LOGFC DISTRIBUTION plots
  create_logfc_distributions <- function(all_data, output_dir) {
    if (!is.null(all_data)) {
      cat("Creating logFC distribution plots...\n")
      
      logfc_data <- data.frame()
      threshold_data <- data.frame()
      
      for (ds_name in names(all_data)) {
        temp_data <- all_data[[ds_name]] %>%
          select(logFC) %>%
          mutate(Dataset = ds_name)
        logfc_data <- rbind(logfc_data, temp_data)
        
        # Add threshold info
        thresholds <- get_thresholds(ds_name)
        threshold_data <- rbind(threshold_data, data.frame(
          Dataset = ds_name,
          Threshold = thresholds$logFC_threshold
        ))
      }
      
      if (nrow(logfc_data) > 0) {
        # Density plot
        p1 <- ggplot(logfc_data, aes(x = logFC, fill = Dataset)) +
          geom_density(alpha = 0.5) +
          labs(
            title = "logFC Distributions Across Datasets",
            x = "log2 Fold Change", y = "Density"
          ) +
          theme_minimal() +
          geom_vline(xintercept = 0, linetype = "dashed") +
          # Add threshold lines
          geom_vline(data = threshold_data, aes(xintercept = Threshold), 
                     linetype = "dashed", color = "red", linewidth = 0.5) +
          geom_vline(data = threshold_data, aes(xintercept = -Threshold), 
                     linetype = "dashed", color = "red", linewidth = 0.5)
        
        # Boxplot
        p2 <- ggplot(logfc_data, aes(x = Dataset, y = logFC, fill = Dataset)) +
          geom_boxplot(alpha = 0.7) +
          labs(
            title = "logFC Distributions - Boxplot",
            x = "Dataset", y = "log2 Fold Change"
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          # Add threshold annotations
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_hline(yintercept = mean(threshold_data$Threshold), 
                     linetype = "dashed", color = "red", linewidth = 0.5) +
          geom_hline(yintercept = -mean(threshold_data$Threshold), 
                     linetype = "dashed", color = "red", linewidth = 0.5)
        
        ggsave(file.path(output_dir, "logfc_density.pdf"), p1, width = 10, height = 6)
        ggsave(file.path(output_dir, "logfc_boxplot.pdf"), p2, width = 10, height = 6)
        ggsave(file.path(output_dir, "logfc_density.png"), p1, width = 10, height = 6)
        ggsave(file.path(output_dir, "logfc_boxplot.png"), p2, width = 10, height = 6)
      }
    }
  }
  
  # 5. SIGNIFICANT GENE SUMMARY plots - COMPLETELY FIXED VERSION
  create_summary_plots <- function(significant_lists, output_dir, get_thresholds_func) {
    cat("Creating summary plots...\n")
    
    # Create summary with threshold information
    summary_data <- data.frame()
    for (ds_name in names(significant_lists)) {
      thresholds <- get_thresholds_func(ds_name)
      sig_data <- significant_lists[[ds_name]]
      
      if (nrow(sig_data) > 0) {
        up_genes <- sum(sig_data$logFC > 0, na.rm = TRUE)
        down_genes <- sum(sig_data$logFC < 0, na.rm = TRUE)
        
        summary_data <- rbind(summary_data, data.frame(
          Dataset = ds_name,
          Significant_Genes = nrow(sig_data),
          Up_regulated = up_genes,
          Down_regulated = down_genes,
          adjP_Threshold = thresholds$adj_pval_threshold,
          logFC_Threshold = thresholds$logFC_threshold
        ))
      } else {
        summary_data <- rbind(summary_data, data.frame(
          Dataset = ds_name,
          Significant_Genes = 0,
          Up_regulated = 0,
          Down_regulated = 0,
          adjP_Threshold = thresholds$adj_pval_threshold,
          logFC_Threshold = thresholds$logFC_threshold
        ))
      }
    }
    
    if (nrow(summary_data) > 0) {
      # Bar plot of significant genes with thresholds in subtitle
      p1 <- ggplot(summary_data, aes(x = Dataset, y = Significant_Genes, fill = Dataset)) +
        geom_bar(stat = "identity", alpha = 0.7) +
        labs(
          title = "Number of Significant Genes per Dataset",
          subtitle = "Dataset-specific thresholds applied",
          x = "Dataset", y = "Number of Significant Genes"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_text(aes(label = Significant_Genes), vjust = -0.3)
      
      # Up/down regulated genes
      regulation_data <- summary_data %>%
        select(Dataset, Up_regulated, Down_regulated) %>%
        pivot_longer(cols = c(Up_regulated, Down_regulated), 
                     names_to = "Regulation", values_to = "Count")
      
      p2 <- ggplot(regulation_data, aes(x = Dataset, y = Count, fill = Regulation)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        labs(
          title = "Up/Down Regulated Genes per Dataset",
          subtitle = "Dataset-specific thresholds applied",
          x = "Dataset", y = "Number of Genes"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = c("Up_regulated" = "red", "Down_regulated" = "blue"),
                          labels = c("Up-regulated", "Down-regulated"))
      
      # Create threshold summary table
      threshold_table <- summary_data %>%
        select(Dataset, adjP_Threshold, logFC_Threshold) %>%
        distinct()
      
      write.csv(threshold_table, file.path(output_dir, "dataset_thresholds.csv"), row.names = FALSE)
      
      ggsave(file.path(output_dir, "regulation_summary.pdf"), p2, width = 10, height = 6)
      ggsave(file.path(output_dir, "regulation_summary.png"), p2, width = 10, height = 6)
      ggsave(file.path(output_dir, "significant_genes_summary.pdf"), p1, width = 10, height = 6)
      ggsave(file.path(output_dir, "significant_genes_summary.png"), p1, width = 10, height = 6)
    }
  }
  
  # 6. HEATMAP of consensus genes across datasets
  create_consensus_heatmap <- function(consensus_results, output_dir) {
    if (length(consensus_results$consensus_genes) > 0 && nrow(consensus_results$presence_matrix) > 0) {
      cat("Creating consensus heatmap...\n")
      
      # Convert presence matrix to data frame for plotting
      heatmap_data <- as.data.frame(consensus_results$presence_matrix)
      heatmap_data$Gene <- rownames(heatmap_data)
      
      # Melt for ggplot
      heatmap_melt <- heatmap_data %>%
        pivot_longer(cols = -Gene, names_to = "Dataset", values_to = "Presence")
      
      p <- ggplot(heatmap_melt, aes(x = Dataset, y = reorder(Gene, consensus_results$gene_counts[Gene]), 
                                   fill = factor(Presence))) +
        geom_tile(color = "white") +
        scale_fill_manual(values = c("0" = "white", "1" = "darkred"), 
                         labels = c("Absent", "Present")) +
        labs(title = "Consensus Genes Across Datasets",
             x = "Dataset", y = "Genes", fill = "Presence") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 6),
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(output_dir, "consensus_heatmap.pdf"), p, width = 10, height = max(6, length(consensus_results$consensus_genes) * 0.1))
      ggsave(file.path(output_dir, "consensus_heatmap.png"), p, width = 10, height = max(6, length(consensus_results$consensus_genes) * 0.1))
    }
  }
  
  # 7. VENN DIAGRAM for dataset overlaps (for up to 5 datasets)
  create_venn_diagrams <- function(significant_lists, output_dir) {
    non_empty_lists <- significant_lists[sapply(significant_lists, nrow) > 0]
    
    if (length(non_empty_lists) >= 2 && length(non_empty_lists) <= 5) {
      cat("Creating Venn diagrams...\n")
      
      # Extract gene lists
      gene_lists <- lapply(non_empty_lists, function(x) unique(x$GeneSymbol_clean))
      
      # Create Venn diagram
      venn_plot <- venn.diagram(
        x = gene_lists,
        category.names = names(gene_lists),
        filename = NULL,
        output = TRUE,
        imagetype = "pdf",
        height = 8, 
        width = 8,
        resolution = 300,
        compression = "lzw",
        lwd = 2,
        lty = 'blank',
        fill = c("red", "blue", "green", "yellow", "purple")[1:length(gene_lists)],
        alpha = 0.5,
        cex = 1.5,
        fontface = "bold",
        fontfamily = "sans",
        cat.cex = 1.2,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135, -135, 0)[1:length(gene_lists)],
        cat.dist = c(0.055, 0.055, 0.055, 0.055, 0.055)[1:length(gene_lists)],
        cat.fontfamily = "sans",
        rotation = 1
      )
      
      # Save Venn diagram
      pdf(file.path(output_dir, "venn_diagram.pdf"), width = 8, height = 8)
      grid.draw(venn_plot)
      dev.off()
      
      png(file.path(output_dir, "venn_diagram.png"), width = 800, height = 800)
      grid.draw(venn_plot)
      dev.off()
    } else if (length(non_empty_lists) > 5) {
      cat("Too many datasets for Venn diagram (max 5), using upset plot instead\n")
    }
  }
  
  # 8. DATASET CORRELATION heatmap
  create_correlation_heatmap <- function(significant_lists, output_dir) {
    if (length(significant_lists) >= 2) {
      cat("Creating dataset correlation heatmap...\n")
      
      # Get all unique genes
      all_genes <- unique(unlist(lapply(significant_lists, function(x) unique(x$GeneSymbol_clean))))
      
      if (length(all_genes) > 0) {
        # Create binary matrix
        binary_matrix <- matrix(0, nrow = length(all_genes), ncol = length(significant_lists))
        rownames(binary_matrix) <- all_genes
        colnames(binary_matrix) <- names(significant_lists)
        
        for (i in 1:length(significant_lists)) {
          dataset_genes <- unique(significant_lists[[i]]$GeneSymbol_clean)
          binary_matrix[rownames(binary_matrix) %in% dataset_genes, i] <- 1
        }
        
        # Calculate Jaccard similarity
        jaccard_similarity <- function(x, y) {
          intersection <- sum(x & y)
          union <- sum(x | y)
          return(intersection / union)
        }
        
        # Create similarity matrix
        n_datasets <- ncol(binary_matrix)
        similarity_matrix <- matrix(0, nrow = n_datasets, ncol = n_datasets)
        rownames(similarity_matrix) <- colnames(binary_matrix)
        colnames(similarity_matrix) <- colnames(binary_matrix)
        
        for (i in 1:n_datasets) {
          for (j in 1:n_datasets) {
            similarity_matrix[i, j] <- jaccard_similarity(binary_matrix[, i], binary_matrix[, j])
          }
        }
        
        # Plot correlation heatmap
        similarity_melt <- as.data.frame(similarity_matrix) %>%
          rownames_to_column("Dataset1") %>%
          pivot_longer(cols = -Dataset1, names_to = "Dataset2", values_to = "Similarity")
        
        p <- ggplot(similarity_melt, aes(x = Dataset1, y = Dataset2, fill = Similarity)) +
          geom_tile(color = "white") +
          scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, limits = c(0, 1)) +
          labs(title = "Dataset Similarity (Jaccard Index)",
               x = "", y = "", fill = "Similarity") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          geom_text(aes(label = round(Similarity, 2)), color = "black", size = 3)
        
        ggsave(file.path(output_dir, "dataset_similarity_heatmap.pdf"), p, width = 8, height = 6)
        ggsave(file.path(output_dir, "dataset_similarity_heatmap.png"), p, width = 8, height = 6)
      }
    }
  }
  
  # Execute all plotting functions - PASS get_thresholds to create_summary_plots
  create_volcano_plots(all_data, output_dir)
  create_ma_plots(all_data, output_dir)
  create_pvalue_distributions(all_data, output_dir)
  create_logfc_distributions(all_data, output_dir)
  create_summary_plots(significant_lists, output_dir, get_thresholds_func = get_thresholds)
  create_consensus_heatmap(consensus_results, output_dir)
  create_venn_diagrams(significant_lists, output_dir)
  create_correlation_heatmap(significant_lists, output_dir)
  
  cat("All visualizations saved to:", output_dir, "\n")
  cat("Dataset thresholds summary saved as: dataset_thresholds.csv\n")
}