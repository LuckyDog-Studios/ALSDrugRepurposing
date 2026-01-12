# Load required libraries
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(data.table)

# Function to identify annotation format - IMPROVED
identify_annotation_format <- function(gene_ids) {
  # Take a sample to identify format
  sample_ids <- head(gene_ids, 100)
  
  # Check for common formats
  format_info <- list()
  
  # Check if they look like Ensembl IDs (starting with ENS)
  ensembl_pattern <- grepl("^ENS[A-Z]*[GT]\\d+|^ENST\\d+|^ENSP\\d+", sample_ids, ignore.case = TRUE)
  format_info$ensembl <- sum(ensembl_pattern, na.rm = TRUE)
  
  # Check if they look like Entrez IDs (all numeric)
  entrez_pattern <- grepl("^\\d+$", sample_ids)
  format_info$entrez <- sum(entrez_pattern, na.rm = TRUE)
  
  # Check if they look like RefSeq IDs (starting with NM_ or NR_)
  refseq_pattern <- grepl("^[NXY]M_|^NR_|^NP_|^XP_|^YP_", sample_ids)
  format_info$refseq <- sum(refseq_pattern, na.rm = TRUE)
  
  # Check if they're already gene symbols
  symbol_pattern <- grepl("^[A-Za-z][A-Za-z0-9-]*$", sample_ids) & 
    nchar(sample_ids) > 1 & nchar(sample_ids) < 30
  format_info$symbol <- sum(symbol_pattern, na.rm = TRUE)
  
  # Check if they look like Affymetrix probes
  affy_pattern <- grepl("_at$|_st$|_x_at$|_s_at$|_a_at$", sample_ids)
  format_info$affy_probe <- sum(affy_pattern, na.rm = TRUE)
  
  # Check if they look like Illumina probes
  illumina_pattern <- grepl("^ILMN_\\d+", sample_ids, ignore.case = TRUE)
  format_info$illumina_probe <- sum(illumina_pattern, na.rm = TRUE)
  
  # Check for other probe patterns
  other_probe_pattern <- grepl("^[A-Z]{2,}\\d+_|^[A-Z]\\d+_|^[A-Z]{2,}\\d+$", sample_ids)
  format_info$other_probe <- sum(other_probe_pattern, na.rm = TRUE)
  
  # Determine the most likely format
  max_hits <- max(unlist(format_info))
  likely_format <- names(format_info)[which(unlist(format_info) == max_hits)][1]
  
  cat("Annotation format detection results:\n")
  for(fmt in names(format_info)) {
    if(format_info[[fmt]] > 0) {
      cat(paste0("  ", fmt, ": ", format_info[[fmt]], " hits\n"))
    }
  }
  cat(paste0("\nMost likely format: ", likely_format, "\n"))
  
  return(likely_format)
}

# Function to map IDs to gene symbols - IMPROVED
map_to_symbols <- function(gene_ids, format_type, original_file_name = "") {
  # Ensure gene_ids are character
  gene_ids <- as.character(gene_ids)
  
  # Remove any duplicates
  unique_ids <- unique(gene_ids)
  cat(paste0("Mapping ", length(unique_ids), " unique IDs...\n"))
  
  # Special handling for probe IDs based on file name patterns
  if(format_type %in% c("affy_probe", "illumina_probe", "other_probe") || 
     grepl("probe", format_type)) {
    return(map_probe_ids(gene_ids, format_type, original_file_name))
  }
  
  if(format_type == "entrez") {
    # Map Entrez IDs to symbols
    mapped <- tryCatch({
      select(org.Hs.eg.db,
             keys = gene_ids,
             columns = c("SYMBOL", "ENTREZID"),
             keytype = "ENTREZID")
    }, error = function(e) {
      cat("Error in mapping Entrez IDs:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(mapped)) {
      mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
      mapped <- mapped[!duplicated(mapped$ENTREZID), ]
      colnames(mapped)[colnames(mapped) == "ENTREZID"] <- "ORIGINAL_ID"
    }
    return(mapped)
    
  } else if(format_type == "ensembl") {
    # First try Ensembl gene IDs
    mapped <- tryCatch({
      # Try without version numbers
      clean_ids <- sub("\\..*$", "", gene_ids)
      select(org.Hs.eg.db,
             keys = clean_ids,
             columns = c("SYMBOL", "ENSEMBL"),
             keytype = "ENSEMBL")
    }, error = function(e) {
      cat("Error mapping as Ensembl gene IDs:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(mapped) && nrow(mapped) > 0) {
      mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
      colnames(mapped)[colnames(mapped) == "ENSEMBL"] <- "ORIGINAL_ID"
      return(mapped)
    }
    
    # If that fails, try other possibilities
    cat("Trying alternative mappings...\n")
    
    # Try if they might be transcript IDs
    if(any(grepl("^ENST", gene_ids))) {
      cat("Detected ENST (transcript) IDs. Trying to map via transcript...\n")
      mapped <- tryCatch({
        select(org.Hs.eg.db,
               keys = sub("\\..*$", "", gene_ids),
               columns = c("SYMBOL", "ENSEMBLTRANS"),
               keytype = "ENSEMBLTRANS")
      }, error = function(e) NULL)
      
      if(!is.null(mapped) && nrow(mapped) > 0) {
        colnames(mapped)[colnames(mapped) == "ENSEMBLTRANS"] <- "ORIGINAL_ID"
        mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
        return(mapped)
      }
    }
    
    return(NULL)
    
  } else if(format_type == "refseq") {
    mapped <- tryCatch({
      select(org.Hs.eg.db,
             keys = gene_ids,
             columns = c("SYMBOL", "REFSEQ"),
             keytype = "REFSEQ")
    }, error = function(e) {
      cat("Error in mapping RefSeq IDs:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(mapped)) {
      mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
      mapped <- mapped[!duplicated(mapped$REFSEQ), ]
      colnames(mapped)[colnames(mapped) == "REFSEQ"] <- "ORIGINAL_ID"
    }
    return(mapped)
    
  } else if(format_type == "symbol") {
    # Already symbols, just create a data frame
    mapped <- data.frame(SYMBOL = gene_ids, 
                         ORIGINAL_ID = as.character(gene_ids),  # Ensure character type
                         stringsAsFactors = FALSE)
    # Remove any rows where SYMBOL is NA or empty
    mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
    return(mapped)
    
  } else {
    cat(paste0("Unknown format type: ", format_type, "\n"))
    
    # Last resort: check if they might be symbols
    symbol_test <- gene_ids[grepl("^[A-Za-z][A-Za-z0-9-]*$", gene_ids) & 
                              nchar(gene_ids) > 1 & nchar(gene_ids) < 30]
    
    if(length(symbol_test) > length(gene_ids) * 0.3) {
      cat("Most IDs look like gene symbols. Treating as symbols.\n")
      mapped <- data.frame(SYMBOL = gene_ids, 
                           ORIGINAL_ID = as.character(gene_ids),  # Ensure character type
                           stringsAsFactors = FALSE)
      mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
      return(mapped)
    }
    
    return(NULL)
  }
}

# Special function for probe IDs
map_probe_ids <- function(gene_ids, format_type, original_file_name = "") {
  cat("Mapping probe IDs...\n")
  cat(paste0("Original file name: ", original_file_name, "\n"))
  
  # Try different mapping strategies
  
  # Strategy 1: Try to extract gene symbols from probe names
  # Some probes have embedded gene symbols
  if(format_type == "affy_probe" || grepl("_at$", gene_ids[1])) {
    cat("Detected Affymetrix-style probes.\n")
    
    # Check for specific datasets
    if(grepl("GSE833", original_file_name)) {
      cat("GSE833 detected - likely Affymetrix U133A.\n")
      
      # Try to map using hgu133a.db or hgu133plus2.db
      if(requireNamespace("hgu133a.db", quietly = TRUE)) {
        cat("Using hgu133a.db for Affymetrix U133A mapping...\n")
        library(hgu133a.db)
        mapped <- tryCatch({
          select(hgu133a.db,
                 keys = gene_ids,
                 columns = c("SYMBOL"),
                 keytype = "PROBEID")
        }, error = function(e) NULL)
        
        if(!is.null(mapped) && nrow(mapped) > 0) {
          mapped$ORIGINAL_ID <- as.character(mapped$PROBEID)  # Ensure character type
          mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
          return(mapped[, c("SYMBOL", "ORIGINAL_ID")])
        }
      }
      
      # Try hgu133plus2.db as fallback
      if(requireNamespace("hgu133plus2.db", quietly = TRUE)) {
        cat("Trying hgu133plus2.db for Affymetrix mapping...\n")
        library(hgu133plus2.db)
        mapped <- tryCatch({
          select(hgu133plus2.db,
                 keys = gene_ids,
                 columns = c("SYMBOL"),
                 keytype = "PROBEID")
        }, error = function(e) NULL)
        
        if(!is.null(mapped) && nrow(mapped) > 0) {
          mapped$ORIGINAL_ID <- as.character(mapped$PROBEID)  # Ensure character type
          mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
          return(mapped[, c("SYMBOL", "ORIGINAL_ID")])
        }
      }
      
      # For GSE833, try to extract numeric part as Entrez IDs
      cat("Trying to extract Entrez IDs from probe names...\n")
      numeric_ids <- gsub("[^0-9]", "", gene_ids)
      numeric_ids <- numeric_ids[nchar(numeric_ids) > 4]
      
      if(length(numeric_ids) > length(gene_ids) * 0.3) {
        cat(paste0("Found ", length(numeric_ids), " numeric IDs. Trying Entrez mapping...\n"))
        mapped <- map_to_symbols(numeric_ids, "entrez")
        if(!is.null(mapped)) {
          # Match back to original IDs
          result <- data.frame(
            SYMBOL = mapped$SYMBOL,
            ORIGINAL_ID = as.character(gene_ids[match(mapped$ORIGINAL_ID, numeric_ids)]),
            stringsAsFactors = FALSE
          )
          result <- result[!is.na(result$ORIGINAL_ID), ]
          return(result)
        }
      }
    }
  }
  
  # Strategy 2: For Illumina probes
  if(format_type == "illumina_probe" || grepl("^ILMN_", gene_ids[1])) {
    cat("Detected Illumina probes.\n")
    
    # Try to map using the illuminaHumanv4.db or similar
    if(requireNamespace("illuminaHumanv4.db", quietly = TRUE)) {
      cat("Using illuminaHumanv4.db for Illumina mapping...\n")
      library(illuminaHumanv4.db)
      mapped <- tryCatch({
        select(illuminaHumanv4.db,
               keys = gene_ids,
               columns = c("SYMBOL"),
               keytype = "PROBEID")
      }, error = function(e) NULL)
      
      if(!is.null(mapped) && nrow(mapped) > 0) {
        mapped$ORIGINAL_ID <- as.character(mapped$PROBEID)  # Ensure character type
        mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
        return(mapped[, c("SYMBOL", "ORIGINAL_ID")])
      }
    }
  }
  
  # Strategy 3: Try all keytypes in org.Hs.eg.db
  cat("Trying all keytypes in org.Hs.eg.db...\n")
  keytypes <- keytypes(org.Hs.eg.db)
  
  for(ktype in keytypes) {
    # Skip some that are unlikely
    if(ktype %in% c("GO", "GOALL", "PATH", "PMID", "ONTOLOGY", "ONTOLOGYALL")) next
    
    mapped <- tryCatch({
      # Test with first 20 IDs
      test_result <- select(org.Hs.eg.db,
                           keys = head(gene_ids, 20),
                           columns = c("SYMBOL"),
                           keytype = ktype)
      if(nrow(test_result) > 0 && sum(!is.na(test_result$SYMBOL)) > 5) {
        # If good match, map all
        select(org.Hs.eg.db,
               keys = gene_ids,
               columns = c("SYMBOL", ktype),
               keytype = ktype)
      } else {
        NULL
      }
    }, error = function(e) NULL)
    
    if(!is.null(mapped) && nrow(mapped) > 0) {
      mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
      cat(paste0("Successfully mapped using keytype: ", ktype, "\n"))
      colnames(mapped)[colnames(mapped) == ktype] <- "ORIGINAL_ID"
      mapped$ORIGINAL_ID <- as.character(mapped$ORIGINAL_ID)  # Ensure character type
      return(mapped[, c("SYMBOL", "ORIGINAL_ID")])
    }
  }
  
  # Strategy 4: If all else fails, return the original IDs as symbols
  cat("Warning: Could not map probe IDs, returning original IDs as symbols\n")
  mapped <- data.frame(
    SYMBOL = gene_ids,
    ORIGINAL_ID = as.character(gene_ids),  # Ensure character type
    stringsAsFactors = FALSE
  )
  return(mapped)
}

# Function to load and clean a DE results file
process_de_file <- function(file_path) {
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat(paste0("Processing: ", basename(file_path), "\n"))
  cat(strrep("=", 60), "\n")
  
  # Extract dataset name from file path
  dataset_name <- basename(file_path)
  dataset_name <- gsub("_DE_results\\.csv$|_annotated_de_.*\\.csv$", "", dataset_name)
  
  # Read the file
  de_data <- tryCatch({
    # Try reading with data.table for better performance
    data <- fread(file_path, data.table = FALSE)
    cat(paste0("Successfully read ", nrow(data), " rows\n"))
    data
  }, error = function(e) {
    cat(paste0("Error reading file with fread: ", e$message, "\n"))
    cat("Trying with read.csv...\n")
    read.csv(file_path, stringsAsFactors = FALSE)
  })
  
  # Standardize column names (case-insensitive)
  colnames(de_data) <- tolower(colnames(de_data))
  
  # Identify gene ID column
  gene_id_col <- NULL
  
  # Common column names for gene identifiers
  possible_gene_cols <- c("gene", "gene_id", "geneid", "id", "gene_symbol", "symbol", 
                          "probe_id", "probeid", "ensembl", "entrez", "refseq")
  
  for(col in possible_gene_cols) {
    if(col %in% colnames(de_data)) {
      gene_id_col <- col
      cat(paste0("Found gene ID column: ", col, "\n"))
      break
    }
  }
  
  # If not found by exact match, try partial matches
  if(is.null(gene_id_col)) {
    for(col in colnames(de_data)) {
      if(grepl("gene|id|symbol|probe|ensembl|entrez", col, ignore.case = TRUE)) {
        gene_id_col <- col
        cat(paste0("Found gene ID column (partial match): ", col, "\n"))
        break
      }
    }
  }
  
  if(is.null(gene_id_col)) {
    cat("Warning: Could not identify gene ID column. Using first column.\n")
    gene_id_col <- colnames(de_data)[1]
  }
  
  # Extract gene IDs and ensure they are character type
  gene_ids <- as.character(de_data[[gene_id_col]])
  
  # Remove any NA values
  gene_ids <- gene_ids[!is.na(gene_ids)]
  
  if(length(gene_ids) == 0) {
    stop("No gene IDs found in the dataset")
  }
  
  # Identify annotation format
  format_type <- identify_annotation_format(gene_ids)
  
  # Map to gene symbols
  mapping <- map_to_symbols(gene_ids, format_type, basename(file_path))
  
  if(is.null(mapping) || nrow(mapping) == 0) {
    cat("Warning: Could not map gene IDs. Using original IDs as symbols.\n")
    mapping <- data.frame(
      SYMBOL = gene_ids,
      ORIGINAL_ID = as.character(gene_ids),  # Ensure character type
      stringsAsFactors = FALSE
    )
  }
  
  # Ensure ORIGINAL_ID is character in the original data
  de_data$ORIGINAL_ID <- as.character(de_data[[gene_id_col]])
  
  # Remove the original gene ID column to avoid duplication
  de_data <- de_data[, !colnames(de_data) %in% gene_id_col]
  
  # Merge mapping with original data
  de_data_annotated <- merge(de_data, mapping, by = "ORIGINAL_ID", all.x = TRUE)
  
  # Remove rows without symbols
  de_data_annotated <- de_data_annotated[!is.na(de_data_annotated$SYMBOL) & 
                                           de_data_annotated$SYMBOL != "", ]
  
  cat(paste0("After annotation: ", nrow(de_data_annotated), " rows with gene symbols\n"))
  
  # Add dataset identifier
  de_data_annotated$Dataset <- dataset_name
  
  # Standardize common column names for easier combination
  # Keep original column names but rename some common ones for consistency
  colnames(de_data_annotated) <- gsub("logfc|log2fc|log2foldchange", "logFC", 
                                      colnames(de_data_annotated), ignore.case = TRUE)
  colnames(de_data_annotated) <- gsub("pval|pvalue|p.val", "PValue", 
                                      colnames(de_data_annotated), ignore.case = TRUE)
  colnames(de_data_annotated) <- gsub("adjpval|adj.pval|padj|fdr", "adj.P.Val", 
                                      colnames(de_data_annotated), ignore.case = TRUE)
  
  # Ensure we have standardized column names
  if(!"logFC" %in% colnames(de_data_annotated)) {
    # Try to find log fold change column
    logfc_col <- grep("log|fc|fold", colnames(de_data_annotated), ignore.case = TRUE, value = TRUE)[1]
    if(!is.na(logfc_col)) {
      colnames(de_data_annotated)[colnames(de_data_annotated) == logfc_col] <- "logFC"
    }
  }
  
  if(!"PValue" %in% colnames(de_data_annotated)) {
    # Try to find p-value column
    pval_col <- grep("pval|p.value", colnames(de_data_annotated), ignore.case = TRUE, value = TRUE)[1]
    if(!is.na(pval_col)) {
      colnames(de_data_annotated)[colnames(de_data_annotated) == pval_col] <- "PValue"
    }
  }
  
  if(!"adj.P.Val" %in% colnames(de_data_annotated)) {
    # Try to find adjusted p-value column
    adjpval_col <- grep("adj|fdr|padj", colnames(de_data_annotated), ignore.case = TRUE, value = TRUE)[1]
    if(!is.na(adjpval_col)) {
      colnames(de_data_annotated)[colnames(de_data_annotated) == adjpval_col] <- "adj.P.Val"
    }
  }
  
  # Convert numeric columns to appropriate types
  if("logFC" %in% colnames(de_data_annotated)) {
    de_data_annotated$logFC <- as.numeric(de_data_annotated$logFC)
  }
  if("PValue" %in% colnames(de_data_annotated)) {
    de_data_annotated$PValue <- as.numeric(de_data_annotated$PValue)
  }
  if("adj.P.Val" %in% colnames(de_data_annotated)) {
    de_data_annotated$adj.P.Val <- as.numeric(de_data_annotated$adj.P.Val)
  }
  
  # Ensure ORIGINAL_ID is character
  de_data_annotated$ORIGINAL_ID <- as.character(de_data_annotated$ORIGINAL_ID)
  
  # Keep only essential columns plus original ones
  essential_cols <- c("SYMBOL", "logFC", "PValue", "adj.P.Val", "Dataset", "ORIGINAL_ID")
  existing_cols <- essential_cols[essential_cols %in% colnames(de_data_annotated)]
  
  # Add any additional columns (keeping original column names)
  other_cols <- setdiff(colnames(de_data_annotated), existing_cols)
  final_cols <- c(existing_cols, other_cols)
  
  de_data_annotated <- de_data_annotated[, final_cols, drop = FALSE]
  
  return(de_data_annotated)
}

# Main processing function
process_all_datasets <- function(de_files) {
  # List to store processed datasets
  all_datasets <- list()
  
  # Process each file
  for(i in seq_along(de_files)) {
    file_path <- de_files[i]
    cat(paste0("\nProcessing file ", i, "/", length(de_files), ": ", basename(file_path), "\n"))
    
    tryCatch({
      processed_data <- process_de_file(file_path)
      all_datasets[[i]] <- processed_data
      cat(paste0("Successfully processed ", basename(file_path), "\n"))
    }, error = function(e) {
      cat(paste0("Error processing ", basename(file_path), ": ", e$message, "\n"))
      all_datasets[[i]] <- NULL
    })
  }
  
  # Remove any NULL entries
  all_datasets <- all_datasets[!sapply(all_datasets, is.null)]
  
  if(length(all_datasets) == 0) {
    stop("No datasets were successfully processed")
  }
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("Combining datasets...\n")
  cat(strrep("=", 60), "\n")
  
  # Combine all datasets with consistent column types
  combined_data <- bind_rows(all_datasets)
  
  # Count how many datasets each gene appears in
  gene_counts <- combined_data %>%
    group_by(SYMBOL) %>%
    summarise(
      n_datasets = n_distinct(Dataset),
      datasets = paste(unique(Dataset), collapse = ";"),
      .groups = "drop"
    )
  
  # Merge the counts back to the combined data
  combined_data <- combined_data %>%
    left_join(gene_counts, by = "SYMBOL")
  
  # Reorder columns to put important ones first
  col_order <- c("SYMBOL", "n_datasets", "datasets", "Dataset", "logFC", 
                 "PValue", "adj.P.Val", "ORIGINAL_ID")
  existing_cols <- col_order[col_order %in% colnames(combined_data)]
  other_cols <- setdiff(colnames(combined_data), existing_cols)
  combined_data <- combined_data[, c(existing_cols, other_cols)]
  
  cat(paste0("Final combined dataset has ", nrow(combined_data), " rows\n"))
  cat(paste0("Unique genes: ", length(unique(combined_data$SYMBOL)), "\n"))
  
  # Summary statistics
  cat("\nDataset summary:\n")
  dataset_summary <- combined_data %>%
    group_by(Dataset) %>%
    summarise(
      n_rows = n(),
      n_genes = n_distinct(SYMBOL),
      .groups = "drop"
    )
  
  for(i in 1:nrow(dataset_summary)) {
    cat(paste0("  ", dataset_summary$Dataset[i], ": ", 
               dataset_summary$n_rows[i], " rows, ", 
               dataset_summary$n_genes[i], " unique genes\n"))
  }
  
  cat("\nGene occurrence across datasets:\n")
  occurrence_table <- table(combined_data$n_datasets)
  for(n in sort(unique(combined_data$n_datasets))) {
    n_genes <- length(unique(combined_data$SYMBOL[combined_data$n_datasets == n]))
    cat(paste0("  Genes in ", n, " dataset(s): ", n_genes, " unique genes\n"))
  }
  
  return(combined_data)
}

# Your list of DE result files
de_files <- c(
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE833\\GSE833_DE_results.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE76220\\GSE76220_DE_results.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE112676\\GSE112676_DE_results.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE118336\\GSE118336_annotated_de_FUS_H517D_Mutant_vs_Control.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE118336\\GSE118336_annotated_de_FUS_Heterozygous_vs_Control.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE124439\\GSE124439_DE_results.csv"
)

# Process all datasets
combined_de_results <- process_all_datasets(de_files)

# Save the combined results
output_file <- "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\combined\\combined_DE_results_annotated.csv"
write.csv(combined_de_results, output_file, row.names = FALSE)
cat(paste0("\nCombined results saved to: ", output_file, "\n"))

# Optional: Save a summary file with only key information
summary_data <- combined_de_results %>%
  select(SYMBOL, n_datasets, datasets, Dataset, logFC, PValue, adj.P.Val) %>%
  arrange(desc(n_datasets), SYMBOL)

summary_file <- "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\combined\\combined_DE_results_summary.csv"
write.csv(summary_data, summary_file, row.names = FALSE)
cat(paste0("Summary file saved to: ", summary_file, "\n"))

# Optional: Save genes by frequency
genes_by_frequency <- combined_de_results %>%
  distinct(SYMBOL, n_datasets, datasets) %>%
  arrange(desc(n_datasets))

frequency_file <- "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\combined\\genes_by_dataset_frequency.csv"
write.csv(genes_by_frequency, frequency_file, row.names = FALSE)
cat(paste0("Genes by frequency saved to: ", frequency_file, "\n"))

cat("\nProcessing complete!\n")