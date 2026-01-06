#!/usr/bin/env Rscript

# Script to process differential expression results, identify gene annotation format,
# map to gene symbols, and filter STRING PPI network

# Load required libraries
library(dplyr)
library(data.table)
library(annotate)
library(org.Hs.eg.db)

# Function to clean file paths (remove quotes and trim whitespace)
clean_file_path <- function(path) {
  # Remove quotes
  path <- gsub('^"|"$', '', path)
  path <- gsub("^'|'$", '', path)
  # Trim whitespace
  path <- trimws(path)
  # Convert to absolute path if relative
  if(!grepl("^[A-Za-z]:|^/", path)) {
    path <- normalizePath(path, mustWork = FALSE)
  }
  return(path)
}

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
  # More permissive pattern for gene symbols
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
map_to_symbols <- function(gene_ids, format_type, file_name = "") {
  # Ensure gene_ids are character
  gene_ids <- as.character(gene_ids)
  
  # Remove any duplicates
  unique_ids <- unique(gene_ids)
  cat(paste0("Mapping ", length(unique_ids), " unique IDs...\n"))
  
  # Special handling for probe IDs based on file name patterns
  if(format_type %in% c("affy_probe", "illumina_probe", "other_probe") || 
     grepl("probe", format_type)) {
    return(map_probe_ids(gene_ids, format_type, file_name))
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
        colnames(mapped)[colnames(mapped) == "ENSEMBLTRANS"] <- "ENSEMBL"
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
    }
    return(mapped)
    
  } else if(format_type == "symbol") {
    # Already symbols, just create a data frame
    mapped <- data.frame(SYMBOL = gene_ids, 
                         ORIGINAL_ID = gene_ids,
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
                           ORIGINAL_ID = gene_ids,
                           stringsAsFactors = FALSE)
      mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
      return(mapped)
    }
    
    return(NULL)
  }
}

# Special function for probe IDs
map_probe_ids <- function(gene_ids, format_type, file_name = "") {
  cat("Mapping probe IDs...\n")
  
  # Try different mapping strategies
  
  # Strategy 1: Try to extract gene symbols from probe names
  # Some probes have embedded gene symbols
  if(format_type == "affy_probe" || grepl("_at$", gene_ids[1])) {
    cat("Detected Affymetrix-style probes.\n")
    
    # Install and try hgu133plus2.db if available
    if(requireNamespace("hgu133plus2.db", quietly = TRUE)) {
      cat("Using hgu133plus2.db for Affymetrix mapping...\n")
      library(hgu133plus2.db)
      mapped <- tryCatch({
        select(hgu133plus2.db,
               keys = gene_ids,
               columns = c("SYMBOL"),
               keytype = "PROBEID")
      }, error = function(e) NULL)
      
      if(!is.null(mapped) && nrow(mapped) > 0) {
        mapped$ORIGINAL_ID <- mapped$PROBEID
        mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
        return(mapped)
      }
    }
    
    # Try to guess from file name
    if(grepl("GSE833", file_name)) {
      cat("GSE833 is likely Affymetrix U133A. Trying alternative mapping...\n")
      # For GSE833, try to extract numeric part
      numeric_ids <- gsub("[^0-9]", "", gene_ids)
      numeric_ids <- numeric_ids[nchar(numeric_ids) > 4]
      
      if(length(numeric_ids) > 0) {
        cat("Trying to map numeric IDs as Entrez...\n")
        mapped <- map_to_symbols(numeric_ids, "entrez")
        if(!is.null(mapped)) {
          # Match back to original IDs
          result <- data.frame(
            SYMBOL = mapped$SYMBOL,
            ORIGINAL_ID = gene_ids[match(mapped$ENTREZID, numeric_ids)],
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
        mapped$ORIGINAL_ID <- mapped$PROBEID
        mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
        return(mapped)
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
      return(mapped[, c("SYMBOL", "ORIGINAL_ID")])
    }
  }
  
  cat("Warning: Could not map probe IDs\n")
  return(NULL)
}

# Main processing function - FIXED STRING COLUMN HANDLING
process_files <- function(file_paths, string_network_path, output_dir = "./results") {
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Clean file paths
  string_network_path <- clean_file_path(string_network_path)
  
  # Load STRING network
  cat(paste0("Loading STRING network from: ", string_network_path, "\n"))
  
  # Check if file exists
  if(!file.exists(string_network_path)) {
    stop(paste0("STRING network file not found: ", string_network_path))
  }
  
  # Read STRING network
  string_network <- fread(string_network_path)
  cat(paste0("STRING network loaded: ", nrow(string_network), " interactions\n"))
  
  # Show column names for debugging
  cat("Columns in STRING network:\n")
  print(colnames(string_network))
  
  # FIXED: Use the correct columns for gene symbols
  # Based on your description, the gene symbols are in '#node1' and 'node2'
  # Remove the '#' from '#node1' for easier handling
  colnames(string_network) <- gsub("^#", "", colnames(string_network))
  
  # Check if we have the expected columns
  if(all(c("node1", "node2") %in% colnames(string_network))) {
    # Rename to gene1 and gene2 for consistency
    colnames(string_network)[colnames(string_network) == "node1"] <- "gene1"
    colnames(string_network)[colnames(string_network) == "node2"] <- "gene2"
    cat("Successfully identified gene symbol columns: node1 -> gene1, node2 -> gene2\n")
  } else {
    # Fallback: look for any columns that might contain gene symbols
    cat("Looking for gene symbol columns...\n")
    
    # Check each column for gene symbol patterns
    gene_cols <- c()
    for(col in colnames(string_network)) {
      # Sample some values
      sample_vals <- head(string_network[[col]], 20)
      # Check if they look like gene symbols
      if(all(grepl("^[A-Z][A-Z0-9-]+$", sample_vals) & nchar(sample_vals) > 1 & nchar(sample_vals) < 20)) {
        gene_cols <- c(gene_cols, col)
      }
    }
    
    if(length(gene_cols) >= 2) {
      colnames(string_network)[colnames(string_network) == gene_cols[1]] <- "gene1"
      colnames(string_network)[colnames(string_network) == gene_cols[2]] <- "gene2"
      cat(paste0("Found gene columns: '", gene_cols[1], "' -> 'gene1', '", gene_cols[2], "' -> 'gene2'\n"))
    } else {
      stop("Could not identify gene symbol columns in STRING network.")
    }
  }
  
  # Show first few rows to verify
  cat("\nFirst few interactions in STRING network:\n")
  print(head(string_network[, c("gene1", "gene2")]))
  
  # Process each file
  all_symbols <- c()
  file_results <- list()
  
  for(i in seq_along(file_paths)) {
    file_path <- clean_file_path(file_paths[i])
    file_name <- basename(file_path)
    
    cat(paste0("\n", strrep("=", 60), "\n"))
    cat(paste0("Processing file ", i, "/", length(file_paths), ": ", file_name, "\n"))
    cat(paste0(strrep("=", 60), "\n"))
    
    # Check if file exists
    if(!file.exists(file_path)) {
      cat(paste0("Warning: File not found: ", file_path, "\nSkipping...\n"))
      next
    }
    
    # Read the file
    tryCatch({
      # Try reading with auto-detection
      cat("Reading file...\n")
      
      # Check file size first
      file_size <- file.info(file_path)$size
      if(file_size > 100000000) { # > 100MB
        cat("Large file detected. Reading first 100,000 rows...\n")
        de_data <- fread(file_path, nrows = 100000)
      } else {
        de_data <- fread(file_path)
      }
      
      cat(paste0("Successfully read file. Columns: ", paste(colnames(de_data), collapse = ", "), "\n"))
      cat(paste0("Rows: ", nrow(de_data), "\n"))
      
      # Get gene IDs from first column
      gene_ids <- de_data[[1]]
      cat(paste0("Total genes in file: ", length(gene_ids), "\n"))
      cat(paste0("Sample IDs: ", paste(head(gene_ids, 10), collapse = ", "), "\n"))
      
      # Identify annotation format
      format_type <- identify_annotation_format(gene_ids)
      
      # Map to symbols
      mapped_symbols <- map_to_symbols(gene_ids, format_type, file_name)
      
      if(is.null(mapped_symbols) || nrow(mapped_symbols) == 0) {
        cat("Warning: No symbols could be mapped for this file\n")
        
        # Try alternative: check other columns
        cat("Checking other columns for gene symbols...\n")
        for(col_idx in 2:min(6, ncol(de_data))) {
          col_name <- colnames(de_data)[col_idx]
          test_values <- head(de_data[[col_idx]], 50)
          
          # Count valid gene symbols
          valid_symbols <- sum(grepl("^[A-Za-z][A-Za-z0-9-]*$", test_values) & 
                                nchar(test_values) > 1 & nchar(test_values) < 30 & 
                                !grepl("^\\d+$", test_values), na.rm = TRUE)
          
          if(valid_symbols > 10) {
            cat(paste0("Found ", valid_symbols, " valid symbols in column '", col_name, "'\n"))
            mapped_symbols <- data.frame(
              SYMBOL = de_data[[col_idx]], 
              ORIGINAL_ID = de_data[[1]],
              stringsAsFactors = FALSE
            )
            # Clean up
            mapped_symbols$SYMBOL[!grepl("^[A-Za-z][A-Za-z0-9-]*$", mapped_symbols$SYMBOL)] <- NA
            mapped_symbols$SYMBOL[nchar(mapped_symbols$SYMBOL) < 2 | nchar(mapped_symbols$SYMBOL) > 30] <- NA
            mapped_symbols <- mapped_symbols[!is.na(mapped_symbols$SYMBOL) & mapped_symbols$SYMBOL != "", ]
            break
          }
        }
      }
      
      if(is.null(mapped_symbols) || nrow(mapped_symbols) == 0) {
        cat("Skipping this file - no mappable gene symbols found\n")
        next
      }
      
      # Clean mapped symbols
      mapped_symbols_clean <- unique(mapped_symbols$SYMBOL[!is.na(mapped_symbols$SYMBOL) & 
                                                             mapped_symbols$SYMBOL != ""])
      cat(paste0("Successfully mapped to ", length(mapped_symbols_clean), " unique gene symbols\n"))
      
      # Add to total symbols list
      all_symbols <- unique(c(all_symbols, mapped_symbols_clean))
      
      # Save mapped results
      safe_name <- gsub("[^A-Za-z0-9._-]", "_", file_name)
      output_file <- file.path(output_dir, paste0("mapped_", sub("\\.(csv|tsv|txt)$", "", safe_name), ".tsv"))
      write.table(mapped_symbols, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat(paste0("Mapped results saved to: ", output_file, "\n"))
      
      # Store for later use
      file_results[[file_name]] <- list(
        original_count = length(gene_ids),
        mapped_count = length(mapped_symbols_clean),
        mapped_symbols = mapped_symbols_clean,
        format_type = format_type
      )
      
    }, error = function(e) {
      cat(paste0("Error processing file ", file_name, ": ", e$message, "\n"))
    })
  }
  
  cat(paste0("\n", strrep("=", 60), "\n"))
  cat("SUMMARY\n")
  cat(paste0(strrep("=", 60), "\n"))
  cat(paste0("Total unique gene symbols across all files: ", length(all_symbols), "\n"))
  
  # Filter STRING network
  if(length(all_symbols) > 0 && all(c("gene1", "gene2") %in% colnames(string_network))) {
    cat("\nFiltering STRING network...\n")
    cat(paste0("Gene symbols available: ", length(all_symbols), "\n"))
    cat(paste0("Sample symbols: ", paste(head(sort(all_symbols)), collapse = ", "), "\n"))
    
    # Show what's in the STRING network
    cat("Sample from STRING network gene1 column: ", paste(head(unique(string_network$gene1)), collapse = ", "), "\n")
    cat("Sample from STRING network gene2 column: ", paste(head(unique(string_network$gene2)), collapse = ", "), "\n")
    
    # Check for overlap
    overlap_gene1 <- sum(string_network$gene1 %in% all_symbols)
    overlap_gene2 <- sum(string_network$gene2 %in% all_symbols)
    cat(paste0("Overlap - gene1 in our symbols: ", overlap_gene1, " (", 
               round(100 * overlap_gene1 / nrow(string_network), 1), "%)\n"))
    cat(paste0("Overlap - gene2 in our symbols: ", overlap_gene2, " (", 
               round(100 * overlap_gene2 / nrow(string_network), 1), "%)\n"))
    
    # Filter for interactions where both genes are in our gene list
    filtered_network <- string_network %>%
      filter(gene1 %in% all_symbols & gene2 %in% all_symbols)
    
    cat(paste0("Filtered STRING network: ", nrow(filtered_network), 
               " interactions (", round(100 * nrow(filtered_network) / nrow(string_network), 2), 
               "% of original)\n"))
    
    if(nrow(filtered_network) > 0) {
      # Save filtered network
      filtered_network_file <- file.path(output_dir, "filtered_string_network.tsv")
      write.table(filtered_network, filtered_network_file, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
      cat(paste0("Filtered STRING network saved to: ", filtered_network_file, "\n"))
      
      # Also save a simple edge list
      edge_list <- filtered_network[, c("gene1", "gene2")]
      edge_list_file <- file.path(output_dir, "ppi_edge_list.tsv")
      write.table(edge_list, edge_list_file, sep = "\t", 
                  row.names = FALSE, quote = FALSE, col.names = FALSE)
      cat(paste0("PPI edge list saved to: ", edge_list_file, "\n"))
      
      # Count unique genes in the filtered network
      unique_genes_filtered <- unique(c(filtered_network$gene1, filtered_network$gene2))
      cat(paste0("Unique genes in filtered network: ", length(unique_genes_filtered), 
                 " (", round(100 * length(unique_genes_filtered) / length(all_symbols), 1), "% of input genes)\n"))
    } else {
      cat("Warning: No interactions found in STRING network for these genes.\n")
      cat("This could mean:\n")
      cat("1. The gene symbols in your data don't match STRING's symbols\n")
      cat("2. There are truly no interactions between these genes in STRING\n")
      
      # Try case-insensitive matching
      cat("\nTrying case-insensitive matching...\n")
      string_genes_upper <- toupper(string_network$gene1)
      our_symbols_upper <- toupper(all_symbols)
      
      filtered_network_case <- string_network %>%
        filter(toupper(gene1) %in% our_symbols_upper & toupper(gene2) %in% our_symbols_upper)
      
      if(nrow(filtered_network_case) > 0) {
        cat(paste0("Found ", nrow(filtered_network_case), " interactions with case-insensitive matching\n"))
        filtered_network_file <- file.path(output_dir, "filtered_string_network_case_insensitive.tsv")
        write.table(filtered_network_case, filtered_network_file, sep = "\t", 
                    row.names = FALSE, quote = FALSE)
        cat(paste0("Case-insensitive filtered network saved to: ", filtered_network_file, "\n"))
      }
    }
    
  } else {
    cat("\nCannot filter STRING network.\n")
    if(!all(c("gene1", "gene2") %in% colnames(string_network))) {
      cat("Missing gene1 or gene2 columns in STRING network.\n")
      cat("Available columns: ", paste(colnames(string_network), collapse = ", "), "\n")
    }
    if(length(all_symbols) == 0) {
      cat("No gene symbols were successfully mapped.\n")
    }
  }
  
  # Save the list of all symbols
  symbols_file <- file.path(output_dir, "all_gene_symbols.txt")
  writeLines(sort(all_symbols), symbols_file)
  cat(paste0("Gene symbols list saved to: ", symbols_file, "\n"))
  
  # Save summary statistics
  summary_file <- file.path(output_dir, "mapping_summary.tsv")
  summary_df <- data.frame(
    File = names(file_results),
    Original_Genes = sapply(file_results, function(x) x$original_count),
    Mapped_Symbols = sapply(file_results, function(x) x$mapped_count),
    Format = sapply(file_results, function(x) x$format_type),
    stringsAsFactors = FALSE
  )
  write.table(summary_df, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste0("Mapping summary saved to: ", summary_file, "\n"))
  
  return(list(
    file_results = file_results,
    all_symbols = all_symbols,
    filtered_network = if(exists("filtered_network") && nrow(filtered_network) > 0) filtered_network else NULL
  ))
}

# Direct execution with your files
cat("Running STRING PPI filtering script...\n")

# Define your files
string_path <- "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\datasets\\ppi\\string_interactions.tsv"
de_files <- c(
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE833\\GSE833_DE_results.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE76220\\GSE76220_DE_results.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE112676\\GSE112676_DE_results.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE118336\\GSE118336_DE_FUS_H517D_Mutant_vs_Control.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE118336\\GSE118336_DE_FUS_Heterozygous_vs_Control.csv",
  "C:\\Users\\noahm\\PycharmProjects\\MarbleProject\\de_results\\GSE124439\\GSE124439_DE_results.csv"
)

output_dir <- "./ppi_results"

# Run the analysis
results <- process_files(de_files, string_path, output_dir)

cat("\nAnalysis complete!\n")