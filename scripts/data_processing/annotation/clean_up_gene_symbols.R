#!/usr/bin/env Rscript

# Script to clean and validate gene symbols in combined DE results

library(data.table)
library(dplyr)

# Function to clean gene symbols
clean_gene_symbol <- function(symbol) {
  if(is.na(symbol) || symbol == "") return(NA)
  
  # Remove whitespace
  symbol <- trimws(symbol)
  
  # Handle concatenated symbols (take first valid one)
  if(grepl(";", symbol)) {
    parts <- strsplit(symbol, ";")[[1]]
    parts <- trimws(parts)
    
    # Filter for valid gene symbols
    valid_parts <- parts[sapply(parts, function(x) is_valid_gene_symbol(x))]
    
    if(length(valid_parts) > 0) {
      # Prefer symbols without "LOC", "LINC", etc.
      preferred <- valid_parts[!grepl("^LOC\\d+|^LINC\\d+|^RP\\d+-|^MIR\\d+|^NCRNA|^OTTHUMG", valid_parts)]
      if(length(preferred) > 0) {
        return(preferred[1])
      }
      return(valid_parts[1])
    }
    return(NA)
  }
  
  # Remove known problematic patterns
  if(grepl("Satellite|:RNA:|^\\(|RNA$|pseudogene$", symbol, ignore.case = TRUE)) {
    return(NA)
  }
  
  # Remove version numbers (e.g., .1, .2)
  symbol <- gsub("\\.\\d+$", "", symbol)
  
  # Remove transcript identifiers
  symbol <- gsub("-\\d+$", "", symbol)
  
  # Convert to uppercase (standard gene symbols are uppercase)
  symbol <- toupper(symbol)
  
  if(is_valid_gene_symbol(symbol)) {
    return(symbol)
  }
  
  return(NA)
}

# Function to validate gene symbols
is_valid_gene_symbol <- function(symbol) {
  if(is.na(symbol) || symbol == "") return(FALSE)
  
  # Common artifacts to exclude
  artifacts <- c("NA", "NONE", "UNDEFINED", "UNKNOWN", "CONTROL", "BLANK", "EMPTY",
                 "CHR", "CHROMOSOME", "SCAFFOLD", "CONTIG", "LOCUS")
  
  if(symbol %in% artifacts) return(FALSE)
  
  # Pattern checks
  pattern_ok <- grepl("^[A-Z][A-Z0-9-]*[A-Z0-9]$", symbol) && 
                nchar(symbol) >= 2 && 
                nchar(symbol) <= 30
  
  # Exclude specific problematic patterns
  not_artifactual <- !grepl("^\\d+$|^[A-Z]{1}$|^[A-Z]{15,}$|^[A-Z]+\\d+[A-Z]+$|^\\d+[A-Z]+\\d+$", symbol)
  
  # Exclude common non-gene annotations
  not_other <- !grepl("^TRN[A-Z]-|^IGH|^IG[LK]|^TRA|^TRB|^TRD|^TRG|^RN[A-Z]|^SNOR|^MIR|^LINCRNA", symbol)
  
  return(pattern_ok && not_artifactual && not_other)
}

# Function to clean the entire DE results file
clean_combined_de_results <- function(input_file, output_file = NULL) {
  
  cat("Loading combined DE results...\n")
  de_data <- fread(input_file)
  
  cat(paste0("Original data: ", nrow(de_data), " rows\n"))
  cat("Columns:\n")
  print(colnames(de_data))
  
  # Check if SYMBOL column exists
  if(!"SYMBOL" %in% colnames(de_data)) {
    stop("SYMBOL column not found in the data")
  }
  
  # Show original symbol statistics
  original_symbols <- de_data$SYMBOL
  unique_original <- unique(original_symbols)
  cat(paste0("\nOriginal unique symbols: ", length(unique_original), "\n"))
  
  # Diagnostic of problematic symbols
  cat("\nDiagnostic of original symbols:\n")
  
  patterns <- list(
    concatenated = sum(grepl(";", original_symbols)),
    satellite = sum(grepl("Satellite", original_symbols, ignore.case = TRUE)),
    parentheses = sum(grepl("\\(", original_symbols)),
    rna_notation = sum(grepl(":RNA:|RNA$", original_symbols, ignore.case = TRUE)),
    looks_like_symbol = sum(sapply(original_symbols, is_valid_gene_symbol))
  )
  
  for(name in names(patterns)) {
    cat(sprintf("  %-20s: %6d (%5.1f%%)\n", 
                name, 
                patterns[[name]], 
                100 * patterns[[name]] / length(original_symbols)))
  }
  
  # Show examples of problematic symbols
  cat("\nExamples of problematic symbols:\n")
  problematic_indices <- which(!sapply(original_symbols, is_valid_gene_symbol))
  if(length(problematic_indices) > 0) {
    problematic <- head(original_symbols[problematic_indices], 10)
    print(problematic)
  }
  
  # Clean symbols
  cat("\nCleaning symbols...\n")
  cleaned_symbols <- sapply(original_symbols, clean_gene_symbol)
  
  # Add cleaned symbols as new column
  de_data$SYMBOL_CLEANED <- cleaned_symbols
  
  # Statistics after cleaning
  valid_rows <- !is.na(de_data$SYMBOL_CLEANED)
  cat(paste0("\nRows with valid cleaned symbols: ", sum(valid_rows), 
             " (", round(100 * sum(valid_rows) / nrow(de_data), 1), "%)\n"))
  
  # Show symbol changes
  changed_symbols <- which(de_data$SYMBOL != de_data$SYMBOL_CLEANED & !is.na(de_data$SYMBOL_CLEANED))
  cat(paste0("Symbols that were changed: ", length(changed_symbols), 
             " (", round(100 * length(changed_symbols) / sum(valid_rows), 1), "%)\n"))
  
  # Show examples of cleaned symbols
  if(length(changed_symbols) > 0) {
    cat("\nExamples of symbol cleaning:\n")
    n_examples <- min(10, length(changed_symbols))
    examples <- data.frame(
      Original = de_data$SYMBOL[changed_symbols[1:n_examples]],
      Cleaned = de_data$SYMBOL_CLEANED[changed_symbols[1:n_examples]],
      stringsAsFactors = FALSE
    )
    print(examples)
  }
  
  # Create cleaned dataset (only rows with valid symbols)
  cleaned_data <- de_data[valid_rows, ]
  
  # Remove the original SYMBOL column and rename cleaned column
  cleaned_data$SYMBOL <- NULL
  colnames(cleaned_data)[colnames(cleaned_data) == "SYMBOL_CLEANED"] <- "SYMBOL"
  
  # Summary by number of datasets
  if("n_datasets" %in% colnames(cleaned_data)) {
    cat("\nCleaned data summary by number of datasets:\n")
    summary_by_datasets <- cleaned_data[, .(count = .N, unique_genes = length(unique(SYMBOL))), 
                                        by = n_datasets][order(-n_datasets)]
    print(summary_by_datasets)
  }
  
  # Summary of unique cleaned symbols
  unique_cleaned <- unique(cleaned_data$SYMBOL)
  cat(paste0("\nUnique cleaned symbols: ", length(unique_cleaned), "\n"))
  
  # Show sample of cleaned symbols
  cat("\nSample cleaned symbols (first 20):\n")
  print(head(sort(unique_cleaned), 20))
  
  # Save cleaned data
  if(is.null(output_file)) {
    output_file <- gsub("\\.(csv|tsv|txt)$", "_cleaned.csv", input_file)
  }
  
  cat(paste0("\nSaving cleaned data to: ", output_file, "\n"))
  fwrite(cleaned_data, output_file)
  
  # Create additional summary files
  output_dir <- dirname(output_file)
  base_name <- tools::file_path_sans_ext(basename(output_file))
  
  # Save list of unique cleaned symbols
  symbols_file <- file.path(output_dir, paste0(base_name, "_symbols.txt"))
  writeLines(sort(unique_cleaned), symbols_file)
  cat(paste0("Unique symbols saved to: ", symbols_file, "\n"))
  
  # Save mapping of original to cleaned symbols
  if(length(changed_symbols) > 0) {
    mapping <- de_data[changed_symbols, .(Original = SYMBOL, Cleaned = SYMBOL_CLEANED)]
    mapping <- mapping[!duplicated(mapping), ]
    mapping_file <- file.path(output_dir, paste0(base_name, "_symbol_mapping.tsv"))
    fwrite(mapping, mapping_file, sep = "\t")
    cat(paste0("Symbol mapping saved to: ", mapping_file, "\n"))
  }
  
  # Save summary statistics
  summary_stats <- data.frame(
    Metric = c("Original rows", "Rows after cleaning", "Retention rate",
               "Original unique symbols", "Cleaned unique symbols",
               "Symbols changed", "Change rate"),
    Value = c(nrow(de_data), nrow(cleaned_data),
              round(100 * nrow(cleaned_data) / nrow(de_data), 1),
              length(unique_original), length(unique_cleaned),
              length(changed_symbols),
              round(100 * length(changed_symbols) / nrow(cleaned_data), 1))
  )
  
  summary_file <- file.path(output_dir, paste0(base_name, "_summary.tsv"))
  fwrite(summary_stats, summary_file, sep = "\t")
  cat(paste0("Summary statistics saved to: ", summary_file, "\n"))
  
  return(list(
    cleaned_data = cleaned_data,
    summary = summary_stats,
    symbol_mapping = if(exists("mapping")) mapping else NULL
  ))
}

# SIMPLE TEST - bypass the problematic test_cleaning function
cat("=== SIMPLE TEST ===\n")
test_symbols <- c("ACTB", "61E3.4;LOC100132247", "(CATTC)n:Satellite:Satellite", "SOD1")
cat("Testing a few symbols:\n")
for(sym in test_symbols) {
  cat(sprintf("  %-40s -> %s\n", sym, clean_gene_symbol(sym)))
}

# Main execution
input_file <- "C:/Users/noahm/PycharmProjects/MarbleProject/de_results/combined/combined_DE_results_annotated.csv"

if(file.exists(input_file)) {
  cat("\n\n=== CLEANING YOUR DE RESULTS ===\n")
  cleaned <- clean_combined_de_results(input_file)
  
  cat("\n\n=== NEXT STEPS ===\n")
  cat("1. Use the cleaned file for PPI filtering\n")
  cat("2. Expected much better overlap with ALS PPI network\n")
  cat("3. Focus on genes that appear in multiple datasets (n_datasets >= 2 or 3)\n")
} else {
  cat(paste0("\nFile not found: ", input_file, "\n"))
  cat("Please update the input_file path.\n")
}