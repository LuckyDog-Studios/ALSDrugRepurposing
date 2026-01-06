#!/usr/bin/env Rscript
# Annotation Format Checker for GSE Differential Expression Results
# Handles unnamed columns (like column0 in Positron)

# Load required libraries
required_packages <- c("dplyr", "stringr", "data.table")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(stringr)
library(data.table)

# Function to check annotation format
check_annotation_format <- function(file_path, column_index = 1) {
  
  cat("Reading file:", file_path, "\n")
  
  # Read first few lines to understand structure
  cat("\n=== File Structure Analysis ===\n")
  
  # Try to read with different approaches
  try_read <- function(file_path) {
    # Read first 5 lines to understand structure
    first_lines <- readLines(file_path, n = 10)
    cat("First 3 lines of the file:\n")
    for(i in seq_len(min(3, length(first_lines)))) {
      cat(i, ":", first_lines[i], "\n")
    }
    
    # Detect separator
    line1 <- first_lines[1]
    tab_count <- str_count(line1, "\t")
    comma_count <- str_count(line1, ",")
    semicolon_count <- str_count(line1, ";")
    
    sep <- ifelse(tab_count > comma_count & tab_count > semicolon_count, "\t",
                  ifelse(comma_count > semicolon_count, ",", ";"))
    
    cat("\nDetected separator:", 
        ifelse(sep == "\t", "TAB", 
               ifelse(sep == ",", "COMMA", "SEMICOLON")), "\n")
    
    # Try reading with data.table (more robust)
    if(sep == "\t") {
      data <- try(fread(file_path, sep = "\t", header = TRUE, fill = TRUE), silent = TRUE)
    } else if(sep == ",") {
      data <- try(fread(file_path, sep = ",", header = TRUE, fill = TRUE), silent = TRUE)
    } else {
      data <- try(fread(file_path, sep = ";", header = TRUE, fill = TRUE), silent = TRUE)
    }
    
    # If fread fails, try base R
    if(inherits(data, "try-error")) {
      cat("fread failed, trying base R...\n")
      data <- try(read.table(file_path, header = TRUE, sep = sep, 
                            stringsAsFactors = FALSE, check.names = FALSE,
                            fill = TRUE, quote = "", comment.char = ""), silent = TRUE)
    }
    
    return(data)
  }
  
  data <- try_read(file_path)
  
  if(inherits(data, "try-error")) {
    # Try reading without header
    cat("\nTrying without header...\n")
    data <- try(fread(file_path, header = FALSE, fill = TRUE), silent = TRUE)
    if(inherits(data, "try-error")) {
      stop("Could not read file. Please check if it's a valid text file.")
    }
    colnames(data) <- paste0("V", 1:ncol(data))
  }
  
  cat("\nFile read successfully!\n")
  cat("Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
  
  # Show column names
  cat("\nColumn names:\n")
  print(colnames(data))
  
  # If column_index is numeric, use it; otherwise try to find column
  if(is.numeric(column_index)) {
    if(column_index > ncol(data)) {
      stop("Column index ", column_index, " exceeds number of columns (", ncol(data), ")")
    }
    target_col <- column_index
    cat("\nUsing column index:", column_index, "\n")
  } else {
    # Try to find column by name pattern
    target_col <- grep(column_index, colnames(data), ignore.case = TRUE)
    if(length(target_col) == 0) {
      cat("\nColumn '", column_index, "' not found by name. Using first column.\n")
      target_col <- 1
    } else {
      target_col <- target_col[1]
    }
  }
  
  # Extract annotations
  annotations <- data[[target_col]]
  col_name <- colnames(data)[target_col]
  
  cat("\nAnalyzing column:", col_name, "(index:", target_col, ")\n")
  cat("First 10 values in this column:\n")
  print(head(annotations, 10))
  
  # Clean annotations
  annotations <- annotations[!is.na(annotations) & annotations != "" & !is.null(annotations)]
  
  if(length(annotations) == 0) {
    stop("No annotations found in the specified column")
  }
  
  cat("\n=== Annotation Format Analysis ===\n")
  cat("Total annotations checked:", length(annotations), "\n")
  
  # Define pattern tests for different annotation formats
  format_tests <- list(
    "AFFYMETRIX_PROBE_ID" = list(
      pattern = "^\\d+_([asxibf])_at$|^AFFX-|^ILMN_|^[A-Za-z]{2,3}\\d{5,}_[asx]_at$",
      description = "Affymetrix/Illumina probe IDs",
      example = c("1007_s_at", "AFFX-r2-Bs-dap-3_at", "ILMN_12345")
    ),
    "ENTREZ_GENE_ID" = list(
      pattern = "^\\d+$",
      description = "Entrez Gene IDs (numeric only)",
      example = c("1", "7157", "100008586")
    ),
    "ENSEMBL_GENE_ID" = list(
      pattern = "^ENS[A-Z]{0,3}G\\d{11}$",
      description = "Ensembl Gene IDs",
      example = c("ENSG00000139618", "ENSMUSG00000000001")
    ),
    "ENSEMBL_TRANSCRIPT_ID" = list(
      pattern = "^ENS[A-Z]{0,3}T\\d{11}$",
      description = "Ensembl Transcript IDs",
      example = c("ENST00000380152", "ENSMUST00000000001")
    ),
    "REFSEQ_RNA" = list(
      pattern = "^N[MR]_\\d+(\\.\\d+)?$|^X[MR]_\\d+(\\.\\d+)?$",
      description = "RefSeq RNA accessions",
      example = c("NM_001126114", "NR_003286", "XM_005249643")
    ),
    "REFSEQ_PROTEIN" = list(
      pattern = "^N[P]_\\d+(\\.\\d+)?$|^X[P]_\\d+(\\.\\d+)?$",
      description = "RefSeq Protein accessions",
      example = c("NP_001119", "XP_005249700")
    ),
    "GENE_SYMBOL" = list(
      pattern = "^[A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]$",
      description = "Gene Symbols",
      example = c("TP53", "BRCA1", "IL-6", "MIR21")
    ),
    "AGILENT_PROBE_ID" = list(
      pattern = "^A_\\d+_P\\d+$|^GE_BrightCorner",
      description = "Agilent probe IDs",
      example = c("A_23_P117082", "GE_BrightCorner")
    ),
    "FLYBASE" = list(
      pattern = "^FB[a-z]{2}\\d+$",
      description = "FlyBase IDs",
      example = c("FBgn0031208")
    ),
    "UNIPROT" = list(
      pattern = "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]",
      description = "UniProt accessions",
      example = c("P12345", "A0A0A0MZQ8")
    ),
    "GENBANK_GI" = list(
      pattern = "^gi\\|\\d+",
      description = "GenBank GI numbers",
      example = c("gi|1234567")
    ),
    "GENBANK_ACCESSION" = list(
      pattern = "^[A-Z]{1,6}_?\\d+(\\.\\d+)?$",
      description = "GenBank accessions",
      example = c("BC053857", "AK123456", "AY123456.1")
    ),
    "ILLUMINA_PROBE_ID" = list(
      pattern = "^ILMN_\\d+$|^GI_\\d+_[-\\w]+",
      description = "Illumina probe IDs",
      example = c("ILMN_1802380", "GI_10047089-S")
    ),
    "UCSC_ID" = list(
      pattern = "^uc[0-9a-z]{3}\\.\\d+$",
      description = "UCSC IDs",
      example = c("uc001aaa.3")
    ),
    "MIRBASE" = list(
      pattern = "^hsa-mir|^mmu-mir|^miR-",
      description = "miRBase miRNA IDs",
      example = c("hsa-mir-21", "miR-21-5p")
    ),
    "LOCUS_ID" = list(
      pattern = "^LOC\\d+$|^AT[0-9CM]G\\d{5}$",
      description = "Locus IDs / TAIR IDs",
      example = c("LOC100288142", "AT1G01010")
    )
  )
  
  # Test each format
  format_results <- data.frame()
  sample_matches <- list()
  
  for(format_name in names(format_tests)) {
    test <- format_tests[[format_name]]
    matches <- grepl(test$pattern, annotations)
    match_count <- sum(matches, na.rm = TRUE)
    match_percent <- round(match_count / length(annotations) * 100, 2)
    
    # Get sample matches
    if(match_count > 0) {
      matched_annots <- annotations[matches]
      sample_matches[[format_name]] <- head(matched_annots, 5)
    }
    
    format_results <- rbind(format_results, data.frame(
      Format = format_name,
      Description = test$description,
      Matches = match_count,
      Percentage = match_percent,
      stringsAsFactors = FALSE
    ))
  }
  
  # Add a check for mixed formats
  format_results$IsMixed <- FALSE
  
  # Sort by match percentage
  format_results <- format_results[order(-format_results$Percentage), ]
  
  # Print results in a nice format
  cat("\nFormat Detection Results:\n")
  cat(rep("-", 80), "\n", sep = "")
  cat(sprintf("%-25s %-30s %-10s %-10s\n", 
              "Format", "Description", "Matches", "Percent"))
  cat(rep("-", 80), "\n", sep = "")
  
  for(i in 1:nrow(format_results)) {
    row <- format_results[i, ]
    cat(sprintf("%-25s %-30s %-10d %-9.1f%%\n", 
                row$Format, row$Description, row$Matches, row$Percentage))
  }
  cat(rep("-", 80), "\n\n", sep = "")
  
  # Identify most likely format
  top_format <- format_results$Format[1]
  top_percent <- format_results$Percentage[1]
  
  cat("=== SUMMARY ===\n")
  cat("Most likely format:", top_format, "\n")
  cat("Match confidence:", top_percent, "%\n\n")
  
  if(top_percent > 90) {
    cat("✅ HIGH confidence: This format appears to be correct.\n")
  } else if(top_percent > 70) {
    cat("⚠️  MEDIUM confidence: This might be the format, but check manually.\n")
  } else if(top_percent > 30) {
    cat("⚠️  LOW confidence: Multiple formats detected. Manual inspection needed.\n")
  } else {
    cat("❓ UNKNOWN format: No clear pattern detected.\n")
  }
  
  # Show top 3 candidate formats with examples
  if(top_percent < 90) {
    cat("\nTop candidate formats:\n")
    for(i in 1:min(3, nrow(format_results))) {
      if(format_results$Percentage[i] > 0) {
        cat(i, ". ", format_results$Format[i], " (", 
            format_results$Percentage[i], "%)\n", sep = "")
        if(!is.null(sample_matches[[format_results$Format[i]]])) {
          cat("   Examples: ", 
              paste(sample_matches[[format_results$Format[i]]], collapse = ", "), "\n")
        }
      }
    }
  }
  
  # Show examples of top format
  if(!is.null(sample_matches[[top_format]])) {
    cat("\nExamples of ", top_format, ":\n", sep = "")
    print(sample_matches[[top_format]])
  }
  
  # Additional diagnostics
  cat("\n=== ADDITIONAL DIAGNOSTICS ===\n")
  
  # Check for version numbers
  versioned <- grepl("\\.\\d+$", annotations)
  cat("Annotations with version numbers (e.g., .1, .2):", 
      sum(versioned), paste0("(", round(sum(versioned)/length(annotations)*100, 1), "%)\n"))
  
  # Check for special characters
  special_chars <- grepl("[^A-Za-z0-9._-]", annotations)
  cat("Annotations with special characters:", 
      sum(special_chars), paste0("(", round(sum(special_chars)/length(annotations)*100, 1), "%)\n"))
  
  # Check lengths
  annot_lengths <- nchar(annotations)
  cat("Average annotation length:", round(mean(annot_lengths), 1), "characters\n")
  cat("Length range:", min(annot_lengths), "-", max(annot_lengths), "characters\n")
  
  # Check for duplicates
  unique_count <- length(unique(annotations))
  cat("Unique annotations:", unique_count, 
      paste0("(", round(unique_count/length(annotations)*100, 1), "% unique)\n"))
  
  # Return results
  return(list(
    data = data,
    annotations = annotations,
    column_name = col_name,
    column_index = target_col,
    format_results = format_results,
    top_format = top_format,
    top_percent = top_percent,
    sample_matches = sample_matches
  ))
}

# Helper function to explore data
explore_data <- function(data) {
  cat("\n=== Data Exploration ===\n")
  cat("Number of rows:", nrow(data), "\n")
  cat("Number of columns:", ncol(data), "\n\n")
  
  cat("Column names:\n")
  for(i in seq_along(colnames(data))) {
    cat(i, ": '", colnames(data)[i], "'\n", sep = "")
  }
  
  cat("\nFirst few rows of each column:\n")
  for(i in 1:min(5, ncol(data))) {
    cat("\nColumn", i, "('", colnames(data)[i], "'):\n", sep = "")
    print(head(data[[i]], 5))
  }
}

# Main execution function
main <- function() {
  cat("GSE Annotation Format Checker\n")
  cat("=============================\n\n")
  
  # Get file path
  if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    file_path <- commandArgs(trailingOnly = TRUE)[1]
  } else {
    cat("Enter the path to your differential expression file:\n")
    file_path <- readline("File path: ")
  }
  
  # Check if file exists
  if(!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # First, explore the data structure
  cat("\nFirst, let's explore the file structure...\n")
  
  # Try to read and show structure
  tryCatch({
    # Read first few lines
    first_lines <- readLines(file_path, n = 5)
    cat("\nFirst few lines of file:\n")
    for(i in seq_along(first_lines)) {
      cat(i, ": ", first_lines[i], "\n", sep = "")
    }
    
    # Ask user which column to check
    cat("\nWhich column contains the annotations?\n")
    cat("Options:\n")
    cat("1. Enter column NUMBER (e.g., 1 for first column)\n")
    cat("2. Enter column NAME if you know it\n")
    cat("3. Enter 'explore' to see all columns first\n")
    
    user_input <- readline("Your choice: ")
    
    if(tolower(user_input) == "explore") {
      # Read the data
      data <- fread(file_path, fill = TRUE)
      explore_data(data)
      
      user_input <- readline("\nNow enter column number or name: ")
    }
    
    # Convert to numeric if possible
    column_spec <- ifelse(grepl("^\\d+$", user_input), 
                         as.numeric(user_input), 
                         user_input)
    
    # Run the analysis
    results <- check_annotation_format(file_path, column_spec)
    
    # Offer to save results
    cat("\nWould you like to save these results? (y/n): ")
    save_choice <- readline()
    
    if(tolower(save_choice) %in% c("y", "yes")) {
      output_file <- paste0(tools::file_path_sans_ext(basename(file_path)), 
                           "_annotation_report.txt")
      
      sink(output_file)
      cat("Annotation Format Analysis Report\n")
      cat("=================================\n\n")
      cat("File:", file_path, "\n")
      cat("Date:", date(), "\n\n")
      
      cat("Column analyzed:", results$column_name, "(index:", results$column_index, ")\n")
      cat("Total annotations:", length(results$annotations), "\n")
      cat("Unique annotations:", length(unique(results$annotations)), "\n\n")
      
      cat("Format Detection Results:\n")
      print(results$format_results)
      cat("\nMost likely format:", results$top_format, "\n")
      cat("Confidence:", results$top_percent, "%\n\n")
      
      cat("Sample annotations (first 20):\n")
      print(head(results$annotations, 20))
      
      sink()
      cat("\nReport saved to:", output_file, "\n")
    }
    
    # Suggest next steps
    cat("\n=== NEXT STEPS ===\n")
    cat("Based on the detected format, you may need to:\n")
    cat("1. Map identifiers using appropriate annotation packages\n")
    cat("2. Use biomaRt for ID conversion\n")
    cat("3. Check GEO for platform annotation files\n")
    
  }, error = function(e) {
    cat("\n❌ Error:", e$message, "\n")
    cat("\nTroubleshooting tips:\n")
    cat("1. Make sure the file is a plain text file (TSV, CSV)\n")
    cat("2. Check if the file has headers\n")
    cat("3. Try opening the file in a text editor first\n")
  })
}

# Run if called from command line
if(!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
} else {
  # For interactive use in RStudio/Positron
  cat("To use this script:\n")
  cat("1. Source this file: source('annotation_checker.R')\n")
  cat("2. Run: main()\n")
  cat("3. Or call directly: check_annotation_format('your_file.tsv', 1)\n")
}