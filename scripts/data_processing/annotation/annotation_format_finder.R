#!/usr/bin/env Rscript
# Enhanced Annotation Checker for GSE833 datasets

# Load required libraries
required_packages <- c("dplyr", "stringr", "data.table", "rentrez")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(stringr)
library(data.table)

# Function to check annotation format
check_gse833_annotation <- function(file_path, column_index = 1) {
  
  cat("GSE833 Annotation Analysis\n")
  cat("==========================\n\n")
  
  # Read the data
  cat("Reading file:", file_path, "\n")
  data <- fread(file_path, fill = TRUE)
  
  if(ncol(data) < column_index) {
    stop("Column index ", column_index, " not found. File has ", ncol(data), " columns.")
  }
  
  # Extract annotations
  annotations <- data[[column_index]]
  col_name <- colnames(data)[column_index]
  
  # Clean annotations
  annotations <- annotations[!is.na(annotations) & annotations != "" & !is.null(annotations)]
  
  if(length(annotations) == 0) {
    stop("No annotations found in the specified column")
  }
  
  cat("\n=== Data Summary ===\n")
  cat("Total rows:", nrow(data), "\n")
  cat("Annotations analyzed:", length(annotations), "\n")
  cat("Column name:", col_name, "\n\n")
  
  # Show examples
  cat("First 20 annotations:\n")
  print(head(annotations, 20))
  
  # Analyze the patterns
  cat("\n=== Pattern Analysis ===\n")
  
  # GenBank accession patterns (most common in GSE833)
  patterns <- list(
    "GenBank_accession_with_suffix" = "^[A-Z][0-9]{5,6}_.*$",  # e.g., U66559_at
    "GenBank_accession_no_suffix" = "^[A-Z][0-9]{5,6}$",       # e.g., U66559
    "RefSeq_accession" = "^N[MR]_[0-9]+$",                    # e.g., NM_123456
    "Ensembl_ID" = "^ENS[A-Z]*[GT][0-9]{11}$",               # e.g., ENSG00000123456
    "Gene_symbol" = "^[A-Z][A-Za-z0-9-]+[0-9]?$",            # e.g., TP53, BRCA1
    "Entrez_ID" = "^[0-9]+$",                                # e.g., 7157
    "Affy_probe_id" = "^[0-9]+_[asxibf]_at$",                # e.g., 1007_s_at
    "Illumina_probe" = "^ILMN_[0-9]+$",                      # e.g., ILMN_12345
    "Agilent_probe" = "^A_[0-9]+_P[0-9]+$",                  # e.g., A_23_P117082
    "Other_accession" = "^[A-Z]{2}[0-9]{5,8}_.*$"            # e.g., AB000816_s_at
  )
  
  # Count matches for each pattern
  results <- data.frame()
  for(pattern_name in names(patterns)) {
    matches <- grepl(patterns[[pattern_name]], annotations)
    match_count <- sum(matches)
    match_percent <- round(match_count / length(annotations) * 100, 2)
    
    results <- rbind(results, data.frame(
      Pattern = pattern_name,
      Matches = match_count,
      Percentage = match_percent,
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by match percentage
  results <- results[order(-results$Percentage), ]
  
  cat("\nPattern Match Results:\n")
  print(results)
  
  # Extract the base accession numbers (remove suffixes)
  cat("\n=== Suffix Analysis ===\n")
  suffixes <- str_extract(annotations, "_(at|s_at|f_at|xpt[0-9]_at|rna[0-9]_at|cds[0-9]_at)$")
  suffix_table <- table(suffixes)
  suffix_table <- suffix_table[order(-suffix_table)]
  
  cat("Common suffixes found:\n")
  print(head(suffix_table, 10))
  
  # Extract prefix patterns (first letters)
  cat("\n=== Prefix Analysis ===\n")
  # Remove suffixes first
  base_ids <- str_replace(annotations, "_(at|s_at|f_at|xpt[0-9]_at|rna[0-9]_at|cds[0-9]_at)$", "")
  
  # Get first 2 characters
  prefixes <- substr(base_ids, 1, 2)
  prefix_table <- table(prefixes)
  prefix_table <- prefix_table[order(-prefix_table)]
  
  cat("Most common 2-letter prefixes:\n")
  print(head(prefix_table, 10))
  
  # Check for specific GenBank accession patterns
  cat("\n=== GenBank Accession Analysis ===\n")
  
  # GenBank accession format rules:
  # 1-letter prefix + 5 digits: A12345, U12345
  # 2-letter prefix + 6 digits: AB123456, AC123456
  # 3-letter prefix + 5 digits: AAA12345
  
  prefix_lengths <- nchar(str_extract(base_ids, "^[A-Z]+"))
  cat("Prefix length distribution:\n")
  print(table(prefix_lengths))
  
  # Check if these look like valid GenBank accessions
  cat("\nAccession number length distribution:\n")
  id_lengths <- nchar(base_ids)
  print(table(id_lengths))
  
  # Get a sample of unique base IDs
  unique_base_ids <- unique(base_ids)
  cat("\nSample base IDs (without suffixes):\n")
  print(head(unique_base_ids, 20))
  
  # Determine the most likely format
  cat("\n=== Conclusion ===\n")
  if(results$Percentage[1] > 80) {
    cat("Primary format detected:", results$Pattern[1], "\n")
    cat("Confidence:", results$Percentage[1], "%\n")
  } else {
    cat("Mixed formats detected. Most common:", results$Pattern[1], "\n")
  }
  
  # Based on your data, it's clearly GenBank accession numbers with Affymetrix suffixes
  cat("\n✅ FOR GSE833: These are GenBank/EMBL/DDBJ accession numbers used as probe IDs.\n")
  cat("   Format: [GenBank Accession]_[suffix]\n")
  cat("   Examples: U66559_at, X82634_at, M95610_at\n")
  cat("   Common suffixes: _at, _s_at, _f_at, _xpt1_at, _rna1_at\n")
  
  # Return annotation mapping function
  cat("\n=== Recommended Next Steps ===\n")
  cat("1. Download the GPL (platform) annotation file from GEO\n")
  cat("2. Use biomaRt or clusterProfiler for ID conversion\n")
  cat("3. Map to current gene symbols using current databases\n")
  
  # Create a mapping helper function
  create_mapping_suggestions <- function(base_ids) {
    cat("\nSuggested mapping approaches:\n")
    cat("a) Direct mapping using GEO platform file (GPL...)\n")
    cat("b) Use rentrez package to convert GenBank to Entrez:\n")
    cat("   library(rentrez)\n")
    cat("   entrez_ids <- entrez_search(db=\"nucleotide\", term=paste(base_ids, \"[ACCN]\", sep=\"\"))\n")
    cat("c) Use biomaRt:\n")
    cat("   library(biomaRt)\n")
    cat("   mart <- useMart(\"ensembl\", dataset=\"hsapiens_gene_ensembl\")\n")
    cat("   getBM(attributes=c('entrezgene_id', 'hgnc_symbol'),\n")
    cat("         filters='refseq_mrna',\n")
    cat("         values=base_ids, mart=mart)\n")
  }
  
  create_mapping_suggestions(head(base_ids, 100))
  
  # Return analysis results
  return(list(
    data = data,
    annotations = annotations,
    base_ids = base_ids,
    suffix_distribution = suffix_table,
    prefix_distribution = prefix_table,
    pattern_results = results,
    format = "GenBank_accession_with_Affymetrix_suffix"
  ))
}

# Function to map GenBank accessions to gene symbols
map_genbank_to_genes <- function(genbank_ids, max_ids = 100) {
  cat("\n=== GenBank to Gene Mapping ===\n")
  cat("Note: This requires internet connection and may take time.\n")
  
  # Limit for demonstration
  if(length(genbank_ids) > max_ids) {
    cat("Limiting to first", max_ids, "IDs for demonstration\n")
    genbank_ids <- head(genbank_ids, max_ids)
  }
  
  # Try to use rentrez if available
  if(requireNamespace("rentrez", quietly = TRUE)) {
    library(rentrez)
    
    cat("\nAttempting to map using NCBI Entrez...\n")
    
    # First, try to get summaries
    tryCatch({
      # Search for each GenBank ID
      results <- list()
      for(i in seq_along(genbank_ids)) {
        cat("Processing", genbank_ids[i], "...\n")
        
        # Search in nucleotide database
        search_res <- try(entrez_search(db="nucleotide", term=paste0(genbank_ids[i], "[ACCN]")), 
                          silent = TRUE)
        
        if(!inherits(search_res, "try-error") && length(search_res$ids) > 0) {
          # Get summary
          summary_res <- try(entrez_summary(db="nucleotide", id=search_res$ids[1]), silent = TRUE)
          
          if(!inherits(summary_res, "try-error")) {
            results[[genbank_ids[i]]] <- list(
              accession = genbank_ids[i],
              title = ifelse(!is.null(summary_res$title), summary_res$title, NA),
              organism = ifelse(!is.null(summary_res$organism), summary_res$organism, NA),
              gi = ifelse(!is.null(summary_res$gi), summary_res$gi, NA)
            )
          }
        }
        
        # Be nice to NCBI
        Sys.sleep(0.1)
      }
      
      cat("\nSample mappings found:\n")
      for(i in 1:min(5, length(results))) {
        res <- results[[i]]
        cat(i, ". Accession:", res$accession, "\n")
        cat("   Title:", substr(res$title, 1, 50), "...\n")
        cat("   Organism:", res$organism, "\n")
        cat("   GI:", res$gi, "\n\n")
      }
      
      return(results)
      
    }, error = function(e) {
      cat("Error with Entrez:", e$message, "\n")
    })
  } else {
    cat("rentrez package not installed. Install with: install.packages('rentrez')\n")
  }
  
  # Alternative: Use pre-built mapping files
  cat("\n=== Alternative Approach ===\n")
  cat("For GSE833, download the platform annotation file from:\n")
  cat("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL80\n")
  cat("or\n")
  cat("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL91\n")
  cat("\nThese contain the mapping from GenBank accessions to genes.\n")
}

# Main execution
main <- function() {
  cat("GSE833 GenBank Accession Analyzer\n")
  cat("==================================\n\n")
  
  # Get file path
  if(length(commandArgs(trailingOnly = TRUE)) > 0) {
    file_path <- commandArgs(trailingOnly = TRUE)[1]
  } else {
    cat("Enter path to your differential expression file:\n")
    file_path <- readline("File path: ")
  }
  
  if(!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Run analysis
  results <- check_gse833_annotation(file_path, 1)
  
  # Offer to try mapping
  cat("\nWould you like to try mapping GenBank IDs to genes? (y/n): ")
  choice <- readline()
  
  if(tolower(choice) %in% c("y", "yes")) {
    map_genbank_to_genes(results$base_ids)
  }
  
  # Save results
  cat("\nSave analysis report? (y/n): ")
  save_choice <- readline()
  
  if(tolower(save_choice) %in% c("y", "yes")) {
    output_file <- paste0(tools::file_path_sans_ext(basename(file_path)), 
                         "_genbank_analysis.txt")
    
    sink(output_file)
    cat("GSE833 GenBank Accession Analysis Report\n")
    cat("=========================================\n\n")
    cat("File:", file_path, "\n")
    cat("Date:", date(), "\n\n")
    
    cat("Total annotations:", length(results$annotations), "\n")
    cat("Unique base IDs:", length(unique(results$base_ids)), "\n\n")
    
    cat("Pattern detection results:\n")
    print(results$pattern_results)
    cat("\n")
    
    cat("Common suffixes:\n")
    print(head(results$suffix_distribution, 10))
    cat("\n")
    
    cat("Common prefixes:\n")
    print(head(results$prefix_distribution, 10))
    cat("\n")
    
    cat("Sample annotations (first 50):\n")
    print(head(results$annotations, 50))
    cat("\n")
    
    cat("Sample base IDs (first 50):\n")
    print(head(unique(results$base_ids), 50))
    
    sink()
    cat("\nReport saved to:", output_file, "\n")
  }
  
  cat("\n=== IMPORTANT FOR GSE833 ===\n")
  cat("Your identifiers are GenBank accession numbers from older Affymetrix arrays.\n")
  cat("To map them to current gene symbols:\n")
  cat("1. Download GPL80 or GPL91 annotation from GEO\n")
  cat("2. Use: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL80\n")
  cat("3. Merge with your data using the GenBank accession (without suffix)\n")
  cat("\nExample R code for merging:\n")
  cat("# After downloading GPL annotation as gpl_data\n")
  cat("library(dplyr)\n")
  cat("your_data$base_id <- gsub('_.*$', '', your_data$column0)\n")
  cat("annotated_data <- left_join(your_data, gpl_data, by=c('base_id'='GB_ACC'))\n")
}

# Run if called from command line
if(!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  main()
} else {
  cat("To use this script:\n")
  cat("1. Source this file: source('gse833_analyzer.R')\n")
  cat("2. Run: main()\n")
  cat("3. Or call directly: check_gse833_annotation('your_file.tsv', 1)\n")
}