#!/usr/bin/env Rscript
# Master DE Merger: Combining 6 ALS datasets without filtering

library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(stringr)

# 1. Setup paths - Change "differential_expression/" to your actual folder
input_folder <- "differential_expression/"
file_list <- list.files(path = input_folder, pattern = "*.csv", full.names = TRUE)

# 2. Processing function for each individual file
process_study <- function(file_path) {
  # Get a clean name for the study from the filename
  study_name <- tools::file_path_sans_ext(basename(file_path))
  
  df <- read_csv(file_path, show_col_types = FALSE) %>%
    filter(!is.na(Gene_Symbol)) %>%
    # Select key columns
    select(Gene_Symbol, logFC, adj.P.Val) %>%
    # Handle duplicates by keeping the most significant probe/entry
    group_by(Gene_Symbol) %>%
    slice_min(order_by = adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Create the Direction/Status column
    mutate(
      Status = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "UP",
        adj.P.Val < 0.05 & logFC < 0 ~ "DOWN",
        TRUE ~ "NS" # Non-Significant
      )
    ) %>%
    # Rename columns to include the study name
    rename_with(~paste0(study_name, "_", .), -Gene_Symbol)
  
  return(df)
}

# 3. Merge all files (Full Join = No Filtering)
cat("Merging all datasets into one Master Table...\n")
master_table <- file_list %>% 
  map(process_study) %>% 
  reduce(full_join, by = "Gene_Symbol")

# 4. Create Summary Columns
# This helps you see at a glance how many datasets a gene appears in
cat("Calculating consensus metrics...\n")

# Identify which columns are "Status" columns
status_cols <- grep("_Status$", colnames(master_table), value = TRUE)

master_table <- master_table %>%
  mutate(
    # Count how many datasets the gene actually appeared in
    Total_Datasets_Present = rowSums(!is.na(select(., all_of(status_cols)))),
    
    # Count how many times it was UP or DOWN across those datasets
    Count_UP = rowSums(select(., all_of(status_cols)) == "UP", na.rm = TRUE),
    Count_DOWN = rowSums(select(., all_of(status_cols)) == "DOWN", na.rm = TRUE),
    
    # List the names of the datasets where it was significant (UP or DOWN)
    Significant_In = apply(select(., all_of(status_cols)), 1, function(x) {
      names <- names(x)[x %in% c("UP", "DOWN")]
      # Clean the names to just show the study title
      paste(gsub("_Status", "", names), collapse = ", ")
    })
  ) %>%
  # Move summary columns to the front for easier reading
  select(Gene_Symbol, Total_Datasets_Present, Count_UP, Count_DOWN, Significant_In, everything())

# 5. Save the full Master List
write_csv(master_table, "differential_expression/consensus/ALS_Master_Consensus_Unfiltered.csv")
cat("Master Table saved! Total unique genes found:", nrow(master_table), "\n")