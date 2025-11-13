library(limma)
library(dplyr)

# Read the file
series_data <- readLines("./datasets/raw/GSE112676/GSE112676_series_matrix.txt")

# Extract the data matrix
data_start <- grep("^!series_matrix_table_begin", series_data) + 1
data_end <- grep("^!series_matrix_table_end", series_data) - 1
data_lines <- series_data[data_start:data_end]
expression_matrix <- read.table(text = data_lines, header = TRUE, sep = "\t",
                               row.names = 1, check.names = FALSE, comment.char = "!")
exprs <- as.matrix(expression_matrix)

# Extract diagnosis information CORRECTLY
diagnosis_lines <- grep("^!Sample_characteristics_ch1", series_data, value = TRUE)

# Parse diagnosis - look for lines containing "diagnosis:"
diagnosis <- character()
for (line in diagnosis_lines) {
  # Split by tab to get individual entries for each sample
  entries <- strsplit(line, "\t")[[1]]
  for (entry in entries) {
    if (grepl("diagnosis:", entry)) {
      # Extract the diagnosis value
      diag_value <- gsub('.*diagnosis:\\s*([^"]+).*', '\\1', entry)
      diagnosis <- c(diagnosis, diag_value)
    }
  }
}

# Create sample data
sample_data <- data.frame(
  Sample = colnames(exprs),
  Diagnosis = diagnosis,
  stringsAsFactors = FALSE
)

# Check what we found
print("Diagnosis distribution:")
print(table(sample_data$Diagnosis))
print(paste("Number of samples:", nrow(sample_data)))

# Clean the diagnosis data
sample_data$Diagnosis <- trimws(sample_data$Diagnosis)  # Remove whitespace
sample_data$Diagnosis <- gsub("^diagnosis:\\s*", "", sample_data$Diagnosis)  # Remove "diagnosis:" prefix

# Convert to factor and check
sample_data$Diagnosis <- factor(sample_data$Diagnosis)
print("Final diagnosis distribution:")
print(table(sample_data$Diagnosis))
print("Levels:")
print(levels(sample_data$Diagnosis))

# Only proceed if we have both ALS and CON
if (length(levels(sample_data$Diagnosis)) >= 2) {
  # Create design matrix (without sex)
  design <- model.matrix(~ 0 + Diagnosis, data = sample_data)
  colnames(design) <- levels(sample_data$Diagnosis)

  print("Design matrix preview:")
  print(head(design))

  # Fit linear model
  fit <- lmFit(exprs, design)

  # Create contrast (ALS vs CON)
  contrast_matrix <- makeContrasts(ALS_vs_CON = ALS - CON, levels = design)

  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # Get results
  results <- topTable(fit2, coef = "ALS_vs_CON", number = Inf, adjust.method = "BH")

  # Add significance flags
  results$significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, TRUE, FALSE)

  # View top results
  print("Top 20 differentially expressed genes:")
  print(head(results, 20))

  # Summary
  print("Summary of results:")
  print(summary(decideTests(fit2, p.value = 0.05, lfc = 0.5)))

  # Save results
  write.csv(results, "./datasets/processed/GSE112676/results_GSE112676.csv", row.names = TRUE)

  # Create visualizations
  volcanoplot(fit2, coef = "ALS_vs_CON", highlight = 10,
              main = "ALS vs Control - Volcano Plot")

  plotMD(fit2, column = "ALS_vs_CON", status = decideTests(fit2)[, "ALS_vs_CON"],
         main = "ALS vs Control - MA Plot")

} else {
  print("ERROR: Need both ALS and CON groups for comparison.")
  print("Current diagnosis levels:")
  print(levels(sample_data$Diagnosis))
}