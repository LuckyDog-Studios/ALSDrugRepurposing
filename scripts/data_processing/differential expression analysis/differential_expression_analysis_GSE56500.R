library(AnnotationDbi)
library(limma)
library(org.Hs.eg.db)
library(dplyr)
library(GEOquery)

# Load series matrix
expr_data <- read.table("../../datasets/raw/GSE56500/GSE56500_series_matrix.txt",
                        header=TRUE, sep="\t", quote="", comment.char="!",
                        row.names=1)

# Create the sample group vector based on your matrix
sample_group <- c(
  "control", "control", "C9ORF72",
  "control", "control", "C9ORF72",
  "sALS", "sALS", "control",
  "control", "C9ORF72", "sALS"
)

# Convert to factor
group <- factor(sample_group, levels = c("control", "sALS", "C9ORF72"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Define contrasts
contrast_matrix <- makeContrasts(
  sALS_vs_control = sALS - control,
  C9ORF72_vs_control = C9ORF72 - control,
  levels = design
)

fit <- lmFit(expr_data, design)

fit2 <- contrasts.fit(fit, contrast_matrix)  # apply contrast
fit2 <- eBayes(fit2)                          # empirical Bayes smoothing

res_sALS <- topTable(fit2, coef="sALS_vs_control", adjust.method="BH", number=Inf)
res_C9   <- topTable(fit2, coef="C9ORF72_vs_control", adjust.method="BH", number=Inf)

# Download annotation table for GPL5188
gpl <- getGEO("GPL5188", AnnotGPL = TRUE)
annot <- Table(gpl)

map_df <- annot[, c("ID", "gene_assignment")]
colnames(map_df) <- c("ProbeID", "GeneSymbol")

# Optional: clean up multiple gene assignments
map_df$GeneSymbol <- sapply(strsplit(map_df$GeneSymbol, " // "), function(x) {
  if (length(x) > 1) x[2] else NA
})

# Add ProbeID column
res_sALS$ProbeID <- rownames(res_sALS)
res_C9$ProbeID   <- rownames(res_C9)

# Merge with annotation to create res_sALS_annot and res_C9_annot
res_sALS_annot <- merge(res_sALS, map_df, by = "ProbeID", all.x = TRUE)
res_C9_annot   <- merge(res_C9, map_df, by = "ProbeID", all.x = TRUE)

# Map to full gene names
res_sALS_annot$GeneName <- mapIds(
  org.Hs.eg.db,
  keys = unique(na.omit(res_sALS_annot$GeneSymbol)),
  column = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first"
)[res_sALS_annot$GeneSymbol]

res_C9_annot$GeneName <- mapIds(
  org.Hs.eg.db,
  keys = unique(na.omit(res_C9_annot$GeneSymbol)),
  column = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first"
)[res_C9_annot$GeneSymbol]

write.csv(res_sALS_annot, "../../datasets/processed/GSE56500/results_sALS_genes.csv", row.names = FALSE)
write.csv(res_C9_annot, "../../datasets/processed/GSE56500/results_C9_genes.csv", row.names = FALSE)

cat("Wrote to csvs")