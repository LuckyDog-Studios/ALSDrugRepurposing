library(DESeq2)

# 1. Read your data
counts <- read.csv("../../datasets/raw/GSE67196/GSE67196_Petrucelli2015_ALS_genes.rawcount.txt", sep="\t", header = TRUE)
rownames(counts) <- make.unique(counts$GeneID)
counts <- counts[, -(1:2)]  # remove GeneID and Chr columns after setting rownames

# 2. Make sample info
condition <- factor(c(rep("fcx", 27), rep("cereb", 26)))

coldata <- data.frame(row.names=colnames(counts), condition)

# 3. Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# 4. Run DESeq2
dds <- DESeq(dds)

# 5. Get results
res <- results(dds)

res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

write.csv(res_df,
          file = "../../datasets/processed/GSE67196/results_GSE67196_cereb_vs_fc.csv",
          row.names = FALSE)

