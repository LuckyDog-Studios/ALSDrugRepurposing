# Set CRAN mirror first
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Then install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery", "Biobase"))