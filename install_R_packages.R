# Installs the R packages used by the SqnPAN scripts (R 4.6.0 was used; versions in README.md).
cran <- c("survival", "ggplot2", "patchwork", "ggrepel", "scales", "ragg", "systemfonts", "readxl", "ActiveDriver", "BiocManager", "remotes")
miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
if (!requireNamespace("MutationalPatterns", quietly = TRUE)) BiocManager::install("MutationalPatterns")
if (!requireNamespace("dndscv", quietly = TRUE)) remotes::install_github("im3sanger/dndscv")
