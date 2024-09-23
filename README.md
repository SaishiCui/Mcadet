# Mcadet 

Mcadet is a feature selection framework for scRNA-seq data (count-based) (UMI-like or full length). Mcadet integrates Multiple Correspondence Analysis (MCA), graph-based community detection, and a novel statistical testing approach.

![Figure1](https://github.com/user-attachments/assets/ffc7e9ab-459d-4680-8dbc-5697b9ca3b4f)

## Usage:

```r
# Step 1: install devtools package (if not installed)
install.packages("devtools")

# Step 2: install Mcadet from Github
devtools::install_github("SaishiCui/Mcadet")

# Step 3: load Mcadet
library(Mcadet)

# Step 4: install required package (if not installed)
install.packages(c("cccd", "igraph", "irlba"))     
BiocManager::install("Seurat")

# Example function call:
mcadet(data, n.comp=60,  run=10, n.feature=NA, nk_percent = 0.005, start_resolution = 0.5,
                 cell_percent=0.005, MC_iter = 50000, fdr = 0.15, seed = 1234)
```

# Function Parameters:

1) "data" is the input data, a raw count matrix. Genes are on rows and cells on columns. The matrix should have row names that represent gene names.

2) "n.comp" is the number of PCs selected from MCA decomposition, default is 60.

3) "run" is the number of iterations you want to run to get a more robust result, default is 10.

4) "n.feature" is the number of genes to select, default is NA, which means the proposed statistical testing method will be used.

5) "start_resolution" is the starting resolution parameter for Leiden algorithm, increasing by 0.1 for each run.

6) "cell_percent": If a gene is only expressed in "cell_percent" (i.e., 0.5%) of cells, we preclude these noninformative genes. Default is 0.005.

7) "MC_iter": Monte-Carlo iteration times, default is 50,000.

8) "fdr": False discovery rate for the proposed statistical testing rate.

9) "seed": Random seed for Leiden algorithm.
