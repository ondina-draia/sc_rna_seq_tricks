# Convert ENSEMBL ID to SYMBOL or other format of gene name inside a Seurat object

## Objective
Learn how to replace Ensembl IDs with Gene Symbols in Seurat objects, ensuring the proper handling of duplicates and NA values, and preparing the object for downstream analysis.

## Steps

### 1. Load required libraries
``` r
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(babelgene)
library(dplyr)
```


### 2. Step 2: Map Ensembl IDs to Gene Symbols

Extract Ensembl IDs from your Seurat object and map them to Gene Symbols using clusterProfiler

```
# Extract Ensembl IDs
human_genes <- rownames(human_main@assays$RNA)

# Convert Ensembl to Gene Symbols
human_gene_symbol <- bitr(human_genes, 
                          fromType = "ENSEMBL", 
                          toType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db)

# Count NA values
na_count <- sum(is.na(human_gene_symbol$SYMBOL))
print(paste("Number of NA values:", na_count))

```
### Step 3: Handle Duplicates and NAs

```
# Remove duplicates
human_gene_symbol_no_dup <- human_gene_symbol[!duplicated(human_gene_symbol$SYMBOL), ]

# Replace NA values with original Ensembl IDs
human_gene_symbol_no_dup$SYMBOL[is.na(human_gene_symbol_no_dup$SYMBOL)] <- 
  human_gene_symbol_no_dup$ENSEMBL

# Update the mapping table
human_gene_symbol <- human_gene_symbol_no_dup
```
### Step 4: Replace Ensembl IDs in Seurat Object

```
# Extract current Ensembl IDs
current_ensembl <- rownames(human_main@assays$RNA@counts)

# Match and replace with Gene Symbols
matched_symbols <- human_gene_symbol$SYMBOL[match(current_ensembl, human_gene_symbol$ENSEMBL)]

# Replace NAs with original Ensembl IDs
matched_symbols[is.na(matched_symbols)] <- current_ensembl

# Update rownames
rownames(human_main@assays$RNA@counts) <- matched_symbols
rownames(human_main@assays$RNA@data) <- matched_symbols
rownames(human_main@assays$RNA@scale.data) <- matched_symbols

```

### Step 5: Validate the Seurat Object

```
# Check for duplicates
dup_counts <- sum(duplicated(rownames(human_main@assays$RNA@counts)))
print(paste("Number of duplicate row names:", dup_counts))

# Check for NA values
na_counts <- sum(is.na(rownames(human_main@assays$RNA@counts)))
print(paste("Number of NA row names:", na_counts))

```

### Step 6: Reformat Meta Features

```
# Sync meta.features with normalized data
features <- rownames(human_main@assays$RNA@data)

meta.features <- data.frame(
  feature_name = features,
  feature_is_filtered = FALSE,
  feature_type = NA,
  row.names = features
)

human_main@assays$RNA@meta.features <- meta.features
```
### Step 7: Perform Downstream Analysis

```
human_main <- NormalizeData(human_main)
human_main <- ScaleData(human_main)
human_main <- RunPCA(human_main)
human_main <- FindNeighbors(human_main)
human_main <- FindClusters(human_main)
human_main <- RunUMAP(human_main, dims = 1:20)
```

Tadam!


