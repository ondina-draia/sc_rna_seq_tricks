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

# Update the mapping table
human_gene_symbol <- human_gene_symbol_no_dup
```
### Step 4: Replace Ensembl IDs in Seurat Object

```
# Extract current Ensembl IDs
current_ensembl <- rownames(human_main@assays$RNA@counts)

# Match and replace with Gene Symbols
matched_symbols <- human_gene_symbol$SYMBOL[match(current_ensembl, human_gene_symbol$ENSEMBL)]

# Update rownames
rownames(human_main@assays$RNA@counts) <- matched_symbols
rownames(human_main@assays$RNA@data) <- matched_symbols
rownames(human_main@assays$RNA@scale.data) <- matched_symbols

```
### Step 5: Update Gene Symbols Across Seurat Data Slots

```
update_gene_names <- function(seurat_object, gene_mapping) {
  # Update counts and data
  for (slot in c("counts", "data")) {
    current_ensembl <- rownames(seurat_object@assays$RNA[[slot]])
    matched_symbols <- gene_mapping$SYMBOL[match(current_ensembl, gene_mapping$ENSEMBL)]
    
    # Replace only valid mappings, keep unmatched as original
    matched_symbols[is.na(matched_symbols)] <- current_ensembl
    
    rownames(seurat_object@assays$RNA[[slot]]) <- matched_symbols
  }

  # Update scale.data (only matching genes)
  scale_genes <- rownames(seurat_object@assays$RNA@scale.data)
  matched_scale_symbols <- gene_mapping$SYMBOL[match(scale_genes, gene_mapping$ENSEMBL)]
  
  # Keep original names for unmatched genes
  matched_scale_symbols[is.na(matched_scale_symbols)] <- scale_genes
  rownames(seurat_object@assays$RNA@scale.data) <- matched_scale_symbols

  # Update variable features
  current_var_features <- seurat_object@assays$RNA@var.features
  matched_var_symbols <- gene_mapping$SYMBOL[match(current_var_features, gene_mapping$ENSEMBL)]
  matched_var_symbols[is.na(matched_var_symbols)] <- current_var_features
  seurat_object@assays$RNA@var.features <- matched_var_symbols

  # Check for duplicates or NAs
  for (slot in c("counts", "data", "scale.data")) {
    duplicate_count <- sum(duplicated(rownames(seurat_object@assays$RNA[[slot]])))
    print(paste("Duplicates in", slot, ":", duplicate_count))
    
    na_count <- sum(is.na(rownames(seurat_object@assays$RNA[[slot]])))
    print(paste("NA values in", slot, ":", na_count))
  }

  return(seurat_object)
}

# Apply the function
human_main <- update_gene_names(human_main, human_gene_symbol)

```

### Step 6: Validate the Seurat Object

```
all(rownames(human_main@assays$RNA@counts) == 
    rownames(human_main@assays$RNA@data)) &&
all(rownames(human_main@assays$RNA@scale.data) %in% 
    rownames(human_main@assays$RNA@counts))

# Check for duplicates
dup_counts <- sum(duplicated(rownames(human_main@assays$RNA@counts)))
print(paste("Number of duplicate row names:", dup_counts))

# Check for NA values
na_counts <- sum(is.na(rownames(human_main@assays$RNA@counts)))
print(paste("Number of NA row names:", na_counts))

```

### Step 7: Reformat Meta Features

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
### Step 8: Perform Downstream Analysis

```
human_main <- NormalizeData(human_main)
human_main <- ScaleData(human_main)
human_main <- RunPCA(human_main)
human_main <- FindNeighbors(human_main)
human_main <- FindClusters(human_main)
human_main <- RunUMAP(human_main, dims = 1:20)
```

Tadam!


