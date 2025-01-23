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


