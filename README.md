# VisiumHD Analysis
Version 1.0 (Jul 2024)

This repository is used to analyze the VisiumHD data in a basic way.

## Pipeline for Basic Analysis Based on Bin-size

### Load the outputs from spaceranger and split by bin-size
```
$ mkdir output_PATH
$ Rscript load_split_visiumHD.R sample_name spaceranger_outs_PATH output_PATH
```

### Pre-process analysis based on the 8um bin-size
```
$ Rscript pp_visiumHD_8um.R spatial_8um.rds spatial_name output_PATH bulk_ref_name bulk_ref.csv
```

### Run basic analysis based on the 8um bin-size
```
$ Rscript analysis_visiumHD_8um.R spatial_8um.rds spatial_name output_PATH bulk_ref_name bulk_ref.csv nCount_RNA_cutoff cluster_resolution
```

## Other Analysis

### Visualization - gene expression and co-expression
```
# R
source("functions_visiumHD.R")
# 1) To plot the gene expression on the spatial
plot_gene_VHD(object, gene)
# 2) To highlight the spots with the gene expressed on the spatial
plot_gene_VHD_highlight(object, gene)
# 3) To highlight the spots with the gene1 expressed, gene2 expressed, and both genes co-expressed on the spatial separately
plot_gene_VHD_highlight_2genes(object, gene1, gene2)
```

### Deconvolution using spacexr
Please refer to the example code `visiumHD_spacexr.R`. - under construction
