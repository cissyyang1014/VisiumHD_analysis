# Load the outputs from spaceranger and split into different bin objects
# Usage: $ Rscript load_split_visiumHD.R sample_name spaceranger_outs_PATH output_PATH

args <-  commandArgs(trailingOnly=TRUE)

library(Seurat)
sample_name <- args[1]
localdir <- args[2]
out_dir <- args[3]
for (i in c(2, 8, 16)){
  object <- Load10X_Spatial(data.dir = localdir, bin.size = i)
  saveRDS(object, paste0(out_dir, "/", sample_name, "_", i, "um.rds"))
}

