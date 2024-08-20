# analysis_visiumHD_8um.R
# Usage: $ Rscript analysis_visiumHD_8um.R spatial_8um.rds spatial_name output_PATH bulk_ref_name bulk_ref.csv nCount_RNA_cutoff

args = commandArgs(trailingOnly=TRUE)
#===============================================================================
# inputs:
spatial_PATH <- args[1]   # spatial data
spatial_name <- args[2]   # spatial sample name (only for plotting)
output_folder <- args[3]  # output folder
nCount_cutoff <- args[6]  # cutoff of nCount_RNA used to filter spots
bulk_name <- args[4]      # sample name of the bulk RNA reference (only for plotting)
bulk_PATH <- args[5]      # bulk RNA reference (built from scRNA-seq) 
#===============================================================================
# libraries and functions
library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

# process:
spatial_obj <- readRDS(spatial_PATH)
#str(spatial_obj)
assay_nms <- Assays(spatial_obj)
mtx <- spatial_obj@assays[[assay_nms]]@layers$counts
colnames(mtx) <- colnames(spatial_obj)
rownames(mtx) <- rownames(spatial_obj)
object <- CreateSeuratObject(counts = mtx, min.cells = 5)
object@meta.data$log10_nCount_RNA <- log10(object@meta.data$nCount_RNA)
object@meta.data$log10_nFeature_RNA <- log10(object@meta.data$nFeature_RNA)
nms <- names(spatial_obj@images)
object@meta.data$x <- spatial_obj@images[[nms]]@boundaries$centroids@coords[, 1]
object@meta.data$y <- spatial_obj@images[[nms]]@boundaries$centroids@coords[, 2]

# plot original spatial w/ nCount_RNA
pdf(paste0(output_folder, "/original_log10_nCount_Spatial.pdf"))
p <- ggplot(data = object@meta.data, aes(x = x, y = y, color = log10_nCount_RNA)) +
  geom_point(size = 0.5) + theme_bw() + ggtitle("Original nCount_RNA") + labs(subtitle = spatial_name) + scale_color_viridis_c(option = "viridis")
print(p)
dev.off()

pdf(paste0(output_folder, "/original_log10_nCount_hist.pdf"))
p <- ggplot(object@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nCount_RNA - original", x = "nCount_RNA", y = "Frequency") +
  theme_classic() + annotate("text", x = mean(object@meta.data$nCount_RNA), y = Inf, 
                             label = paste0("mean=", round(mean(object@meta.data$nCount_RNA), 3)), 
                             hjust = 1, vjust = 1, color = "red")
#print(p)
p_box <- ggplot(object@meta.data, aes(x="", y = nCount_RNA)) +
  geom_boxplot(fill = "blue", color = "black") + 
  coord_flip() +
  theme_classic() +
  xlab("") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(p + p_box + plot_layout(nrow = 2, heights = c(2, 1)))

p <- ggplot(object@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nCount_RNA - original", x = "nCount_RNA", y = "Frequency") +
  theme_bw() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept = nCount_cutoff, linetype = "dashed", color = "red")
print(p)
dev.off()

pdf(paste0(output_folder, "/original_log10_nFeature_hist.pdf"))
p <- ggplot(object@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nFeature_RNA - original", x = "nFeature_RNA", y = "Frequency") +
  theme_bw() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
print(p)
p <- ggplot(object@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nFeature_RNA - original", x = "nFeature_RNA", y = "Frequency") +
  theme_classic() + annotate("text", x = mean(object@meta.data$nFeature_RNA), y = Inf, 
                             label = paste0("mean=", round(mean(object@meta.data$nFeature_RNA), 3)), 
                             hjust = 1, vjust = 1, color = "red")
p_box <- ggplot(object@meta.data, aes(x="", y = nFeature_RNA)) +
  geom_boxplot(fill = "blue", color = "black") + 
  coord_flip() +
  theme_classic() +
  xlab("") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(p + p_box + plot_layout(nrow = 2, heights = c(2, 1)))
dev.off()

# subset by nCount_RNA
object <- subset(object, subset = nCount_RNA > nCount_cutoff)

pdf(paste0(output_folder, "/subset_log10_nCount_Spatial.pdf"))
p <- ggplot(data = object@meta.data, aes(x = x, y = y, color = log10_nCount_RNA)) +
  geom_point(size = 0.5) + theme_bw() + ggtitle("Subset nCount_RNA") + labs(subtitle = spatial_name) + scale_color_viridis_c(option = "viridis")
print(p)
dev.off()

pdf(paste0(output_folder, "/subset_log10_nCount_hist.pdf"))
p <- ggplot(object@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nCount_RNA - subset", x = "nCount_RNA", y = "Frequency") +
  theme_bw()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 
print(p)

p <- ggplot(object@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nCount_RNA - subset", x = "nCount_RNA", y = "Frequency") +
  theme_classic() + annotate("text", x = mean(object@meta.data$nCount_RNA), y = Inf, 
                             label = paste0("mean=", round(mean(object@meta.data$nCount_RNA), 3)), 
                             hjust = 1, vjust = 1, color = "red")
#print(p)
p_box <- ggplot(object@meta.data, aes(x="", y = nCount_RNA)) +
  geom_boxplot(fill = "blue", color = "black") + 
  coord_flip() +
  theme_classic() +
  xlab("") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(p + p_box + plot_layout(nrow = 2, heights = c(2, 1)))
dev.off()

pdf(paste0(output_folder, "/subset_log10_nFeature_hist.pdf"))
p <- ggplot(object@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nFeature_RNA - subset", x = "nFeature_RNA", y = "Frequency") +
  theme_classic() + annotate("text", x = mean(object@meta.data$nFeature_RNA), y = Inf, 
                             label = paste0("mean=", round(mean(object@meta.data$nFeature_RNA), 3)), 
                             hjust = 1, vjust = 1, color = "red")
p_box <- ggplot(object@meta.data, aes(x="", y = nFeature_RNA)) +
  geom_boxplot(fill = "blue", color = "black") + 
  coord_flip() +
  theme_classic() +
  xlab("") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(p + p_box + plot_layout(nrow = 2, heights = c(2, 1)))

p <- ggplot(object@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of nFeature_RNA - subset", x = "nFeature_RNA", y = "Frequency") +
  theme_bw() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
print(p)
dev.off()


# transcripts level
bulk_data <- read.csv(bulk_PATH)[, -1]
colnames(bulk_data) <- c("gene", "sc_bulk")
bulk_data$sc_bulk <- as.numeric(bulk_data$sc_bulk)

Spatial_trancripts <- as.data.frame(rowSums(object@assays$RNA$counts))
colnames(Spatial_trancripts) <- "spatial_counts"
Spatial_trancripts$gene <- rownames(Spatial_trancripts)
Spatial_trancripts$spatial_counts <- as.numeric(Spatial_trancripts$spatial_counts)

Spatial_trancripts <- Spatial_trancripts %>% left_join(bulk_data)
Spatial_trancripts <- Spatial_trancripts[Spatial_trancripts$spatial_counts != 0,]
Spatial_trancripts <- Spatial_trancripts[Spatial_trancripts$sc_bulk != 0,]
Spatial_trancripts$log_spatial_counts <- log10(Spatial_trancripts$spatial_counts)
Spatial_trancripts$log_sc_bulk <- log10(Spatial_trancripts$sc_bulk)

model <- lm(log_spatial_counts ~ log_sc_bulk, data = Spatial_trancripts)
model_summary <- summary(model)
r_squared <- model_summary$r.squared
r_squared
pcc <- cor(Spatial_trancripts$log_spatial_counts, Spatial_trancripts$log_sc_bulk)

pdf(paste0(output_folder, "/transcripts_level.pdf"))
p <- ggplot(data = Spatial_trancripts) +
  geom_point(aes(x = sc_bulk, y = spatial_counts)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()+ annotation_logticks()  +
  ggtitle("Transcripts Level") +
  labs(subtitle = paste0(spatial_name, " vs. ", bulk_name)) + xlab("bulk RNA") + ylab("Spatial RNA") +
  geom_abline(slope=1, intercept = 0, col = "red", size = 1.5) + 
  annotate("text", x = min(Spatial_trancripts$sc_bulk), y = max(Spatial_trancripts$spatial_counts), 
           label = (as.expression(bquote(R^2 == .(round(r_squared, 3))))), 
           hjust = 0, vjust = 1, color = "red") + 
  annotate("text", x = min(Spatial_trancripts$sc_bulk), y = (max(Spatial_trancripts$spatial_counts) - 200000),
           label = paste0("PCC = ", round(pcc, 3)), 
           hjust = 0, vjust = 1, color = "red")
print(p)
dev.off()


#===============================================================================
# normalization and clustering
object <- NormalizeData(object)

object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(object) <- "sketch"

object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 0.05)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

object <- ProjectData(
  object = object,
  assay = "RNA",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# plot umaps
DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")
# switch to full dataset
DefaultAssay(object) <- "RNA"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

pdf(paste0(output_folder, "/clustering_umaps.pdf"), width = 10)
print(p1 | p2)
dev.off()

pdf(paste0(output_folder, "/seurat_clusters.pdf"))
p <- ggplot(data = object@meta.data, aes(x = x, y = y, color = seurat_cluster.projected)) +
  geom_point(size = 0.5) + theme_bw() + ggtitle("Seurat clusters") + labs(subtitle = spatial_name)
print(p)
dev.off()

pdf(paste0(output_folder, "/umaps_nCount_nFeature.pdf"))
FeaturePlot(object, features = "nCount_RNA")
FeaturePlot(object, features = "log10_nCount_RNA")
FeaturePlot(object, features = "log10_nFeature_RNA")
FeaturePlot(object, features = "nFeature_RNA")
dev.off()


# plot cluster sperately
df <- as.data.frame(cbind(object@meta.data[, c("x", "y", "seurat_cluster.projected")]))
pdf(paste0(output_folder, "/seurat_clusters_sep.pdf"))
for (i in c(0:(length(unique(object@meta.data$seurat_cluster.projected))-1))){
  df$highlight <- as.character(ifelse(df$seurat_cluster.projected == i, 1, 0))
  p <- ggplot(data = df, aes(x = x, y = y, color = highlight)) +
    geom_point(size = 0.1) + ggtitle(paste0("Cluster ", i)) + theme_bw()+ scale_color_manual(values= c("grey50", "#FFFF00")) + theme(legend.position = "none")
  print(p)
}
dev.off()

# find markers
# Crete downsampled object to make visualization either
DefaultAssay(object) <- "RNA"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["RNA"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "RNA"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "RNA", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "RNA", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
write.csv(markers, paste0(output_folder, "/makers.csv"))

object_subset <- ScaleData(object_subset, assay = "RNA", features = top5$gene)
pdf(paste0(output_folder, "/heatmap_seurat_clusters.pdf"))
p <- DoHeatmap(object_subset, assay = "RNA", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
print(p)
dev.off()

# save object
saveRDS(object, paste0(output_folder, "/object.rds"))
