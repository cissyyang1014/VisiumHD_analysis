# pp_visiumHD_8um.R
# Usage: $ Rscript pp_visiumHD_8um.R spatial_8um.rds spatial_name output_PATH bulk_ref_name bulk_ref.csv

args = commandArgs(trailingOnly=TRUE)
#===============================================================================
# inputs:
spatial_PATH <- args[1]   # spatial data
spatial_name <- args[2]   # spatial sample name (only for plotting)
output_folder <- args[3]  # output folder
nCount_cutoff <- 30       # cutoff of nCount_RNA used to filter spots, default = 30
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


# mkdir output folder
if (!file.exists(output_folder)){
  dir.create(output_folder)
}

# process:
spatial_obj <- readRDS(spatial_PATH)
#str(spatial_obj)
assay_nms <- Assays(spatial_obj)
mtx <- spatial_obj@assays[[assay_nms]]@layers$counts
colnames(mtx) <- colnames(spatial_obj)
rownames(mtx) <- rownames(spatial_obj)
object <- CreateSeuratObject(mtx)
object@meta.data$log10_nCount_RNA <- log10(object@meta.data$nCount_RNA)
object@meta.data$log10_nFeature_RNA <- log10(object@meta.data$nFeature_RNA)
nms <- names(spatial_obj@images)
object@meta.data$x <- spatial_obj@images[[nms]]@boundaries$centroids@coords[, 1]
object@meta.data$y <- spatial_obj@images[[nms]]@boundaries$centroids@coords[, 2]

saveRDS(object, paste0(output_folder, "/original_object.rds"))

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


# library complexity for cutoff selection
cutoffs <- seq(0, 200, 10)
n_feature_1 <- as.numeric()
n_feature_3 <- as.numeric()
n_feature_5 <- as.numeric()
n_spots <- as.numeric()
for (i in c(1:length(cutoffs))){
  cutoff <- cutoffs[i]
  ncounts <- colSums(mtx)
  if (cutoff != 0){
    rm_index <- which(ncounts < cutoff)
    #dim(mtx)
    mtx <- mtx[, -rm_index]
  }
  n_spots[i] <- dim(mtx)[2]

  n_feature_1[i] <- sum(rowSums(mtx > 0) > 0)
  n_feature_3[i] <- sum(rowSums(mtx > 0) > 2)
  n_feature_5[i] <- sum(rowSums(mtx > 0) > 4)
  
}

pdf(paste0(output_folder, "/lib_complexity.pdf"))
compx_df <- as.data.frame(cbind(cutoffs, n_spots, n_feature_1, n_feature_3, n_feature_5))
df_long <- pivot_longer(compx_df, cols = c(n_feature_1, n_feature_3, n_feature_5), names_to = "variable", values_to = "value")
p <- ggplot(data = df_long, aes(x = cutoffs, y = value, color = variable)) +
  geom_line() +
  theme_bw() +
  labs(title = "Library Complexity", x = "Cutoff of nCount_RNA", y = "Number of Features Detected", color = "", subtitle = spatial_name) + 
  scale_color_manual(values = c("n_feature_1" = "pink", "n_feature_3" = "darkgreen", "n_feature_5" = "darkblue"),
                     labels = c("n_feature_1" = "min.cells=1", "n_feature_3" = "min.cells=3", "n_feature_5" = "min.cells=5"))
print(p)

p <- ggplot(data = compx_df, aes(x = cutoffs, y = n_spots)) +
  geom_line() +
  theme_bw() +
  labs(title = "Number of Spots", x = "Cutoff of nCount_RNA", y = "Number of Spots", subtitle = spatial_name)
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


# subset by nCount_RNA - by default cutoff = 30
object <- subset(object, subset = nCount_RNA > nCount_cutoff)

pdf(paste0(output_folder, "/subset_log10_nCount_Spatial_default30.pdf"))
p <- ggplot(data = object@meta.data, aes(x = x, y = y, color = log10_nCount_RNA)) +
  geom_point(size = 0.5) + theme_bw() + ggtitle("Subset nCount_RNA") + labs(subtitle = spatial_name) + scale_color_viridis_c(option = "viridis")
print(p)
dev.off()

pdf(paste0(output_folder, "/subset_log10_nCount_hist_default30.pdf"))
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

pdf(paste0(output_folder, "/subset_log10_nFeature_hist_default30.pdf"))
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

