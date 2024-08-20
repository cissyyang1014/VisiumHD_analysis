# libraries
library(ggplot2)
library(patchwork)

# functions
plot_gene_VHD <- function(object, gene){
  df <- as.data.frame(cbind(object@meta.data[, c("x", "y")], object@assays$RNA@layers$data[which(object@assays$RNA@features@dimnames[[1]] == gene),]))
  colnames(df) <- c("x", "y", "feature")
  p <- ggplot(data = df, aes(x = x, y = y, color = feature)) + scale_color_viridis_c(option = "viridis") + 
    geom_point(size = 0.1) + ggtitle(gene) + theme_bw()
  return(p)
}

plot_gene_VHD_highlight <- function(object, gene){
  df <- as.data.frame(cbind(object@meta.data[, c("x", "y")], object@assays$RNA@layers$data[which(object@assays$RNA@features@dimnames[[1]] == gene),]))
  colnames(df) <- c("x", "y", "feature")
  df$highlight <- as.character(ifelse(df$feature == 0, 0, 1))
  p <- ggplot(data = df, aes(x = x, y = y, color = highlight)) +
    geom_point(size = 0.1) + ggtitle(gene) + theme_bw()+ scale_color_manual(values= c("grey50", "#FFFF00")) + theme(legend.position = "none")
  return(p)
}

plot_gene_VHD_highlight_2genes <- function(object, gene1, gene2){
  df <- as.data.frame(cbind(object@meta.data[, c("x", "y")], 
                            object@assays$RNA@layers$data[which(object@assays$RNA@features@dimnames[[1]] == gene1),], 
                            object@assays$RNA@layers$data[which(object@assays$RNA@features@dimnames[[1]] == gene2),]))
  colnames(df) <- c("x", "y", "feature1", "feature2")
  df$highlight1 <- as.character(ifelse(df$feature1 == 0, 0, 1))
  #df$size1 <- ifelse(df$highlight1 == 1, 5, 1)
  df$highlight2 <- as.character(ifelse(df$feature2 == 0, 0, 1))
  df$highlight <- as.character(ifelse(df$highlight1 == 1 & df$highlight2 == 1, 1, 0))
  
  p1 <- ggplot(data = df, aes(x = x, y = y, color = highlight1)) +
    geom_point(size = 0.01, alpha = 0.4) + ggtitle(gene1) + theme_void() + scale_color_manual(values= c("lightgrey", "#FFFF00")) + theme(legend.position = "none")
  p2 <- ggplot(data = df, aes(x = x, y = y, color = highlight2)) +
    geom_point(size = 0.01, alpha = 0.4) + ggtitle(gene2) + theme_void()+ scale_color_manual(values= c("lightgrey", "#4FC3F7")) + theme(legend.position = "none")
  p3 <- ggplot(data = df, aes(x = x, y = y, color = highlight)) +
    geom_point(size = 0.01, alpha = 0.4) + ggtitle("co-express") + theme_void()+ scale_color_manual(values= c("lightgrey", "red")) + theme(legend.position = "none")
  p <- p1 + p2 + p3 + plot_layout(ncol = 3, nrow = 1)
  return(p)
}
#plot_gene_VHD_highlight_2genes(object = object, gene1 = "POU4F2", "NEFM")

get_sparsity <- function(object, gene){
  df <- object@assays$RNA@layers$data[which(object@assays$RNA@features@dimnames[[1]] == gene),]
  sparsity <- sum(df == 0)/length(df)
  return(sparsity)
}
