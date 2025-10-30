library(tidyverse)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

dat <- read.delim("/Users/cmdb/qb25-answers/week7/read_matrix.tsv", header = TRUE, check.names = FALSE)
#Chatgpt suggested that check.names creates a cleaner loading pattern
rownames(dat) <- paste0("gene_", seq_len(nrow(dat)))
mat <- as.matrix(dat)
storage.mode(mat) <- "double"
gene_sd <- rowSds(mat) #top 500 most variant genes
top_idx <- order(gene_sd, decreasing = TRUE)[1:500]
mat_top <- mat[top_idx, ]
#prcomp() does principal component analysis but collecting samples in rows
pca <- prcomp(t(mat_top), scale. = TRUE, center = TRUE)
var_expl <- pca$sdev^2
var_expl <- var_expl / sum(var_expl)
#Chatgpt suggested using tidyr::seperate() to extract tissue and replicate info from sample names
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$sample <- rownames(pca$x)
pca_df <- separate(pca_df, sample, into = c("tissue","replicate"), sep = "_", remove = FALSE)
pca_df$replicate <- gsub("^Rep", "", pca_df$replicate)
#Exercise 1 visualization
# PCA scatter
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue, shape = replicate)) +
  geom_point(size = 3) +
  labs(title = "PCA of top 500 variable genes",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_classic()
#Chatgpt helped me create a cleaner structure for the ggplot call

# Scree 
scree_df <- data.frame(PC = seq_along(var_expl), VarExpl = var_expl)
p_scree <- ggplot(scree_df, aes(x = PC, y = VarExpl)) +
  geom_col() +
  labs(title = "% Variance Explained",
       x = "Principal Component",
       y = "Proportion of Variance") +
  theme_classic()

outdir <- "/Users/cmdb/qb25-answers/week7"
dir.create(outdir, showWarnings = FALSE)

ggsave(file.path(outdir, "pca_plot_week6.png"),   plot = p_pca,   width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "scree_plot_week6.png"), plot = p_scree, width = 7, height = 5, dpi = 300)
#Chatgpt showed me how to save the plots safely to the week 7 folder. "outdir" is a variable that encodes holding the folder path

#Exercise 2 script
samples <- colnames(mat)
meta <- tibble(sample = samples) |>
  tidyr::separate(sample, into = c("tissue","replicate"), sep = "_", remove = FALSE)
meta$replicate <- gsub("^Rep", "", meta$replicate)
meta$replicate <- as.integer(meta$replicate)
#tissue replicates are ordered together
o <- order(meta$tissue, meta$replicate)
meta <- meta[o, ]
mat  <- mat[, o, drop = FALSE]
if (ncol(mat) == 21 && all(table(meta$tissue) == 3)) {
  combined <- mat[, seq(1, 21, 3), drop = FALSE]
  combined <- combined + mat[, seq(2, 21, 3), drop = FALSE]
  combined <- combined + mat[, seq(3, 21, 3), drop = FALSE]
  combined <- combined / 3
  colnames(combined) <- meta$tissue[seq(1, 21, 3)]
#I was suggested by Chatgpt to create a script for fallback to average tissue labels instead of fixed positions
} else { 
  tissues <- unique(meta$tissue)
  combined <- sapply(tissues, function(tt) {
    cols <- meta$sample[meta$tissue == tt]
    rowMeans(mat[, cols, drop = FALSE])
  })
  combined <- as.matrix(combined)
  colnames(combined) <- tissues
}
g_sd <- matrixStats::rowSds(combined)
filt <- combined[g_sd > 1, , drop = FALSE]

#k-means clustering, sed.seed(42) fixes randomness and nstart=100 suns k-means multiple times with different starts
set.seed(42)
km <- kmeans(filt, centers = 12, nstart = 100)
labels <- km$cluster
#order rows
ord_rows  <- order(labels)
filt_ord  <- filt[ord_rows, , drop = FALSE]
labels_ord <- labels[ord_rows]
#heatmap
hm_path <- file.path(outdir, "heatmap_week6.png")
png(hm_path, width = 1800, height = 2200, res = 200)
heatmap(filt_ord,
        Rowv = NA, Colv = NA,
        # Chatgpt mentioned `scale="row"` makes patterns clearer;
        scale = "row",
        RowSideColors = RColorBrewer::brewer.pal(12, "Paired")[labels_ord],
        labRow = NA,
        margins = c(8, 6),
        xlab = "Tissues (replicate-averaged)",
        ylab = "Genes (SD > 1)")
legend("right", inset = c(-0.2, 0), xpd = TRUE, bty = "n",
       fill = RColorBrewer::brewer.pal(12, "Paired"),
       legend = paste("Cluster", 1:12), cex = 0.8)
dev.off()
message("Saved heatmap: ", hm_path)

#Exercise 3
cl_dir <- file.path(outdir, "clusters_for_GO")
dir.create(cl_dir, showWarnings = FALSE)
cluster_ids <- split(rownames(filt), labels)
for (k in sort(as.integer(names(cluster_ids)))) {
  fn <- file.path(cl_dir, sprintf("cluster_%02d_genes.txt", k))
  writeLines(cluster_ids[[as.character(k)]], con = fn)
}
write.csv(
  tibble(cluster = as.integer(names(cluster_ids)),
         n_genes = vapply(cluster_ids, length, integer(1))) |>
    arrange(cluster),
  file.path(outdir, "cluster_sizes.csv"),
  row.names = FALSE
)
message("cluster gene", cl_dir)

