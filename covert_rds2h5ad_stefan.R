

Sys.setenv(LANG="en")
library(spacexr)
library(Matrix)
library(Seurat)

seurat_object  <- readRDS("seu_liver_ref.rds")

# 1. Extract Gene Names
gene_names <- rownames(seurat_object)
# 2. Extract Cell Names (Barcodes)
cell_names <- colnames(seurat_object)

# Save gene names to a CSV file
write.csv(gene_names, file = "gene_names.csv", row.names = FALSE,quote=FALSE)



cluster_names <- Idents(seurat_object)

# Alternatively, if the cluster names are stored in a specific metadata column (e.g., 'seurat_clusters'):
# cluster_names <- seurat_object@meta.data$seurat_clusters

# Combine everything into a data frame for easy export
# Ensure the cluster names are matched to the cell names properly
cell_metadata <- data.frame(Cell_Name = cell_names, Cluster_Name = cluster_names)

# Save cell names and cluster names to a CSV file
write.csv(cell_metadata, file = "cell_metadata.csv", row.names = FALSE, quote=FALSE)



counts_matrix <- GetAssayData(seurat_object, slot = "counts")
counts_matrix <- as(counts_matrix, "sparseMatrix")
writeMM(counts_matrix, file = "counts_matrix.mtx")
umap_coordinates <- Embeddings(seurat_object, "umap")
write.csv(umap_coordinates, file = "umap_coordinates.csv", row.names = TRUE,quote=FALSE)
