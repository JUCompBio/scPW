# Standard Seurat clustering workflow, run as-is with documented defaults.
#
# The canonical pipeline with nothing added, and every parameter set explicitly
# rather than left implicit:
#
#   NormalizeData      LogNormalize, scale.factor 10000   (Seurat default)
#   FindVariableFeatures vst, nfeatures 4000              matched to the HVG count
#   ScaleData          all selected features
#   RunPCA             npcs 50                            (Seurat default)
#   FindNeighbors      dims 1:50
#   FindClusters       resolution 0.55                    matched across arms
#
# The input is the same QC-filtered matrix the Python arms use, exported to MTX, so
# the comparison is not confounded by different preprocessing. Only the clustering
# method differs.
#
# Metrics are computed in R on the PCA embedding, which is Seurat's own
# representation. Each method is scored in the space it clusters on, which is also
# why cross-method Silhouette values are not strictly comparable: the same partition
# scores differently in a 5-d latent space than in 50-d PCA.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Metric definitions live in one place so all three R arms compute them identically
# and cannot drift apart. The helper also reports cosine Silhouette alongside
# euclidean, because the same partition can score very differently under the two and
# a table mixing them would be meaningless.
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "metrics_common.R"))

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]     # directory holding matrix.mtx, barcodes.tsv, features.tsv
out_dir <- args[2]
dataset <- args[3]
seed <- as.integer(args[4])
resolution <- as.numeric(args[5])

set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mat <- readMM(file.path(in_dir, "matrix.mtx"))
genes <- readLines(file.path(in_dir, "features.tsv"))
cells <- readLines(file.path(in_dir, "barcodes.tsv"))
rownames(mat) <- genes
colnames(mat) <- cells

obj <- CreateSeuratObject(counts = mat, project = dataset)
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = min(4000, nrow(obj)), verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

npcs <- min(50, ncol(obj) - 1, length(VariableFeatures(obj)) - 1)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs,
              seed.use = seed, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:npcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = resolution, random.seed = seed,
                    verbose = FALSE)

emb <- Embeddings(obj, "pca")[, 1:npcs, drop = FALSE]
labels <- as.integer(as.character(Idents(obj))) + 1L

m <- cluster_metrics(emb, labels, seed = seed)

write.csv(
  cbind(
    data.frame(dataset = dataset, method = "seurat", seed = seed,
               resolution = resolution,
               clustering = "Seurat FindClusters (Louvain)",
               native_distance = "euclidean",
               npcs = npcs, stringsAsFactors = FALSE),
    m
  ),
  file.path(out_dir, paste0(dataset, "_metrics.csv")), row.names = FALSE
)

# Cluster assignments, so the partition itself can be compared against the Python
# arms rather than only its summary statistics.
write.csv(data.frame(cell = colnames(obj), cluster = Idents(obj)),
          file.path(out_dir, paste0(dataset, "_clusters.csv")), row.names = FALSE)

# Markers, matching the 50 per cluster the Python arms save.
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1,
                          logfc.threshold = 0.25, verbose = FALSE)
# Seurat returns a 'gene' column in some versions and rownames in others.
if (!"gene" %in% colnames(markers) && nrow(markers) > 0) {
  markers$gene <- rownames(markers)
}
if (nrow(markers) > 0) {
  markers <- do.call(rbind, lapply(split(markers, markers$cluster),
                                   function(x) head(x[order(-x$avg_log2FC), ], 50)))
  write.csv(markers, file.path(out_dir, paste0(dataset, "_markers.csv")),
            row.names = FALSE)
}

# UMAP figure, so the directory holds a PDF like the other methods.
obj <- RunUMAP(obj, dims = 1:npcs, seed.use = seed, verbose = FALSE)
pdf(file.path(out_dir, paste0("umap_", dataset, ".pdf")), width = 6, height = 5)
print(DimPlot(obj, reduction = "umap", label = TRUE) +
        ggplot2::ggtitle(paste0(dataset, " - Seurat")))
dev.off()

cat(sprintf("DONE %s sil_e=%.4f sil_c=%.4f ch=%.1f db=%.4f k=%d cells=%d\n",
            dataset, m$silhouette, m$silhouette_cosine, m$calinski_harabasz,
            m$davies_bouldin, m$n_clusters, m$n_cells))
