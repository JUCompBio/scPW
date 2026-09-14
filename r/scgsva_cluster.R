# scGSVA: pathway scoring, then clustering on the pathway scores.
#
# scGSVA is not a clustering method. It produces a pathway x cell score matrix, so
# something has to turn those scores into cluster labels before any clustering metric can
# be computed. That step is a choice made on scGSVA's behalf, and the published table
# reports scGSVA Silhouette values without recording what the choice was.
#
# The choice made here is the minimal one: cluster the pathway-score matrix with the same
# neighbours-plus-Louvain procedure and the same resolution used by every other arm, so
# the only difference from the Seurat row is the feature space (pathways instead of
# genes). It is stated in the output and belongs in the caption.
#
# Consequence worth stating plainly: the metrics come out in pathway space, whose
# dimensionality is the number of gene sets rather than 50 PCs. Silhouette is a distance
# ratio and moves with dimensionality, so this row is not directly comparable with the
# gene-space rows. Both euclidean and cosine Silhouette are recorded, as for every arm,
# and n_dims is in the CSV so the reader can see what space produced the number.
#
# Pipeline, all parameters explicit:
#   NormalizeData    LogNormalize, 1e4                 same as the Seurat arm
#   buildAnnot       species human, KEGG               pathway source
#   scgsva           method ssgsea                     fast, and the GSVA default family
#   PCA on scores    npcs 50 (capped by pathway count)
#   FindNeighbors    dims 1:npcs
#   FindClusters     resolution as passed               matched to every other arm

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scGSVA)
})

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "metrics_common.R"))

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]
out_dir <- args[2]
dataset <- args[3]
seed <- as.integer(args[4])
resolution <- as.numeric(args[5])
n_cores <- if (length(args) >= 6) as.integer(args[6]) else 1L

set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mat <- readMM(file.path(in_dir, "matrix.mtx"))
rownames(mat) <- readLines(file.path(in_dir, "features.tsv"))
colnames(mat) <- readLines(file.path(in_dir, "barcodes.tsv"))
mat <- as(mat, "dgCMatrix")

cat(sprintf("[%s] %d genes x %d cells\n", dataset, nrow(mat), ncol(mat)))

obj <- CreateSeuratObject(counts = mat, project = dataset)
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

  # KEGG gene sets for human, the method's own pathway source.
annot <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

sc <- scgsva(obj, annot, method = "ssgsea", cores = n_cores, verbose = FALSE)

# scGSVA returns an object carrying the score matrix; accept either a matrix or a slot,
# since the accessor has changed across versions and a silent failure here would be
# reported as a clustering result.
scores <- if (is.matrix(sc) || inherits(sc, "Matrix")) {
  as.matrix(sc)
} else if (!is.null(sc@gsva)) {
  as.matrix(sc@gsva)
} else {
  stop("cannot locate the score matrix in the scgsva result")
}

# Orient as cells x pathways.
if (ncol(scores) == ncol(obj) && nrow(scores) != ncol(obj)) {
  scores <- t(scores)
}
stopifnot(nrow(scores) == ncol(obj))
cat(sprintf("[%s] %d cells x %d pathways\n", dataset, nrow(scores), ncol(scores)))

# Drop pathways with no variance: they carry no information and would make the PCA
# rank-deficient.
keep <- apply(scores, 2, function(v) sd(v) > 0)
scores <- scores[, keep, drop = FALSE]
cat(sprintf("[%s] %d pathways with non-zero variance\n", dataset, ncol(scores)))

# Cluster in pathway space via PCA, matching the Seurat arm's procedure so only the
# feature space differs.
n_pcs <- min(50, ncol(scores) - 1, nrow(scores) - 1)
pca <- prcomp(scores, center = TRUE, scale. = FALSE, rank. = n_pcs)
emb <- pca$x
rownames(emb) <- rownames(scores)

nn <- FindNeighbors(emb, dims = 1:n_pcs, verbose = FALSE)
cl <- FindClusters(nn$snn, resolution = resolution, random.seed = seed,
                   verbose = FALSE)
labels <- as.integer(as.character(cl[[1]])) + 1L

m <- cluster_metrics(emb, labels, seed = seed)

out <- cbind(
  data.frame(dataset = dataset, method = "scgsva", seed = seed,
             resolution = resolution,
             clustering = "PCA + SNN + Louvain on pathway scores",
             native_distance = "euclidean",
             feature_space = "KEGG pathway scores (ssgsea)",
             n_pathways = ncol(scores), npcs = n_pcs,
             stringsAsFactors = FALSE),
  m
)
write.csv(out, file.path(out_dir, paste0(dataset, "_metrics.csv")),
          row.names = FALSE)

write.csv(data.frame(cell = rownames(emb), cluster = labels),
          file.path(out_dir, paste0(dataset, "_clusters.csv")), row.names = FALSE)

# Which pathways separate which cluster. The analogue of markers for this arm, and what
# makes the row interpretable rather than only scoreable.
means <- t(apply(scores, 2, function(v) tapply(v, labels, mean)))
rows <- do.call(rbind, lapply(sort(unique(labels)), function(cl_i) {
  own <- means[, as.character(cl_i)]
  rest <- rowMeans(means[, colnames(means) != as.character(cl_i), drop = FALSE])
  d <- own - rest
  o <- order(-d)[seq_len(min(50L, length(d)))]
  data.frame(cluster = cl_i, rank = seq_along(o),
             pathway = rownames(means)[o],
             mean_in_cluster = own[o], mean_elsewhere = rest[o],
             difference = d[o], stringsAsFactors = FALSE)
}))
write.csv(rows, file.path(out_dir, paste0(dataset, "_pathways.csv")),
          row.names = FALSE)

cat(sprintf("DONE %s sil_e=%.4f sil_c=%.4f ch=%.1f db=%.4f k=%d paths=%d cells=%d\n",
            dataset, m$silhouette, m$silhouette_cosine, m$calinski_harabasz,
            m$davies_bouldin, m$n_clusters, ncol(scores), m$n_cells))
