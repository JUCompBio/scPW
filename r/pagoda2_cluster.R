# Pagoda2 clustering, run as its authors intend.
#
# This runs the canonical Pagoda2 walkthrough on the same QC-filtered matrix the
# Python arms use, so only the method differs.
#
# Two deliberate differences from the other arms, both reported rather than hidden.
#
# Resolution. Pagoda2 clusters with igraph community detection (multilevel by
# default) and exposes no resolution parameter, so the resolution used by the other
# arms cannot be matched. Forcing a foreign routine onto Pagoda2 would make it a
# worse representative of itself than accepting its default does, so the default is
# used and the resulting k is reported next to every metric. k matters because these
# indices respond to partition granularity as well as to structure.
#
# Distance. Pagoda2 builds its kNN graph with cosine distance. Scoring its partition
# with euclidean distance measures it in a metric it never used, and that choice
# alone can move Silhouette substantially. Both are therefore computed and both are
# written out; see r/metrics_common.R.
#
# Pipeline, all parameters explicit:
#   basicP2proc     n.odgenes 4000        matched to our HVG count
#                   nPcs      50          matched to Seurat's npcs
#                   k         30          Pagoda2 default
#                   log.scale TRUE        Pagoda2 default
#                   trim      10          Pagoda2 default
#                   get.largevis/tsne FALSE   embeddings unused here, and both are slow
#   makeKnnGraph    k 30, distance cosine, center TRUE   Pagoda2 defaults
#   getKnnClusters  igraph::multilevel.community          Pagoda2 default

suppressPackageStartupMessages({
  library(pagoda2)
  library(Matrix)
  library(igraph)
})

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                               value = TRUE)[1]))
source(file.path(script_dir, "metrics_common.R"))

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]
out_dir <- args[2]
dataset <- args[3]
seed <- as.integer(args[4])
# args[5] is resolution, accepted for a uniform calling convention across the R arms but
# unused here: Pagoda2's community detection has no resolution parameter.
n_cores <- if (length(args) >= 6) as.integer(args[6]) else 1L

set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mat <- readMM(file.path(in_dir, "matrix.mtx"))
rownames(mat) <- readLines(file.path(in_dir, "features.tsv"))
colnames(mat) <- readLines(file.path(in_dir, "barcodes.tsv"))
mat <- as(mat, "dgCMatrix")

cat(sprintf("[%s] %d genes x %d cells\n", dataset, nrow(mat), ncol(mat)))

n_od <- min(4000, nrow(mat))
n_pcs <- min(50, ncol(mat) - 1, n_od - 1)

# Pagoda2's winsorization (trim) writes NAs into the count matrix when a cell has very
# few detected genes, and because the NAs land in the gene-mean calculation a single
# bad cell can remove every gene from the variance fit. basicP2proc then dies in
# lm.fit with "0 (non-NA) cases".
#
# Two fixes are possible: disable trim everywhere, or drop the offending cells.
# Disabling trim changes the result on datasets where trim works, so it would alter
# every dataset to accommodate one. Dropping the cells changes only the dataset that
# fails, and usually by a handful of cells, so the shared cell set is preserved almost
# exactly. That is the smaller intervention and the count is recorded in the output.
#
# The cells are identified by running Pagoda2's own normalization once and looking for
# NAs, rather than by guessing a count threshold, so the criterion is exactly "this cell
# breaks the variance fit".
probe <- Pagoda2$new(mat, n.cores = n_cores, trim = 10, log.scale = TRUE,
                     min.cells.per.gene = 0, min.transcripts.per.cell = 0,
                     verbose = FALSE)
bad <- rownames(probe$counts)[Matrix::rowSums(is.na(probe$counts)) > 0]
rm(probe)
n_dropped <- length(bad)
if (n_dropped > 0) {
  cat(sprintf("[%s] dropping %d cell(s) that make winsorization emit NAs: %s\n",
              dataset, n_dropped, paste(head(bad, 5), collapse = ", ")))
  mat <- mat[, setdiff(colnames(mat), bad), drop = FALSE]
}

# min.cells.per.gene and min.transcripts.per.cell are 0 on purpose: the input is already
# QC-filtered identically to every other arm, and letting Pagoda2 filter again would
# change the cell and gene set and break the comparison.
p2 <- basicP2proc(
  mat, n.cores = n_cores, n.odgenes = n_od, nPcs = n_pcs, k = 30,
  log.scale = TRUE, trim = 10,
  min.cells.per.gene = 0, min.transcripts.per.cell = 0,
  get.largevis = FALSE, get.tsne = FALSE, make.geneknn = FALSE
)

p2$makeKnnGraph(k = 30, type = "PCA", center = TRUE, distance = "cosine",
                n.cores = n_cores)
p2$getKnnClusters(method = igraph::multilevel.community, type = "PCA",
                  name = "community")

emb <- p2$reductions$PCA
# Align by name rather than trusting both objects to be in input order.
clusters <- p2$clusters$PCA$community[rownames(emb)]
stopifnot(!any(is.na(clusters)))
labels <- as.integer(as.factor(as.character(clusters)))

m <- cluster_metrics(emb, labels, seed = seed)

# Modularity on Pagoda2's own graph: metric-free, and the fairest single statement about
# whether the partition is good in the structure the method actually used.
g <- p2$graphs$PCA
mod <- if (!is.null(g)) {
  modularity(g, labels[match(V(g)$name, rownames(emb))])
} else NA_real_

out <- cbind(
  data.frame(dataset = dataset, method = "pagoda2", seed = seed,
             resolution = NA_real_,
             clustering = "igraph::multilevel.community",
             native_distance = "cosine",
             npcs = n_pcs, n_odgenes = n_od,
             cells_dropped = n_dropped,
             modularity = mod, stringsAsFactors = FALSE),
  m
)
write.csv(out, file.path(out_dir, paste0(dataset, "_metrics.csv")),
          row.names = FALSE)

write.csv(data.frame(cell = rownames(emb), cluster = labels),
          file.path(out_dir, paste0(dataset, "_clusters.csv")), row.names = FALSE)

# Differential expression per cluster, so this arm can be compared on markers and not
# only on cluster geometry.
ok <- tryCatch({
  p2$getDifferentialGenes(type = "PCA", clusterType = "community", verbose = FALSE)
  TRUE
}, error = function(e) { cat("  DE failed:", conditionMessage(e), "\n"); FALSE })

if (ok) {
  de <- p2$diffgenes$PCA$community
  if (!is.null(de) && length(de)) {
    rows <- do.call(rbind, lapply(names(de), function(cl) {
      d <- de[[cl]]
      if (is.null(d) || !nrow(d)) return(NULL)
      d <- d[order(-d$Z), , drop = FALSE]
      keep <- seq_len(min(200L, nrow(d)))
      data.frame(
        cluster = cl, rank = keep,
        gene = if ("Gene" %in% names(d)) as.character(d$Gene)[keep]
               else rownames(d)[keep],
        Z = d$Z[keep],
        M = if ("M" %in% names(d)) d$M[keep] else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
    write.csv(rows, file.path(out_dir, paste0(dataset, "_markers.csv")),
              row.names = FALSE)
  }
}

cat(sprintf("DONE %s sil_e=%.4f sil_c=%.4f ch=%.1f db=%.4f k=%d mod=%.3f cells=%d\n",
            dataset, m$silhouette, m$silhouette_cosine, m$calinski_harabasz,
            m$davies_bouldin, m$n_clusters, mod, m$n_cells))
