# Cluster-quality metrics, shared by every R arm.
#
# Sourced by r/seurat_cluster.R, r/pagoda2_cluster.R and r/scgsva_cluster.R so all three
# compute metrics identically. Definitions match sklearn's, so the values are comparable
# with the Python arms. clusterCrit is absent from these images and installing it at run
# time would add a network dependency for two formulas that are a dozen lines each.
#
# Both euclidean and cosine Silhouette are reported. This is not redundancy: the same
# partition in the same embedding can swing by more than an order of magnitude between
# the two, produced by the distance function alone. Methods also differ in which
# distance they use internally: Pagoda2 builds its kNN graph with cosine, Seurat and
# scVI work in euclidean space. Reporting one column silently privileges whichever
# method happens to match the chosen metric, so both go in the CSV.

suppressPackageStartupMessages({
  library(cluster)
})

calinski_harabasz <- function(x, labels) {
  n <- nrow(x); k <- length(unique(labels))
  if (k < 2 || k >= n) return(NA_real_)
  overall <- colMeans(x)
  between <- 0; within <- 0
  for (g in unique(labels)) {
    pts <- x[labels == g, , drop = FALSE]
    ctr <- colMeans(pts)
    between <- between + nrow(pts) * sum((ctr - overall)^2)
    within <- within + sum((sweep(pts, 2, ctr))^2)
  }
  if (within == 0) return(NA_real_)
  (between / (k - 1)) / (within / (n - k))
}

davies_bouldin <- function(x, labels) {
  groups <- unique(labels); k <- length(groups)
  if (k < 2) return(NA_real_)
  ctrs <- t(sapply(groups, function(g) colMeans(x[labels == g, , drop = FALSE])))
  spread <- sapply(seq_along(groups), function(i) {
    pts <- x[labels == groups[i], , drop = FALSE]
    mean(sqrt(rowSums((sweep(pts, 2, ctrs[i, ]))^2)))
  })
  worst <- sapply(seq_along(groups), function(i) {
    others <- setdiff(seq_along(groups), i)
    max(sapply(others, function(j) {
      d <- sqrt(sum((ctrs[i, ] - ctrs[j, ])^2))
      if (d == 0) Inf else (spread[i] + spread[j]) / d
    }))
  })
  mean(worst)
}

cosine_dist <- function(x) {
  nrm <- sqrt(rowSums(x^2)); nrm[nrm == 0] <- 1
  xn <- x / nrm
  # Clamp to [-1, 1] before subtracting: floating point can push a self-similarity
  # marginally above 1 and produce a negative distance.
  s <- pmin(pmax(xn %*% t(xn), -1), 1)
  as.dist(1 - s)
}

# Returns a one-row data frame. `subsample` caps the silhouette computation, since the
# full distance matrix is O(n^2) and the 32k-cell datasets would otherwise dominate
# runtime and memory. The cap that was applied is recorded in the output so a full and a
# subsampled value are never silently compared.
cluster_metrics <- function(emb, labels, subsample = 10000L, seed = 420L) {
  n <- nrow(emb)
  set.seed(seed)
  if (n > subsample) {
    idx <- sample.int(n, subsample)
    scope <- paste0("subsampled_", subsample)
  } else {
    idx <- seq_len(n)
    scope <- "full"
  }
  sub_emb <- emb[idx, , drop = FALSE]
  sub_lab <- labels[idx]

  ok <- length(unique(sub_lab)) > 1
  sil_e <- if (ok) mean(silhouette(sub_lab, dist(sub_emb))[, 3]) else NA_real_
  sil_c <- if (ok) mean(silhouette(sub_lab, cosine_dist(sub_emb))[, 3]) else NA_real_

  data.frame(
    silhouette = sil_e,              # euclidean, the primary reported column
    silhouette_cosine = sil_c,
    silhouette_scope = scope,
    calinski_harabasz = calinski_harabasz(emb, labels),
    davies_bouldin = davies_bouldin(emb, labels),
    n_clusters = length(unique(labels)),
    n_cells = n,
    n_dims = ncol(emb),
    stringsAsFactors = FALSE
  )
}
