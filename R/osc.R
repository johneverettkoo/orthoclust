#' Orthogonal Spectral Clustering
#'
#' @param A Adjacency matrix describing an undirected graph.
#' @param K Number of communities.
#' @param model Type of block model (SBM, DCBM, or PABM).
#'   The dimensionality of the embedding of A depends on this model.
#'   If the model is SBM or DCBM, then the embedding dimension is K.
#'   If the model is PABM, then the embedding dimension is \eqn{K^2}
#'   and \eqn{p = K (K + 1) / 2} and \eqn{q = K (K - 1) / 2}.
#' @param assortative Whether the block model is assortative or disassortative.
#'   Only valid for SBM and DCBM. This argument is ignored for PABM.
#'   If assortative is TRUE, then p = K and q = 0.
#'   If assortative is FALSE, then p = 1 and q = K - 1.
#'   If assortative is an integer, then p is the value of that integer and
#'   q = K - p.
#' @return clustering, the cluster assignments, and B, the affinity matrix.
#'   affinity.matrix, the affinity matrix of OSC.
#'   embedding, the unscaled ASE of A.
#' @import mclust
#' @export
#' @examples
#' # generate a small homogeneous SBM with two communities
#' n1 <- 20
#' n2 <- 20
#' n <- n1 + n2
#' p.within <- .5
#' p.between <- .1
#' P <- matrix(p.between, nrow = n, ncol = n)
#' P[seq(n1), seq(n1)] <- p.within
#' P[seq(n1 + 1, n), seq(n1 + 1, n)] <- p.within
#' A <- orthoclust::draw.graph(P)
#'
#' # cluster
#' out <- orthoclust::osc(A, K = 2, model = 'sbm')
osc <- function(A,
                K = 2,
                model = 'pabm',
                assortative = TRUE) {
  # compute the affinity matrix of pairwise dot products
  affinity.out <- compute.affinity.matrix(A, K, model, assortative)
  B <- affinity.out$B
  V <- affinity.out$V
  p <- affinity.out$p
  q <- affinity.out$q
  n <- ncol(A)

  # compute the lapacian of the affinity matrix
  # and apply spectral clustering
  L <- .normalized.laplacian(B)
  eigenmap <- eigen(L, symmetric = TRUE)$vectors[, seq(n, n - p - q + 1)]
  clustering <- mclust::Mclust(eigenmap, K)$classification

  return(list(clustering = clustering,
              affinity.matrix = B,
              embedding = V))
}

.normalized.laplacian <- function(W) {
  # compute the normalized laplacian
  n <- nrow(W)
  d.neg.sqrt <- colSums(W) ** -.5
  d.neg.sqrt[is.infinite(d.neg.sqrt)] <- 0
  D.neg.sqrt <- diag(d.neg.sqrt)
  I <- diag(n)
  return(I - D.neg.sqrt %*% W %*% D.neg.sqrt)
}
