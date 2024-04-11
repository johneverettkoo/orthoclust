#' Compute the affinity matrix for orthogonal spectral clustering
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
#' @return B, the affinity matrix.
#'   \eqn{B[i, j]} approaches 0 asymptotically
#'   if i and j are in different communities.
#'   V, first p and last q eigenvectors of A.
#'   The rows of V that correspond to vertices in two different communities
#'   approach orthogonality asymptotically.
#'   p, q, the embedding dimensions.
#' @export
#'
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
#' # compute the affinity matrix
#' B <- orthoclust::compute.affinity.matrix(A, K = 2, model = 'sbm')
compute.affinity.matrix <- function(A, K = 2L,
                                    model = 'pabm',
                                    assortative = TRUE) {
  # verify that the model is correct
  # and make necessary transformations
  model.params <- .model.checks(model, K, assortative)
  K <- model.params$K
  p <- model.params$p
  q <- model.params$q

  # verify that A is a valid adjacency matrix
  # also compute n, the sample size
  A.n <- .check.adj.matrix(A)
  A <- A.n$A
  n <- A.n$n

  # construct the unscaled embedding
  V <- .ase(A, n, p, q, scale = FALSE)

  # compute the affinity matrix from the embedding
  B <- abs(tcrossprod(V))

  return(list(B = B,
              V = V,
              p = p, q = q))
}

.model.checks <- function(model, K, assortative) {
  # helper function for verifying that parameters are correct
  # also determines the dimensionality of the embedding
  # which depends on the type of model (SBM, DCBM, PABM)
  # and whether the model is assortative

  # verify K is an integer greater than or equal to 2
  if (!is.integer(K)) {
    if (is.numeric(K)) {
      K <- as.integer(K)
    } else {
      stop('K must be an integer')
    }
  }

  # compute p and q based on the model
  model <- tolower(model)
  if (model %in% c('sbm', 'dcbm', 'dcsbm')) {
    if (is.logical(assortative)) {
      if (assortative) {
        p <- K
        q <- 0L
      } else {
        p <- 1L
        q <- K - p
      }
    } else if (is.numeric(assortative)) {
      if (assortative > K) {
        stop('number of assortative dimensions must be <= K for SBM and DCBM')
      }
      p <- as.integer(assortative)
      q <- K - p
    } else {
      stop('assortative must be a logical or integer')
    }
  } else if (model %in% c('pabm', 'pasbm')) {
    p <- (K * (K + 1L)) %/% 2L
    q <- (K * (K - 1L)) %/% 2L
  } else {
    stop('model must be SBM, DCBM, or PABM')
  }

  return(list(K = K, p = p, q = q))
}
