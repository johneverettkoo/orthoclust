#' Adjacency Spectral Embedding
#'
#' Adjacency spectral embedding of undirected graph with adjacency matrix A.
#'
#'
#' @param A The adjacency matrix of the graph to be embedded
#' @param p Number of assortative dimensions
#' @param q Number of disassortative dimensions
#' @param scale Whether to scale the embedding by the eigenvalues
#'
#' @return An nrow(A) by (p+q) matrix of embedded points
#' @export
#'
#' @examples
#' # generate a small homogeneous SBM with two communities
#' n1 <- 20
#' n2 <- 20
#' n <- n1 + n2
#' p.within <- .5
#' q.between <- .1
#' P <- matrix(p.between, nrow = n, ncol = n)
#' P[seq(n1), seq(n1)] <- p.within
#' P[seq(n1 + 1, n), seq(n1 + 1, n)] <- p.within
#' A <- osc::draw.graph(P)
#'
#' # embed A in two assortative dimensions
#' X <- ase(A, 2, 0)
#' # the embedding will consist of points around two point masses
#' plot(X)
ase <- function(A, p = 2L, q = 0L,
                scale = TRUE) {

  # make sure function arguments are valid
  # and make changes if necessary
  Apq <- .ase.checks(A, p, q)
  A <- Apq$A
  n <- Apq$n
  p <- Apq$p
  q <- Apq$q

  # embed
  X <- .ase(A, n, p, q, scale)

  return(X)
}

.check.adj.matrix <- function(A) {
  # helper function for verifying that A is valid

  # verify that A is a matrix
  if (!('matrix' %in% class(A))) {
    stop('A must be a matrix')
  }
  if (!(typeof(A) %in% c('integer', 'double'))) {
    stop('A must be numeric')
  }

  # verify A is square
  dim.A <- dim(A)
  if (diff(dim.A) != 0) {
    stop('A must be a square matrix')
  }

  # verify A is symmetric
  if (any(A != t(A))) {
    stop('A must be symmetric')
  }

  # verify entries of A are nonnegative
  if (any(A < 0)) {
    stop('A must contain only nonnegative entries')
  }

  # force the graph to have no self loops
  if (any(diag(A) != 0)) {
    warning('A is not hollow')
  }

  return(list(A = A,
              n = dim.A[1]))
}

.ase.checks <- function(A, p, q) {
  # helper function for verifying that A, p, and q are valid

  # verify A
  A.n <- .check.adj.matrix(A)
  A <- A.n$A
  n <- A.n$n

  # verify p and q are integers
  if (class(p) != 'integer') {
    if (class(p) == 'numeric') {
      p <- as.integer(p)
    } else {
      stop('p must be an integer')
    }
  }
  if (class(q) != 'integer') {
    if (class(q) == 'numeric') {
      q <- as.integer(q)
    } else {
      stop('q must be an integer')
    }
  }
  if (length(p) > 1) {
    warning('using only the first entry of p')
    p <- p[1]
  }
  if (length(q) > 1) {
    warning('using only the first entry of q')
    q <- q[1]
  }

  # verify the embedding dimensions are valid
  if (p * q < 0 | p + q <= 0) {
    stop('One of p or q must be > nonzero, and both p and q must be positive')
  }

  # make sure p + q < n
  if (n <= p + q) {
    stop('embedding dimensions must be less than the sample size')
  }

  return(list(A = A, n = n, p = p, q = q))
}

.ase <- function(A, n, p, q, scale) {
  # base adjacency spectral embedding function without checks

  # compute the spectral decomposition of A
  eigen.A <- eigen(A, symmetric = TRUE)

  # take the p most positive and q most negative eigenvectors
  if (p * q > 0) {
    keep <- c(seq(p), seq(n, n - q + 1))
  } else if (p == 0) {
    keep <- seq(n, n - q + 1)
  } else {
    keep <- seq(p)
  }
  U <- eigen.A$vectors[, keep]

  # if we are scaling the embedding,
  # multiply by the square roots of the corresponding eigenvalues
  # otherwise normalize by the square root of the sample size
  if (scale) {
    S <- diag(sqrt(abs(eigen.A$values[keep])))
    X <- U %*% S
  } else {
    # this is to ensure that the entries of U do not vanish as n -> Inf
    X <- U * sqrt(n)
  }

  return(X)
}
