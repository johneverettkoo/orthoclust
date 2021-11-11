ase <- function(A, p = 2, q = 0,
                scale = TRUE,
                eps = 1e-6) {

  # make sure function arguments are valid
  # and make changes if necessary
  Apq <- .ase.checks(A, p, q)
  A <- Apq$A
  p <- Apq$p
  q <- Apq$q

  # determine the number of elements of A
  n <- nrow(A)

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

.ase.checks <- function(A, p, q) {
  # helper function for verifying that A, p, and q are valid

  # verify A is a matrix and p and q are integers
  if (!('matrix' %in% class(A))) {
    stop('A must be a matrix')
  }
  if (!(typeof(A) %in% c('integer', 'double'))) {
    stop('A must be numeric')
  }
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
  if (p * q <= 0 | p + q <= 0) {
    stop('One of p or q must be > nonzero, and both p and q must be positive')
  }

  # verify A is square
  if (diff(dim(A)) != 0) {
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
    warning('forcing A to be hollow')
    diag(A) <- 0
  }

  # make sure p + q < n
  if (nrow(A) <= p + q) {
    stop('embedding dimensions must be less than the sample size')
  }

  return(list(A = A, p = p, q = q))
}
