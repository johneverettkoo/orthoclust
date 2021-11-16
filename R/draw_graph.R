#' Draw the adjacency matrix of a graph from edge probability matrix P.
#'
#' Independently draw \eqn{A[i, j] ~ F(P[i, j])} for distribution F and i > j.
#' A is then symmetrized.
#'
#' @param P The symmetric edge probability matrix from which
#'   we are sampling the graph
#' @param distribution A string indicating the distribution to draw from.
#'   Currently supports bernoulli, poisson, and exponential.
#'
#' @return A the adjacency matrix sampled from P.
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
#' A <- osc::draw.graph(P)
#'
#' # visualize
#' qgraph::qgraph(A, layout = 'spring')
draw.graph <- function(P,
                       distribution = 'bernoulli') {
  # verify that the parameters are valid
  valid.params <- .draw.graph.checks(P, distribution)
  P <- valid.params$P
  n <- valid.params$n
  rdist <- valid.params$rdist

  # draw the adjacency matrix
  A <- .draw.graph(P, n, rdist)

  return(A)
}

.draw.graph.checks <- function(P, distribution) {
  # helper function for verifying that P and distribution are valid
  # also converts distribution (character) to an appropriate function

  # verify P is a matrix
  if (!('matrix') %in% class(P)) {
    stop('P must be a matrix')
  }

  # verify P is square
  dim.P <- dim(P)
  if (diff(dim.P) != 0) {
    stop('P must be a square matrix')
  }

  # verify P is symmetric
  if (any(P != t(P))) {
    warning('P is not symmetric--ignoring the lower triangular portion of P')
  }

  # verify entries of P are nonnegative
  if (any(P < 0)) {
    stop('P must contain only nonnegative entries')
  }

  # verify that distribution is valid
  # and choose the corresponding rng
  if (!('character' %in% class(distribution))) {
    stop('distribution must be of type character')
  }
  if (length(distribution) > 1) {
    distribution <- distribution[1]
    warning('distribution has length > 1; only the first element will be used')
  }
  if (distribution == 'bernoulli') {
    rdist <- .draw.graph.bernoulli
  } else if (distribution == 'poisson') {
    rdist <- .draw.graph.poisson
  } else if (distribution == 'exponential') {
    rdist <- .draw.graph.exponential
  } else {
    stop('choose a supported distribution')
  }

  # obtain the sample size
  n <- dim.P[1]

  return(list(P = P, n = n, rdist = rdist))
}

.draw.graph <- function(P, n, rdist) {
  # initialize the adjacency matrix with zeros
  A <- matrix(0, nrow = n, ncol = n)

  # fill each A[i, j] by drawing from each P[i, j]
  A[upper.tri(A)] <- rdist(P, n)

  # symmeterize
  A <- A + t(A)

  return(A)
}

# the following are helper functions for drawing random numbers
# from the upper triangular portion of matrix P

.draw.graph.bernoulli <- function(P, n) {
  rbinom(n * (n - 1) / 2, 1, P[upper.tri(P)])
}

.draw.graph.poisson <- function(P, n) {
  rpois(n * (n - 1) / 2, P[upper.tri(P)])
}

.draw.graph.exponential <- function(P, n) {
  rexp(n * (n - 1) / 2, 1 / P[upper.tri(P)])
}
