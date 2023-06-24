#' Construct the edge probability matrix for a PABM
#'
#' @param groups A vector of group labels.
#' @param popularity.params An n by K matrix of popularity parameters.
#'   \expr{popularity.params[i, k]} is vertex i's affinity to community k.
#'
#' @return The n by n edge probability matrix.
#'
#' @examples
#' # set the size of each group
#' n1 <- 100
#' n2 <- 100
#' n3 <- 100
#' n <- n1 + n2 + n3
pabm.edge.probability.matrix <- function(groups,
                                         popularity.params) {
  # check that the group labels are valid
  group.params <- .check.group.labels(groups)
  groups <- group.params$groups
  n <- group.params$n
  K <- group.params$K
  n.group <- group.params$n.group

  # check that the popularity parameters are valid
  popularity.params <- .check.popularity.params(popularity.params, n, K)

  # construct edge probability matrix
  P <- .pabm.edge.probability.matrix(groups, popularity.params, K, n, n.group)
  return(P)
}

.check.group.labels <- function(groups) {
  # helper function for verifying that the group labels are valid

  # groups must be a vector
  if (!is.vector(groups)) {
    stop('groups must be a vector')
  }

  n <- length(groups)

  # groups must be integers
  if (!is.integer(groups)) {
    groups <- as.integer(groups)
  }

  # make sure the labels are 1...K
  K <- max(groups)
  unique.groups <- seq(K)
  n.group <- table(groups)  # size of each group
  if (length(unique.groups) != length(n.group)) {
    stop('groups must be labeled 1...K')
  }
  if (any(unique.groups != as.integer(names(n.group)))) {
    stop('groups must be labeled 1...K')
  }

  return(list(groups = groups,
              n = n,
              K = K,
              n.group = n.group))
}

.check.popularity.params <- function(popularity.params, n, K) {
  # helper function for verifying that the popularity params are valid

  # verify that popularity.params is a positive numeric n by K matrix
  if (!is.matrix(popularity.params)) {
    stop('popularity.params must be a matrix')
  }
  if (!is.numeric(popularity.params)) {
    stop('popularity.params must be numeric')
  }
  if (any(popularity.params < 0)) {
    warning('setting all negative popularity parameters to 0')
    popularity.params <- popularity.params[popularity.params < 0] <- 0
  }
  if (any(dim(popularity.params) != c(n, K))) {
    stop(paste('popularity parameters must be of dimension',
               'sample size by number of groups'))
  }

  return(popularity.params)
}

.pabm.edge.probability.matrix <- function(groups,
                                          popularity.params,
                                          K, n, n.group) {
  # construct edge probability matrix for the PABM

  # initialize
  P <- matrix(NA, n, n)

  # fill in P by group
  for (k in seq(K)) {
    for (l in seq(k)) {
      # if this is within-group
      if (k == l) {
        # extract the popularity params
        lambda <- popularity.params[groups == k, k]

        # construct the edge probabiliy block
        P.kk <- tcrossprod(lambda)
        P[groups == k, groups == k] <- P.kk
      # between-group
      } else {
        # extract the popularity params
        lambda.kl <- popularity.params[groups == k, l]
        lambda.lk <- popularity.params[groups == l, k]

        # construct the edge probability block
        P.kl <- tcrossprod(lambda.kl, lambda.lk)
        P[groups == k, groups == l] <- P.kl
        P[groups == l, groups == k] <- t(P.kl)
      }
    }
  }
  return(P)
}
