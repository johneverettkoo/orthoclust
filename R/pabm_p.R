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

  # sort the group labels to make construction of P more straightforward
  sort.groups <- sort(groups)

  # fill in P by group
  for (k in seq(K)) {
    # determine the size of this group
    n.k <- n.group[k]

    # index range of P for this group (after sorting)
    low.ind.k <- ifelse(k == 1, 1, sum(n.group[seq(k - 1)]) + 1)
    high.ind.k <- sum(n.group[seq(k)])

    for (l in seq(k)) {
      n.l <- n.group[l]

      # if this is within-group
      if (k == l) {
        lambda <- popularity.params[sort.groups == k, k]
        P.kk <- tcrossprod(lambda)
        P[seq(low.ind.k, high.ind.k), seq(low.ind.k, high.ind.k)] <- P.kk
        # between-group
      } else {
        low.ind.l <- ifelse(l == 1, 1, sum(n.group[seq(l - 1)]) + 1)
        high.ind.l <- sum(n.group[seq(l)])
        lambda.kl <- popularity.params[sort.groups == k, l]
        lambda.lk <- popularity.params[sort.groups == l, k]
        P.kl <- tcrossprod(lambda.kl, lambda.lk)
        P[seq(low.ind.k, high.ind.k), seq(low.ind.l, high.ind.l)] <- P.kl
        P[seq(low.ind.l, high.ind.l), seq(low.ind.k, high.ind.k)] <- t(P.kl)
      }
    }
  }

  P <- P[order(groups), order(groups)]

  return(P)
}
