% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_graph.R
\name{draw.graph}
\alias{draw.graph}
\title{Draw the adjacency matrix of a graph from edge probability matrix P.}
\usage{
draw.graph(P, distribution = "bernoulli")
}
\arguments{
\item{P}{The symmetric edge probability matrix from which
we are sampling the graph}

\item{distribution}{A string indicating the distribution to draw from.
Currently supports bernoulli, poisson, and exponential.}
}
\value{
A the adjacency matrix sampled from P.
}
\description{
Independently draw \eqn{A[i, j] ~ F(P[i, j])} for distribution F and i > j.
A is then symmetrized.
}
\examples{
# generate a small homogeneous SBM with two communities
n1 <- 20
n2 <- 20
n <- n1 + n2
p.within <- .5
p.between <- .1
P <- matrix(p.between, nrow = n, ncol = n)
P[seq(n1), seq(n1)] <- p.within
P[seq(n1 + 1, n), seq(n1 + 1, n)] <- p.within
A <- orthoclust::draw.graph(P)

# visualize
qgraph::qgraph(A, layout = 'spring')
}
