library(testthat)
# library(osc)

# test_check("osc")

set.seed(12345)
n1 <- 20
n2 <- 20
n <- n1 + n2
p.within <- .5
p.between <- .1
P <- matrix(p.between, nrow = n, ncol = n)
P[seq(n1), seq(n1)] <- p.within
P[seq(n1 + 1, n), seq(n1 + 1, n)] <- p.within
A <- osc::draw.graph(P)

testthat::test_that('draw.graph outputs an adjacency matrix', {
  testthat::expect_equal(nrow(A), n)
  testthat::expect_equal(ncol(A), n)
  testthat::expect_equal(A, t(A))
  testthat::expect_equal(sum(abs(diag(A))), 0)
  testthat::expect_true(all(A >= 0))
})

out <- osc::osc(A, K = 2, model = 'sbm')

testthat::test_that('check for outputs of osc', {
  testthat::expect_true(all(out$clustering %in% seq(2)))
})
