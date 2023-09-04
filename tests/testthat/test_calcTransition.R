library(MigConnectivity)
context('Calculate transition probabilities')

test_that('telemetry data produces right transition probabilities', {
  psiMats <- list(matrix(0.25, 2, 4, dimnames = list(LETTERS[1:2], as.character(1:4))),
                  matrix(c(0.2, 0.8, 0.875, 0.125, 1/9, 8/9, 0.8, 0.2, 0.2, 0.8),
                         5, 2, TRUE, list(LETTERS[1:5], 1:2)))

  expect_equal(calcTransition(originAssignment = rep(1:2, each = 8),
                              targetAssignment = rep(1:4, 4),
                              originNames = LETTERS[1:2],
                              targetNames = as.character(1:4))$psi,
               psiMats[[1]])
  expect_equal(calcTransition(counts = matrix(c(2, 7, 1, 8, 2, 8, 1, 8, 2, 8),
                                              5, 2, FALSE,
                                              list(LETTERS[1:5], 1:2)))$psi,
               psiMats[[2]])
})
