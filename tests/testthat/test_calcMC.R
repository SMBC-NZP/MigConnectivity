library(MigConnectivity)
context('Calculate migratory connectivity')

test_that('psi examples generate right MC values', {
  psiMats <- list(matrix(0.25, 4, 4), #"Full Mix"
                matrix(c(rep(0.32, 12), rep(0.04, 4)), 4, 4), #"Avoid One Site"
                diag(nrow=4),  #"Full Connectivity"
                matrix(c(rep(c(0.5, 0), 2, each=2), rep(c(0, 0.5), 2, each=2)),
                       4, 4), #"Half Mix"
                matrix(c(0.55, 0.2, 0.15, 0.1, 0.1, 0.55, 0.2, 0.15, 0.15, 0.1,
                         0.55, 0.2, 0.2, 0.15, 0.1, 0.55), 4, 4, byrow=T), #"Low"
                matrix(c(rep(c(0.75, 0.15, rep(0.05, 3)), 3), 0.75), 4, 4,
                       byrow=T), #"Medium"
                matrix(c(rep(0.25, 12), rep(0, 3), 1), 4, 4, byrow=T), #"Site Pref"
                matrix(c(0.01, 0.49, 0.49, 0.01, 0.49, 0.01, 0.01, 0.49, 0.49,
                         0.01, 0.01, 0.49, 0.01, 0.49, 0.49, 0.01), 4, 4,
                       byrow=T)) #Negative
  psiMats[[6]][4, 1] <- 0.15
  nBreeding <- sapply(psiMats, nrow)
  nWintering <- sapply(psiMats, ncol)
  nScenarios <- length(psiMats)

  # Relative abundances
  breedingRelN <- vector("list", nScenarios)
  # Distances
  genericD <- matrix(c(0:3, 1, 0, 1, 2, 2, 1, 0, 1, 3:0), 4, 4)
  breedingD <- winteringD <- vector("list", nScenarios)
  for (i in 1:nScenarios) {
    breedingRelN[[i]] <- rep(1/nBreeding[i], nBreeding[i])
    breedingD[[i]] <- genericD[1:nBreeding[i], 1:nBreeding[i]]
    winteringD[[i]] <- genericD[1:nWintering[i], 1:nWintering[i]]
  }

  expect_equal(calcMC(breedingD[[1]], winteringD[[1]], psiMats[[1]],
                      breedingRelN[[1]]), 0)
  expect_equal(calcMC(breedingD[[2]], winteringD[[2]], psiMats[[2]],
                      breedingRelN[[2]]), 0)
  expect_equal(calcMC(breedingD[[3]], winteringD[[3]], psiMats[[3]],
                      breedingRelN[[3]]), 1)
  expect_equal(calcMC(breedingD[[4]], winteringD[[4]], psiMats[[4]],
                      breedingRelN[[4]]), 0.6)
  expect_equal(calcMC(breedingD[[5]], winteringD[[5]], psiMats[[5]],
                      breedingRelN[[5]]), 0.196)
  expect_equal(calcMC(breedingD[[6]], winteringD[[6]], psiMats[[6]],
                      breedingRelN[[6]]), 0.504)
  expect_equal(calcMC(breedingD[[7]], winteringD[[7]], psiMats[[7]],
                      breedingRelN[[7]]), 0.164144856)
  expect_equal(calcMC(breedingD[[8]], winteringD[[8]], psiMats[[8]],
                      breedingRelN[[8]]), -6.656e-02)
})
