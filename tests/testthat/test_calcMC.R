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
  breedingRelN <- breedingN <- vector("list", nScenarios)
  # Distances
  genericD <- matrix(c(0:3, 1, 0, 1, 2, 2, 1, 0, 1, 3:0), 4, 4)
  breedingD <- winteringD <- vector("list", nScenarios)
  for (i in 1:nScenarios) {
    breedingN[[i]] <- rep(10000, nBreeding[i])
    breedingRelN[[i]] <- breedingN[[i]] / sum(breedingN[[i]])
    breedingD[[i]] <- genericD[1:nBreeding[i], 1:nBreeding[i]]
    winteringD[[i]] <- genericD[1:nWintering[i], 1:nWintering[i]]
  }

  expect_equal(calcMC(originDist = breedingD[[1]], targetDist = winteringD[[1]],
                      originRelAbund = breedingRelN[[1]], psi = psiMats[[1]]),
               0)
  expect_equal(calcMC(originDist = breedingD[[2]], targetDist = winteringD[[2]],
                      originRelAbund = breedingRelN[[2]], psi = psiMats[[2]]),
               0)
  expect_equal(calcMC(originDist = breedingD[[3]], targetDist = winteringD[[3]],
                      originRelAbund = breedingRelN[[3]], psi = psiMats[[3]]),
               1)
  expect_equal(calcMC(originDist = breedingD[[4]], targetDist = winteringD[[4]],
                      originRelAbund = breedingRelN[[4]], psi = psiMats[[4]]),
               0.6)
  expect_equal(calcMC(originDist = breedingD[[5]], targetDist = winteringD[[5]],
                      originRelAbund = breedingRelN[[5]], psi = psiMats[[5]]),
               0.196)
  expect_equal(calcMC(originDist = breedingD[[6]], targetDist = winteringD[[6]],
                      originRelAbund = breedingRelN[[6]], psi = psiMats[[6]]),
               0.504)
  expect_equal(calcMC(originDist = breedingD[[7]], targetDist = winteringD[[7]],
                      originRelAbund = breedingRelN[[7]], psi = psiMats[[7]]),
               0.164144856)
  expect_equal(calcMC(originDist = breedingD[[8]], targetDist = winteringD[[8]],
                      originRelAbund = breedingRelN[[8]], psi = psiMats[[8]]),
               -6.656e-02)
  expect_equal(calcMC(originDist = breedingD[[1]], targetDist = winteringD[[1]],
                      originRelAbund = breedingRelN[[1]], psi = psiMats[[1]],
                      sampleSize = sum(breedingN[[1]])),
               -0.000041669)
  expect_equal(calcMC(originDist = breedingD[[2]], targetDist = winteringD[[2]],
                      originRelAbund = breedingRelN[[2]], psi = psiMats[[2]],
                      sampleSize = sum(breedingN[[2]])),
               -0.000039222)
  expect_equal(calcMC(originDist = breedingD[[3]], targetDist = winteringD[[3]],
                      originRelAbund = breedingRelN[[3]], psi = psiMats[[3]],
                      sampleSize = sum(breedingN[[1]])),
               1)
  expect_equal(calcMC(originDist = breedingD[[4]], targetDist = winteringD[[4]],
                      originRelAbund = breedingRelN[[4]], psi = psiMats[[4]],
                      sampleSize = sum(breedingN[[1]])),
               0.599983332)
  expect_equal(calcMC(originDist = breedingD[[5]], targetDist = winteringD[[5]],
                      originRelAbund = breedingRelN[[5]], psi = psiMats[[5]],
                      sampleSize = sum(breedingN[[1]])),
               0.195966498)
  expect_equal(calcMC(originDist = breedingD[[6]], targetDist = winteringD[[6]],
                      originRelAbund = breedingRelN[[6]], psi = psiMats[[6]],
                      sampleSize = sum(breedingN[[1]])),
               0.503979332)
  expect_equal(calcMC(originDist = breedingD[[7]], targetDist = winteringD[[7]],
                      originRelAbund = breedingRelN[[7]], psi = psiMats[[7]],
                      sampleSize = sum(breedingN[[1]])),
               0.164112566)
  expect_equal(calcMC(originDist = breedingD[[8]], targetDist = winteringD[[8]],
                      originRelAbund = breedingRelN[[8]], psi = psiMats[[8]],
                      sampleSize = sum(breedingN[[1]])),
               -0.066604443)
})
