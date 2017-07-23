set.seed(75)

# Transfer between true populations and researcher defined ones
# (only for breeding, as not messing with winter populations here)
breedingSiteTrans16 <- list(1:100, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)), 1:100,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))

#lapply(breedingSiteTrans16, matrix, nrow=10, ncol=10)

nScenarios16 <- length(breedingSiteTrans16)
nSims16 <- 100
# Basing positions of researcher defined breeding populations on above
breedingPos16 <- breedingPos15
winteringPos16 <- winteringPos15

breedDist16 <- breedDist15
nonbreedDist16 <- nonbreedDist15
nBreeding16 <- nBreeding15
nWintering16 <- rep(100, nScenarios16)

# Highest abundance in lower right corner, lowest in top left
# In fact basing on distance from top left population

breedingN16base <- breedingN15base
breedingN16 <- breedingN15
breedingRelN16 <- lapply(breedingN16, "/", sum(breedingN16base))

# Set up psi matrix
# Each quadrant of breeding range has different MC

MC.levels16 <- seq(0.15, 0.6, 0.15)

nLevels16 <- 4

psi16 <- matrix(NA, nBreeding16[1], nWintering16[1])


for (i in 1:nLevels16) {
  cat("MC", MC.levels16[i])
  # Find a psi matrix that produces the given MC (for whole species)
  o16a <- optimize(mlogitMC, MC.in = MC.levels16[i],
                   origin.dist = breedDist16[[1]],
                   target.dist = nonbreedDist16[[1]],
                   origin.abund = breedingN16[[1]]/sum(breedingN16[[1]]),
                   sample.size = sum(breedingN16[[1]]),
                   interval=c(0,10), tol=.Machine$double.eps^0.5)

  slope16a <- o16a$minimum

  cat(" slope", slope16a, "\n")

  psi16a <- mlogitMat(slope16a, breedDist16[[1]])

  # Then use the rows of that psi matrix only for the one breeding quadrant
  rows <- 50*(i %/% 3) + rep(1:5, 5) + rep(seq(0, 40, 10), each=5) + ((i-1) %% 2) * 5
  psi16[rows, ] <- psi16a[rows, ]
}


# Baseline strength of migratory connectivity
MC16 <- calcMC(originDist = breedDist16[[1]],
               targetDist = nonbreedDist16[[1]],
               psi = psi16,
               originRelAbund = breedingN16[[1]]/sum(breedingN16[[1]]),
               sampleSize = sum(breedingN16[[1]]))


# Set up sampling regimes (different number than number of scenarios)
nSample16 <- 100
sampleBreeding16 <- list(round(breedingRelN16[[1]]*nSample16),
                         c(rep(0, 22), round(breedingRelN16[[3]][1]*nSample16),
                           rep(0, 4), round(breedingRelN16[[3]][2]*nSample16),
                           rep(0, 44), round(breedingRelN16[[3]][3]*nSample16),
                           rep(0, 4), round(breedingRelN16[[3]][4]*nSample16),
                           rep(0, 22)),
                         round(breedingRelN16[[1]]*nSample16)[100:1],
                         c(rep(0, 22), round(breedingRelN16[[3]][1]*nSample16),
                           rep(0, 4), round(breedingRelN16[[3]][2]*nSample16),
                           rep(0, 44), round(breedingRelN16[[3]][3]*nSample16),
                           rep(0, 4), round(breedingRelN16[[3]][4]*nSample16),
                           rep(0, 22))[100:1])

#sapply(sampleBreeding16, sum)

# Run sampling regimes
scenarioToSampleMap16 <- c(1, 1, 2, 3, 3, 4)
animalLoc16 <- vector("list", nScenarios16)
results16 <- vector("list", nScenarios16)
compare16 <- data.frame(Scenario = c("True",
                                     "Base",
                                     "Breeding4",
                                     "CentroidSampleBreeding4",
                                     "BiasedSample",
                                     "BiasedSampleBreeding4",
                                     "BiasedCentroidSampleBreeding4"),
                        MC = c(MC16, rep(NA, nScenarios16)),
                        MCA = c(MC16, rep(NA, nScenarios16)),
                        MCss1 = c(MC16, rep(NA, nScenarios16)),
                        MCss2 = c(MC16, rep(NA, nScenarios16)),
                        Mantel = c(MC16, rep(NA, nScenarios16)))

compare16.array <- array(NA, c(nSims16, nScenarios16, 5),
                         dimnames = list(1:nSims16,
                                         c("Base", "Breeding4",
                                           "CentroidSampleBreeding4",
                                           "BiasedSample",
                                           "BiasedSampleBreeding4",
                                           "BiasedCentroidSampleBreeding4"),
                                         c("MC", "MCA", "MCss1", "MCss2", "Mantel")))

compare16.array2 <- array(NA, c(nSims16, 5),
                          dimnames = list(1:nSims16,
                                          c("MC", "MCA", "MCss1", "MCss2", "Mantel")))
set.seed(80)
for (sim in 1:nSims16) {
  cat("Simulation", sim, "of", nSims16, '\n')

  sim16 <- lapply(sampleBreeding16, simMove, breedingDist = breedDist16[[1]],
                  winteringDist=nonbreedDist16[[1]], psi=psi16, nYears=nYears,
                  nMonths=nMonths)

  for (i in nScenarios16) {
    cat("\tScenario", i, "\n")
    animalLoc16[[i]] <- changeLocations(sim16[[scenarioToSampleMap16[i]]]$animalLoc,
                                        breedingSiteTrans16[[i]], 1:nWintering16[i])

    results16[[i]] <- calcPsiMC(originDist = breedDist16[[i]],
                                targetDist = nonbreedDist16[[i]],
                                locations = animalLoc16[[i]],
                                originRelAbund = breedingRelN16[[i]],
                                verbose = FALSE)

    compare16.array2[sim, 'MC'] <- results16[[i]]$MC

    compare16.array2[sim, 'MCA'] <- calcMC(originDist = breedDist16[[i]],
                                           targetDist = nonbreedDist16[[i]],
                                           psi = results16[[i]]$psi,
                                           originRelAbund = breedingN16[[i]]/sum(breedingN16[[i]]))

    compare16.array2[sim, 'MCss1'] <- calcMC(originDist = breedDist16[[i]],
                                             targetDist = nonbreedDist16[[i]],
                                             psi = results16[[i]]$psi,
                                             originRelAbund = table(animalLoc16[[i]][ , 1, 1, 1])/
                                               sum(table(animalLoc16[[i]][ , 1, 1, 1])))

    compare16.array2[sim, 'MCss2'] <- calcMC(originDist = breedDist16[[i]],
                                             targetDist = nonbreedDist16[[i]],
                                             psi = results16[[i]]$psi,
                                             originRelAbund = (breedingRelN16[[i]] * nSample16)/
                                               sum((breedingRelN16[[i]] * nSample16))) # pmax(breedingRelN16[[i]] * nSample16, 1))
    compare16.array2[sim, 'Mantel'] <- calcStrengthInd(breedDist16[[1]],
                                                       nonbreedDist16[[1]],
                                                       sim16[[scenarioToSampleMap16[i]]]$animalLoc,
                                                       resamp=0)$correlation
  }
}

compare16$MC[1:nScenarios16 + 1] <- apply(compare16.array2[,,'MC'], 2, mean, na.mean = TRUE)
compare16$MCA[1:nScenarios16 + 1] <- apply(compare16.array2[,,'MCA'], 2, mean, na.mean = TRUE)
compare16$MCss1[1:nScenarios16 + 1] <- apply(compare16.array2[,,'MCss1'], 2, mean, na.mean = TRUE)
compare16$MCss2[1:nScenarios16 + 1] <- apply(compare16.array2[,,'MCss2'], 2, mean, na.mean = TRUE)
compare16$Mantel[1:nScenarios16 + 1] <- apply(compare16.array2[,,'Mantel'], 2, mean, na.mean = TRUE)
compare16 <- transform(compare16, MC.diff=MC - MC[1], MCA.diff=MCA - MCA[1],
                       MCss1.diff = MCss1 - MCss1[1],
                       MCss2.diff = MCss2 - MCss2[1],
                       Mantel.diff=Mantel - Mantel[1])

compare16a <- as.matrix(compare16[1 + 1:nScenarios16, c('MC', 'MC.diff', "Mantel", "Mantel.diff")])
rownames(compare16a) <- compare16$Scenario[1 + 1:nScenarios16]
