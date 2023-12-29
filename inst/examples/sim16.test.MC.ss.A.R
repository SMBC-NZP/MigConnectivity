calcPsiMC0 <- function(originDist, targetDist, originRelAbund, locations,
                      years = 1, months = 1, verbose=F) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  psiMat <- matrix(0, nOrigin, nTarget)
  nInd <- dim(locations)[1]
  nYears <- dim(locations)[3]
  nMonths <- dim(locations)[4]
  for (i in 1:nInd) {
    if (i %% 1000 == 0 && verbose) #
      cat("Individual", i, "of", nInd, "\n")
    originMat <- locations[i, 1, years, months]
    targetMat <- locations[i, 2, years, months]
    bIndices <- which(!is.na(originMat))
    wIndices <- which(!is.na(targetMat))
    if (length(bIndices) && length(wIndices))
      for (bi in bIndices)
        for (wi in wIndices)
          psiMat[originMat[bi], targetMat[wi]] <- psiMat[originMat[bi], targetMat[wi]] + 1
  }
  psiMat <- apply(psiMat, 2, "/", rowSums(psiMat))
  MC <- calcMC(originDist, targetDist, psi = psiMat,
               originRelAbund = originRelAbund)
  return(list(psi=psiMat, MC=MC))
}

calcPsiMC <- function(originDist, targetDist, originRelAbund, locations,
                      years = 1, months = 1, verbose=F) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  psiMat <- matrix(0, nOrigin, nTarget)
  nInd <- dim(locations)[1]
  nYears <- dim(locations)[3]
  nMonths <- dim(locations)[4]
  for (i in 1:nInd) {
    if (i %% 1000 == 0 && verbose) #
      cat("Individual", i, "of", nInd, "\n")
    originMat <- locations[i, 1, years, months]
    targetMat <- locations[i, 2, years, months]
    bIndices <- which(!is.na(originMat))
    wIndices <- which(!is.na(targetMat))
    if (length(bIndices) && length(wIndices))
      for (bi in bIndices)
        for (wi in wIndices)
          psiMat[originMat[bi], targetMat[wi]] <- psiMat[originMat[bi], targetMat[wi]] + 1
  }
  psiMat <- apply(psiMat, 2, "/", rowSums(psiMat))
  MC <- calcMC(originDist, targetDist, psi = psiMat,
               originRelAbund = originRelAbund, sampleSize = nInd)
  return(list(psi=psiMat, MC=MC))
}

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

sapply(sampleBreeding16, sum)

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
                                     "BiasedCentroidSampleBreeding4"))

compare16.array <- array(NA, c(nSims16, nScenarios16, 5),
                         dimnames = list(1:nSims16,
                                         c("Base", "Breeding4",
                                           "CentroidSampleBreeding4",
                                           "BiasedSample",
                                           "BiasedSampleBreeding4",
                                           "BiasedCentroidSampleBreeding4"),
                                         c("MCR", "MCA", "MC", "MCss", "Mantel")))

set.seed(80)
for (sim in 1:nSims16) {
  cat("Simulation", sim, "of", nSims16, '\n')

  sim16 <- lapply(sampleBreeding16, simMove, breedingDist = breedDist16[[1]],
                  winteringDist=nonbreedDist16[[1]], psi=psi16, nYears=nYears,
                  nMonths=nMonths)

  for (i in c(2, 3, 5, 6)) {
    cat("\tScenario", i, "\n")
    animalLoc16[[i]] <- changeLocations(sim16[[scenarioToSampleMap16[i]]]$animalLoc,
                                        breedingSiteTrans16[[i]], 1:nWintering16[i])

    results16[[i]] <- calcPsiMC0(originDist = breedDist16[[i]],
                                 targetDist = nonbreedDist16[[i]],
                                 locations = animalLoc16[[i]],
                                 originRelAbund = breedingRelN16[[i]],
                                 verbose = FALSE)

    compare16.array[sim, i, 'MCR'] <- results16[[i]]$MC

    compare16.array[sim, i, 'MC'] <- calcMC(originDist = breedDist16[[i]],
                                          targetDist = nonbreedDist16[[i]],
                                          psi = results16[[i]]$psi,
                                          originRelAbund = breedingRelN16[[i]],
                                          sampleSize = dim(animalLoc16[[i]])[1])

    compare16.array[sim, i, 'MCA'] <- calcMC(originDist = breedDist16[[i]],
                                          targetDist = nonbreedDist16[[i]],
                                          psi = results16[[i]]$psi,
                                          originRelAbund = breedingRelN16[[i]],
                                          sampleSize = sum(breedingN16[[i]]))

    compare16.array[sim, i, 'MCss'] <- calcMC(originDist = breedDist16[[i]],
                                           targetDist = nonbreedDist16[[i]],
                                           psi = results16[[i]]$psi,
                                           originRelAbund = table(animalLoc16[[i]][,1,1,1])/dim(animalLoc16[[i]])[1],
                                           sampleSize = dim(animalLoc16[[i]])[1])
    compare16.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist16[[1]],
                                                       nonbreedDist16[[1]],
                                                       sim16[[scenarioToSampleMap16[i]]]$animalLoc,
                                                       resamp=0)$correlation
  }
}

compare16$MCR <- c(MC16, apply(compare16.array[,,'MCR'], 2, mean, na.mean = TRUE))
compare16$MCA <- c(MC16, apply(compare16.array[,,'MCA'], 2, mean, na.mean = TRUE))
compare16$MC <- c(MC16, apply(compare16.array[,,'MC'], 2, mean, na.mean = TRUE))
compare16$MCss <- c(MC16, apply(compare16.array[,,'MCss'], 2, mean, na.mean = TRUE))
compare16$Mantel <- c(MC16, apply(compare16.array[,,'Mantel'], 2, mean, na.mean = TRUE))
compare16 <- transform(compare16, MC.diff=MC - MC[1], MCA.diff=MCA - MCA[1],
                       MCss.diff = MCss - MCss[1],
                       MCR.diff = MCR - MCR[1],
                       Mantel.diff=Mantel - Mantel[1])





set.seed(75)

# Transfer between true populations and researcher defined ones (only for
# breeding, as not messing with winter populations here)

breedingSiteTrans15 <- list(1:100,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            1:100,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))

nScenarios15 <- length(breedingSiteTrans15)

nSims15 <- 100

# Basing positions of researcher defined breeding populations on above
breedingPos15 <- list(breedingPos,
                      breedingPos14[[5]],
                      breedingPos14[[5]],
                      breedingPos,
                      breedingPos14[[5]],
                      breedingPos14[[5]])

winteringPos15 <- rep(list(winteringPos), nScenarios15)

breedDist15 <- lapply(breedingPos15, distFromPos)

nonbreedDist15 <- lapply(winteringPos15, distFromPos)

nBreeding15 <- rep(c(100, 4, 4), 2)

nWintering15 <- rep(100, nScenarios15)

# Highest abundance in lower right corner, lowest in top left
# Making symmetrical

breedingN15base <- rep(NA, 100)
for (i in 1:10) #row
  for (j in 1:10)  #column
    breedingN15base[i+10*(j-1)] <- 500 + 850*i*j

sum(breedingN15base)

# For researcher defined populations
breedingN15 <- lapply(breedingSiteTrans15, rowsum, x=breedingN15base)

breedingRelN15 <- lapply(breedingN15, "/", sum(breedingN15base))

nSample15 <- 100 # Total number sampled per simulation

# Number sampled per natural population
sampleBreeding15 <- list(round(breedingRelN15[[1]]*nSample15),
                         c(rep(0, 22), round(breedingRelN15[[3]][1]*nSample15),
                           rep(0, 4), round(breedingRelN15[[3]][2]*nSample15),
                           rep(0, 44), round(breedingRelN15[[3]][3]*nSample15),
                           rep(0, 4), round(breedingRelN15[[3]][4]*nSample15),
                           rep(0, 22)),
                         round(breedingRelN15[[1]]*nSample15)[100:1],
                         c(rep(0, 22), round(breedingRelN15[[3]][1]*nSample15),
                           rep(0, 4), round(breedingRelN15[[3]][2]*nSample15),
                           rep(0, 44), round(breedingRelN15[[3]][3]*nSample15),
                           rep(0, 4), round(breedingRelN15[[3]][4]*nSample15),
                           rep(0, 22))[100:1])

# Set up psi matrix
o15 <- optimize(mlogitMC, MC.in = 0.25, origin.dist = breedDist15[[1]],
                target.dist = nonbreedDist15[[1]],
                origin.abund = breedingN15[[1]]/sum(breedingN15[[1]]),
                sample.size = sum(breedingN15[[1]]),
                interval = c(0,10),
                tol = .Machine$double.eps^0.5)

slope15 <- o15$minimum

psi15 <- mlogitMat(slope15, breedDist15[[1]])

# Baseline strength of migratory connectivity
MC15 <- calcMC(originDist = breedDist15[[1]],
               targetDist = nonbreedDist15[[1]],
               psi = psi15,
               originRelAbund = breedingN15[[1]]/sum(breedingN15[[1]]),
               sampleSize = sum(breedingN15[[1]]))

# Run sampling regimes
scenarioToSampleMap15 <- c(1, 1, 2, 3, 3, 4)

animalLoc15 <- vector("list", nScenarios15)

results15 <- vector("list", nScenarios15)

compare15 <- data.frame(Scenario = c("True",
                                     "Base",
                                     "Breeding4",
                                     "CentroidSampleBreeding4",
                                     "BiasedSample",
                                     "BiasedSampleBreeding4",
                                     "BiasedCentroidSampleBreeding4"))

compare15.array <- array(NA, c(nSims15, nScenarios15, 5),
                         dimnames = list(1:nSims15,
                                         c("Base", "Breeding4", "CentroidSampleBreeding4",
                                           "BiasedSample", "BiasedSampleBreeding4",
                                           "BiasedCentroidSampleBreeding4"),
                                         c("MCR", "MCA", "MC", "MCss", "Mantel")))


for (sim in 1:nSims15) {
  cat("Simulation", sim, "of", nSims15, '\n')

  sim15 <- lapply(sampleBreeding15, simMove, breedingDist = breedDist15[[1]],
                  winteringDist=nonbreedDist15[[1]], psi=psi15, nYears=nYears,
                  nMonths=nMonths)

  for (i in c(2, 3, 5, 6)) {

    cat("\tScenario", i, "\n")
    animalLoc15[[i]] <- changeLocations(animalLoc = sim15[[scenarioToSampleMap15[i]]]$animalLoc,
                                        breedingSiteTrans = breedingSiteTrans15[[i]],
                                        winteringSiteTrans = 1:nWintering15[i])

    results15[[i]] <- calcPsiMC(originDist = breedDist15[[i]],
                                targetDist = nonbreedDist15[[i]],
                                originRelAbund = breedingRelN15[[i]],
                                locations = animalLoc15[[i]],
                                verbose = F)

    compare15.array[sim, i, 'MC'] <- results15[[i]]$MC
    compare15.array[sim, i, 'MCR'] <- calcMC(originDist = breedDist15[[i]],
                                             targetDist = nonbreedDist15[[i]],
                                             psi = results15[[i]]$psi,
                                             originRelAbund = breedingRelN15[[i]])

    compare15.array[sim, i, 'MCA'] <- calcMC(originDist = breedDist15[[i]],
                                             targetDist = nonbreedDist15[[i]],
                                             psi = results15[[i]]$psi,
                                             originRelAbund = breedingRelN15[[i]],
                                             sampleSize = sum(breedingN15[[i]]))

    compare15.array[sim, i, 'MCss'] <- calcMC(originDist = breedDist15[[i]],
                                              targetDist = nonbreedDist15[[i]],
                                              psi = results15[[i]]$psi,
                                              originRelAbund = table(animalLoc15[[i]][,1,1,1])/dim(animalLoc15[[i]])[1],
                                              sampleSize = dim(animalLoc15[[i]])[1])
    compare15.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist15[[1]],
                                                         nonbreedDist15[[1]],
                                                         sim15[[scenarioToSampleMap15[i]]]$animalLoc,
                                                         resamp=0)$correlation
  }
}

compare15$MCR <- c(MC15, apply(compare15.array[,,'MCR'], 2, mean, na.mean = TRUE))
compare15$MCA <- c(MC15, apply(compare15.array[,,'MCA'], 2, mean, na.mean = TRUE))
compare15$MC <- c(MC15, apply(compare15.array[,,'MC'], 2, mean, na.mean = TRUE))
compare15$MCss <- c(MC15, apply(compare15.array[,,'MCss'], 2, mean, na.mean = TRUE))
compare15$Mantel <- c(MC15, apply(compare15.array[,,'Mantel'], 2, mean, na.mean = TRUE))
compare15 <- transform(compare15, MC.diff=MC - MC[1], MCA.diff=MCA - MCA[1],
                       MCss.diff = MCss - MCss[1],
                       MCR.diff = MCR - MCR[1],
                       Mantel.diff=Mantel - Mantel[1])



