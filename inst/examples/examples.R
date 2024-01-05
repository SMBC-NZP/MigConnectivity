library(MigConnectivity)
set.seed(75)

###############################################################################
# Utility functions for use in simulations
###############################################################################
# Calculates probability matrix based on exponential decline with distance
mlogitMat <- function(slope, dist) {
  preMat <- exp(-slope/mean(dist)*dist)
  diag(preMat) <- 0
  nr <- nrow(dist)
  nc <- ncol(dist)
  outMat <- matrix(0, nr, nc)
  for (b in 1:nr) {
    outMat[b,] <- preMat[b,]/(1+sum(preMat[b, ]))
    outMat[b,b] <- 1 - sum(outMat[b, ])
  }
  return(outMat)
}

# Crude optimizable function for developing MC pattern based on MC strength
mlogitMC <- function(slope, MC.in, origin.dist, target.dist, origin.rel.abund) {
  nBreeding <- nrow(origin.dist)
  nWintering <- nrow(target.dist)
  psi <- mlogitMat(slope, origin.dist)
  if (any(psi<0))
    return(5*slope^2)
  MC <- calcMC(origin.dist, target.dist, psi, origin.rel.abund)
  return((MC.in - MC)^2)
}

# rho function for individuals
calcStrengthInd <- function(originDist, targetDist, locations, resamp=1000, verbose = 0) {
  nInd <- dim(locations)[1]
  originDist2 <- targetDist2 <- matrix(0, nInd, nInd)
  for (i in 1:(nInd-1)) {
    for (j in (i+1):nInd) {
      originDist2[i,j] <- originDist2[j,i] <- originDist[locations[i,1,1,1], locations[j,1,1,1]]
      targetDist2[i,j] <- targetDist2[j,i] <- targetDist[locations[i,2,1,1], locations[j,2,1,1]]
    }
  }
  return(ncf::mantel.test(originDist2, targetDist2, resamp=resamp, quiet = !verbose))
}

###############################################################################
# rho function for individuals
###############################################################################
calcStrengthInd <- function(originDist, targetDist, locations, resamp=1000,
                            verbose = 0) {
  nInd <- dim(locations)[1]
  originDist2 <- targetDist2 <- matrix(0, nInd, nInd)
  for (i in 1:(nInd-1)) {
    for (j in (i+1):nInd) {
      originDist2[i,j] <- originDist2[j,i] <- originDist[locations[i,1,1,1],
                                                         locations[j,1,1,1]]
      targetDist2[i,j] <- targetDist2[j,i] <- targetDist[locations[i,2,1,1],
                                                         locations[j,2,1,1]]
    }
  }
  return(ncf::mantel.test(originDist2, targetDist2, resamp=resamp, quiet = !verbose))
}

###############################################################################
# Simple approach to estimate psi matrix and MC from simulated (or real) data
# (doesn't include uncertainty)
###############################################################################
calcPsiMC <- function(originDist, targetDist, originRelAbund, locations,
                      verbose=FALSE) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  psiMat <- matrix(0, nOrigin, nTarget)
  nInd <- dim(locations)[1]
  nYears <- dim(locations)[3]
  nMonths <- dim(locations)[4]
  for (i in 1:nInd) {
    if (i %% 1000 == 0 && verbose) #
      cat("Individual", i, "of", nInd, "\n")
    originMat <- locations[i,1,,]
    targetMat <- locations[i,2,,]
    bIndices <- which(!is.na(originMat))
    wIndices <- which(!is.na(targetMat))
    if (length(bIndices) && length(wIndices))
      for (bi in bIndices)
        for (wi in wIndices)
          psiMat[originMat[bi], targetMat[wi]] <- psiMat[originMat[bi],
                                                         targetMat[wi]] + 1
  }
  psiMat <- apply(psiMat, 2, "/", rowSums(psiMat))
  MC <- calcMC(originDist, targetDist, psiMat, originRelAbund)
  return(list(psi=psiMat, MC=MC))
}


###############################################################################
# Parameters for simulations
###############################################################################
\donttest{
nSeasons <- 2
nYears <- 10
nMonths <- 4 # Each season

nBreeding <- 100
nWintering <- 100
breedingPos <- matrix(c(rep(seq(-99,-81,2), each=sqrt(nBreeding)),
                        rep(seq(49,31,-2), sqrt(nBreeding))), nBreeding, 2)
winteringPos <- matrix(c(rep(seq(-79,-61,2), each=sqrt(nWintering)),
                         rep(seq(9,-9,-2), sqrt(nWintering))), nWintering, 2)
head(breedingPos)
tail(breedingPos)
head(winteringPos)
tail(winteringPos)

breedDist <- distFromPos(breedingPos, 'ellipsoid')
nonbreedDist <- distFromPos(winteringPos, 'ellipsoid')
breedDist[1:12, 1:12]
breedDist[1:12, c(1,91,100)]

# Breeding Abundance
breedingN <- rep(500, nBreeding)
breedingRelN <- breedingN/sum(breedingN)

# Set up psi matrix
o <- optimize(mlogitMC, MC.in = 0.25, origin.dist = breedDist,
              target.dist = nonbreedDist, origin.rel.abund = breedingRelN,
              interval = c(0, 10), tol = .Machine$double.eps^0.5)
o
slope <- o$minimum
psi <- mlogitMat(slope, breedDist)
round(psi[1:12, 1:12],5)
rowSums(psi)

# Baseline strength of migratory connectivity
MC <- calcMC(breedDist, nonbreedDist, psi, breedingRelN)
MC

# Simulation
sim <- simMove(breedingN, breedDist, nonbreedDist, psi, nYears, nMonths)

###############################
# Sampling regime 1 of 3
# Researchers divide populations differently than reality
# Delineation of seasonal ranges into regions
#I) Breeding range divided along equal longitudinal breaks into ten regions
#II) Non-breeding range divided along equal longitudinal breaks into ten regions
#III) Breeding and non-breeding ranges divided along equal longitudinal breaks into ten regions
#IV) Breeding range divided along the longitudinal and latitudinal midpoint into four regions
#V) Non-breeding range divided along the longitudinal and latitudinal midpoint into four regions
#VI) Breeding range divided along the longitudinal and latitudinal midpoint into four regions and
#non-breeding range divided along equal longitudinal breaks into ten regions
#VII) Breeding range divided along equal longitudinal breaks into ten regions and
#non-breeding range divided along the longitudinal and latitudinal midpoint into four regions
###############################
#Run functions and parameters above first

set.seed(75)

#each element is for a scenario (see above 1-8)
breedingSiteTrans14 <- list(1:nBreeding, rep(1:10, each=10), 1:nBreeding, rep(1:10, each=10),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)), 1:nBreeding,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)), rep(1:nBreeding, each=10))
winteringSiteTrans14 <- list(1:nWintering, 1:nWintering, rep(1:10, each=10), rep(1:10, each=10),
                             1:nWintering, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)), rep(1:10, each=10),
                             c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))
breedingSiteTrans14
lapply(breedingSiteTrans14, matrix, nrow=10, ncol=10)
lapply(winteringSiteTrans14, matrix, nrow=10, ncol=10)
breedingPos14 <- list(breedingPos, rowsum(breedingPos, rep(1:10, each=10))/10, #positions of the human defined populations
                      breedingPos, rowsum(breedingPos, rep(1:10, each=10))/10,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))/25, breedingPos,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))/25,
                      rowsum(breedingPos, rep(1:10, each=10))/10)
breedingPos14
winteringPos14 <- list(winteringPos, winteringPos,
                       rowsum(winteringPos, rep(1:10, each=10))/10,
                       rowsum(winteringPos, rep(1:10, each=10))/10, winteringPos,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))/25,
                       rowsum(winteringPos, rep(1:10, each=10))/10,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))/25)
winteringPos14

breedDist14 <- lapply(breedingPos14, distFromPos, surface = 'ellipsoid')
breedDist14
lapply(breedDist14, max)
nonbreedDist14 <- lapply(winteringPos14, distFromPos, surface = 'ellipsoid')
nonbreedDist14
lapply(nonbreedDist14, max)
nBreeding14 <- c(100, 10, 100, 10, 4, 100, 4, 10)
nWintering14 <- c(100, 100, 10, 10, 100, 4, 10, 4)
breedingRelN14 <- lapply(nBreeding14, function(x) rep(1/x, x))
breedingRelN14
nSample14 <- 1000 # Total number sampled per simulation

# for the baseline use the simulation from above, sim
animalLoc14base <- sim$animalLoc
#transferring the simulated bird locations from the true populations to the researcher defined populations
changeLocations <- function(animalLoc, breedingSiteTrans, winteringSiteTrans) {
  animalLoc[,1,,] <- breedingSiteTrans[animalLoc[,1,,]]
  animalLoc[,2,,] <- winteringSiteTrans[animalLoc[,2,,]]
  return(animalLoc)
}

nScenarios14 <- length(breedingSiteTrans14)
animalLoc14 <- vector("list", nScenarios14) #making an empty list to fill
for (i in 1:nScenarios14)
  animalLoc14[[i]] <- changeLocations(animalLoc14base, breedingSiteTrans14[[i]], winteringSiteTrans14[[i]])
animalLoc14[[2]][101:110,1,1,1]
animalLoc14[[6]][101:110,1,1,1]
animalLoc14[[6]][101:110,2,1,1]
results14 <- vector("list", nScenarios14)
compare14 <- data.frame(Scenario = c("True", "Base", "Breeding10", "Wintering10",
                                     "Breeding10Wintering10", "Breeding4",
                                     "Wintering4", "Breeding4Wintering10",
                                     "Breeding10Wintering4"),
                        MC = c(MC, rep(NA, nScenarios14)))
for (i in 1:nScenarios14) {
  cat("\nScenario", i, "\n")
  results14[[i]] <- calcPsiMC(breedDist14[[i]], nonbreedDist14[[i]],
                              breedingRelN14[[i]], animalLoc14[[i]], TRUE)
  compare14$MC[i+1] <- results14[[i]]$MC
}
compare14 <- transform(compare14, diff=MC - MC[1], prop=MC/MC[1])
compare14


###############################################################################
# Sampling regime 2  of 3
# Researchers divide populations differently than reality PLUS
# Different distributions of sampling animals across breeding range PLUS
# Sample sizes don't always match relative abundances PLUS
# Compare our approach and simple Mantel approach
#   1. Base (10 years, uneven abundances but matches sampling)
#   2. Breeding pops divided into 4 squares, sample across breeding range
#   3. Breeding pops divided into 4 squares, sample at centroid of each square
#   4. Sampling high in low abundance populations plus base
#   5. Scenarios 2 plus 4
#   6. Scenarios 3 plus 4
#May want to create another set because we have varied these two things together
###############################################################################

set.seed(75)

# Transfer between true populations and researcher defined ones (only for
# breeding, as not messing with winter populations here)
breedingSiteTrans15 <- list(1:100, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)), 1:100,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))
breedingSiteTrans15
lapply(breedingSiteTrans15, matrix, nrow=10, ncol=10)
nScenarios15 <- length(breedingSiteTrans15)
nSims15 <- 100
# Basing positions of researcher defined breeding populations on above
breedingPos15 <- list(breedingPos, breedingPos14[[5]],
                      breedingPos14[[5]], breedingPos, breedingPos14[[5]], breedingPos14[[5]])
breedingPos15
winteringPos15 <- rep(list(winteringPos), nScenarios15)
winteringPos15
breedDist15 <- lapply(breedingPos15, distFromPos)
breedDist15[[1]][1,]
nonbreedDist15 <- lapply(winteringPos15, distFromPos)
nBreeding15 <- rep(c(100, 4, 4), 2)
nBreeding15
nWintering15 <- rep(100, nScenarios15)
# Highest abundance in lower right corner, lowest in top left
# Making symmetrical
breedingN15base <- rep(NA, 100)
for (i in 1:10) #row
  for (j in 1:10)  #column
    breedingN15base[i+10*(j-1)] <- 500 + 850*i*j
breedingN15base
matrix(breedingN15base, 10, 10)
sum(breedingN15base)
breedingN15 <- lapply(breedingSiteTrans15, rowsum, x=breedingN15base) # For researcher defined populations
breedingRelN15 <- lapply(breedingN15, "/", sum(breedingN15base))
nSample15 <- 1000 # Total number sampled per simulation
# Four sampling regimes/simulations, because repeat a couple of them in different scenarios
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
lapply(sampleBreeding15, matrix, nrow=10, ncol=10)
lapply(sampleBreeding15, sum)

# Set up psi matrix
o15 <- optimize(mlogitMC, MC.in = 0.25, origin.dist = breedDist15[[1]],
                target.dist = nonbreedDist15[[1]],
                origin.rel.abund = breedingRelN15[[1]], interval = c(0,10),
                tol = .Machine$double.eps^0.5)
o15
slope15 <- o15$minimum
psi15 <- mlogitMat(slope15, breedDist15[[1]])
round(psi15[1:12, 1:12],5)
rowSums(psi15)
# Baseline strength of migratory connectivity
MC15 <- calcMC(breedDist15[[1]], nonbreedDist15[[1]], psi15, breedingRelN15[[1]])
MC15
# Run sampling regimes
scenarioToSampleMap15 <- c(1, 1, 2, 3, 3, 4)
animalLoc15 <- vector("list", nScenarios15)

results15 <- vector("list", nScenarios15)
compare15 <- data.frame(Scenario = c("True", "Base", "Breeding4",
                                     "CentroidSampleBreeding4", "BiasedSample", "BiasedSampleBreeding4",
                                     "BiasedCentroidSampleBreeding4"),
                        MC = c(MC15, rep(NA, nScenarios15)), Mantel = c(MC15, rep(NA, nScenarios15)))
compare15.array <- array(NA, c(nSims15, nScenarios15, 2),
                         dimnames = list(1:nSims15,
                                         c("Base", "Breeding4", "CentroidSampleBreeding4",
                                           "BiasedSample", "BiasedSampleBreeding4",
                                           "BiasedCentroidSampleBreeding4"),
                                         c("MC", "Mantel")))
for (sim in 1:nSims15) {
  cat("Simulation", sim, "of", nSims15, '\n')
  sim15 <- lapply(sampleBreeding15, simMove, breedingDist = breedDist15[[1]],
                  winteringDist=nonbreedDist15[[1]], psi=psi15, nYears=nYears,
                  nMonths=nMonths)
  for (i in 1:nScenarios15) {
    cat("\tScenario", i, "\n")
    animalLoc15[[i]] <- changeLocations(sim15[[scenarioToSampleMap15[i]]]$animalLoc,
                                      breedingSiteTrans15[[i]], 1:nWintering15[i])
    results15[[i]] <- calcPsiMC(breedDist15[[i]], nonbreedDist15[[i]],
                                   breedingRelN15[[i]], animalLoc15[[i]], FALSE)
    compare15.array[sim, i, 'MC'] <- results15[[i]]$MC
    compare15.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist15[[1]],
                                                         nonbreedDist15[[1]],
                                                         sim15[[scenarioToSampleMap15[i]]]$animalLoc,
                                                         resamp=0)$correlation
  }
}

compare15$MC[1:nScenarios15 + 1] <- apply(compare15.array[,,'MC'], 2, mean)
compare15$Mantel[1:nScenarios15 + 1] <- apply(compare15.array[,,'Mantel'], 2,
                                              mean)
compare15 <- transform(compare15, MC.diff=MC - MC[1],
                       Mantel.diff=Mantel - Mantel[1],
                       MC.prop=MC/MC[1], Mantel.prop=Mantel/Mantel[1])
compare15
compare15a <- as.matrix(compare15[2:7,c(2,4,3,5)])
rownames(compare15a) <- compare15$Scenario[2:7]
round(compare15a, 3)
round(compare15a, 2)


###############################################################################
# Sampling regime 3 of 3
# Researchers divide populations differently than reality (simulations) PLUS
# Different distributions of sampled animals across breeding range PLUS
# Sample sizes don't always match relative abundances PLUS
# Compare our approach and simple Mantel approach PLUS
# MC not same across subsections of range
#   1. Base (uneven MC (0.15 for NW breeding, 0.3 for SW, 0.45 for NE, and 0.6
#     for SE), uneven abundances (lowest in NW, highest in SE), sampling
#     proportional to abundance
#   2. Breeding pops divided into 4 squares, sample across breeding range
#   3. Breeding pops divided into 4 squares, sample at centroid of each square
#   4. Sampling high in low abundance populations
#   5. Scenarios 2 plus 4
#   6. Scenarios 3 plus 4
###############################################################################

set.seed(75)


# Transfer between true populations and researcher defined ones
# (only for breeding, as not messing with winter populations here)
breedingSiteTrans16 <- list(1:100, c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)), 1:100,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))

lapply(breedingSiteTrans16, matrix, nrow=10, ncol=10)
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
lapply(breedingRelN16, sum)

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
                   origin.rel.abund = breedingRelN16[[1]],
                   interval=c(0,10), tol=.Machine$double.eps^0.5)
  slope16a <- o16a$minimum
  cat(" slope", slope16a, "\n")
  psi16a <- mlogitMat(slope16a, breedDist16[[1]])
  # Then use the rows of that psi matrix only for the one breeding quadrant
  rows <- 50*(i %/% 3) + rep(1:5, 5) + rep(seq(0, 40, 10), each=5) + ((i-1) %% 2) * 5
  psi16[rows, ] <- psi16a[rows, ]
}
round(psi16[1:12, 1:12],5)
rowSums(psi16)

# Baseline strength of migratory connectivity
MC16 <- calcMC(breedDist16[[1]], nonbreedDist16[[1]], psi16, breedingRelN16[[1]])
MC16

# Set up sampling regimes (different number than number of scenarios)
nSample16 <- 1000
sampleBreeding16 <- sampleBreeding15
lapply(sampleBreeding16, matrix, nrow=10, ncol=10)
lapply(sampleBreeding16, sum)



# Run sampling regimes
scenarioToSampleMap16 <- c(1, 1, 2, 3, 3, 4)
animalLoc16 <- vector("list", nScenarios16)
results16 <- vector("list", nScenarios16)
compare16 <- data.frame(Scenario = c("True", "Base", "Breeding4",
                                     "CentroidSampleBreeding4", "BiasedSample",
                                     "BiasedSampleBreeding4",
                                     "BiasedCentroidSampleBreeding4"),
                        MC = c(MC16, rep(NA, nScenarios16)),
                        Mantel = c(MC16, rep(NA, nScenarios16)))
compare16.array <- array(NA, c(nSims16, nScenarios16, 2),
                         dimnames = list(1:nSims16,
                                         c("Base", "Breeding4",
                                           "CentroidSampleBreeding4",
                                           "BiasedSample",
                                           "BiasedSampleBreeding4",
                                           "BiasedCentroidSampleBreeding4"),
                                         c("MC", "Mantel")))
for (sim in 1:nSims15) {
  cat("Simulation", sim, "of", nSims15, '\n')
  sim15 <- lapply(sampleBreeding15, simMove, breedingDist = breedDist15[[1]],
                  winteringDist=nonbreedDist15[[1]], psi=psi15, nYears=nYears,
                  nMonths=nMonths)
  for (i in 1:nScenarios15) {
    cat("\tScenario", i, "\n")
    animalLoc15[[i]] <- changeLocations(sim15[[scenarioToSampleMap15[i]]]$animalLoc,
                                        breedingSiteTrans15[[i]], 1:nWintering15[i])
    results15[[i]] <- calcPsiMC(breedDist15[[i]], nonbreedDist15[[i]],
                                breedingRelN15[[i]], animalLoc15[[i]], FALSE)
    compare15.array[sim, i, 'MC'] <- results15[[i]]$MC
    compare15.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist15[[1]],
                                                         nonbreedDist15[[1]],
                                                         sim15[[scenarioToSampleMap15[i]]]$animalLoc,
                                                         resamp=0)$correlation
  }
}
for (sim in 1:nSims16) {
  cat("Simulation", sim, "of", nSims16, '\n')
  sim16 <- lapply(sampleBreeding16, simMove, breedingDist = breedDist16[[1]],
                  winteringDist=nonbreedDist16[[1]], psi=psi16, nYears=nYears,
                  nMonths=nMonths)
  for (i in 1:nScenarios16) {
    cat("\tScenario", i, "\n")
    animalLoc16[[i]] <- changeLocations(sim16[[scenarioToSampleMap16[i]]]$animalLoc,
                                      breedingSiteTrans16[[i]], 1:nWintering16[i])
    results16[[i]] <- calcPsiMC(breedDist16[[i]], nonbreedDist16[[i]],
                                breedingRelN16[[i]], animalLoc16[[i]], FALSE)
    compare16.array[sim, i, 'MC'] <- results16[[i]]$MC
    compare16.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist16[[1]],
                                                         nonbreedDist16[[1]],
                                                         sim16[[scenarioToSampleMap16[i]]]$animalLoc,
                                                         resamp=0)$correlation
  }
}
compare16$MC[1:nScenarios16 + 1] <- apply(compare16.array[,,'MC'], 2, mean)
compare16$Mantel[1:nScenarios16 + 1] <- apply(compare16.array[,,'Mantel'], 2, mean)
compare16 <- transform(compare16, MC.diff=MC - MC[1], Mantel.diff=Mantel - Mantel[1],
                       MC.prop=MC/MC[1], Mantel.prop=Mantel/Mantel[1])
compare16
compare16a <- as.matrix(compare16[2:7,c(2,4,3,5)])
rownames(compare16a) <- compare16$Scenario[2:7]
round(compare16a, 3)
round(compare16a, 2)
}
