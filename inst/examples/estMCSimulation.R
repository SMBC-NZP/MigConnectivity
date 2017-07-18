set.seed(75)
nBreeding <- 100 # Number of populations
nWintering <- 100 # Number of populations
breedingPos <- matrix(c(rep(seq(-99,-81,2), each=sqrt(nBreeding)),
                        rep(seq(49,31,-2), sqrt(nBreeding))), nBreeding, 2)
winteringPos <- matrix(c(rep(seq(-79,-61,2), each=sqrt(nWintering)),
                         rep(seq(9,-9,-2), sqrt(nWintering))), nWintering, 2)

capLocs<-SpatialPoints(breedingPos, CRS(WGS84))
targLocs<-SpatialPoints(winteringPos, CRS(WGS84))
#capLocsM<-spTransform(capLocs, CRS(Lambert))
#targLocsM<-spTransform(targLocs, CRS(Lambert))

breedDist <- distFromPos(breedingPos, 'ellipsoid') # calculate distance between populations
nonbreedDist <- distFromPos(winteringPos, 'ellipsoid') # calculate distance between populations


scenarios14 <- c("Base",
                 "Breeding10",
                 "Wintering10",
                 "Breeding10Wintering10",
                 "Breeding4",
                 "Wintering4",
                 "Breeding4Wintering10",
                 "Breeding10Wintering4",
                 "Breeding4Wintering4",
                 "CentroidSampleBreeding10",
                 "CentroidSampleBreeding10Wintering10",
                 "CentroidSampleBreeding10Wintering4",
                 "CentroidSampleBreeding4",
                 "CentroidSampleBreeding4Wintering10",
                 "CentroidSampleBreeding4Wintering4")

# Each element is for a scenario (see above 1-8), transferring from natural
# breeding populations to defined ones
breedingSiteTrans14 <- list(1:nBreeding,
                            rep(1:10, each=10),
                            1:nBreeding,
                            rep(1:10, each=10),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            1:nBreeding,
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            rep(1:10, each=10),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            rep(1:10, each=10),
                            rep(1:10, each=10),
                            rep(1:10, each=10),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                            c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))

# Same for non-breeding populations
winteringSiteTrans14 <- list(1:nWintering,
                             1:nWintering,
                             rep(1:10, each=10),
                             rep(1:10, each=10),
                             1:nWintering,
                             c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                             rep(1:10, each=10),
                             c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                             c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                             1:nWintering,
                             rep(1:10, each=10),
                             c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)),
                             1:nWintering,
                             rep(1:10, each=10),
                             c(rep(1:2, 5, each=5), rep(3:4, 5, each=5)))

# Examine the transfers in matrix form
lapply(breedingSiteTrans14, matrix, nrow=10, ncol=10)
lapply(winteringSiteTrans14, matrix, nrow=10, ncol=10)

#positions of the human defined populations
breedingPos14 <- list(breedingPos,
                      rowsum(breedingPos, rep(1:10, each=10))/10,
                      breedingPos,
                      rowsum(breedingPos, rep(1:10, each=10))/10,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5),
                                            rep(3:4, 5, each=5)))/25,
                      breedingPos,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5),
                                            rep(3:4, 5, each=5)))/25,
                      rowsum(breedingPos, rep(1:10, each=10))/10,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5),
                                            rep(3:4, 5, each=5)))/25,
                      rowsum(breedingPos, rep(1:10, each=10))/10,
                      rowsum(breedingPos, rep(1:10, each=10))/10,
                      rowsum(breedingPos, rep(1:10, each=10))/10,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5),
                                            rep(3:4, 5, each=5)))/25,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5),
                                            rep(3:4, 5, each=5)))/25,
                      rowsum(breedingPos, c(rep(1:2, 5, each=5),
                                            rep(3:4, 5, each=5)))/25)

winteringPos14 <- list(winteringPos,
                       winteringPos,
                       rowsum(winteringPos, rep(1:10, each=10))/10,
                       rowsum(winteringPos, rep(1:10, each=10))/10,
                       winteringPos,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5),
                                              rep(3:4, 5, each=5)))/25,
                       rowsum(winteringPos,rep(1:10, each=10))/10,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5),
                                              rep(3:4, 5, each=5)))/25,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5),
                                              rep(3:4, 5, each=5)))/25,
                       winteringPos,
                       rowsum(winteringPos, rep(1:10, each=10))/10,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5),
                                              rep(3:4, 5, each=5)))/25,
                       winteringPos,
                       rowsum(winteringPos, rep(1:10, each=10))/10,
                       rowsum(winteringPos, c(rep(1:2, 5, each=5),
                                              rep(3:4, 5, each=5)))/25)

capLocs14 <- lapply(breedingPos14, SpatialPoints, proj4string = CRS(WGS84))
targLocs14 <- lapply(winteringPos14, SpatialPoints, proj4string = CRS(WGS84))


# Calculate distances between defined breeding populations
breedDist14 <- lapply(breedingPos14, distFromPos, surface = 'ellipsoid')

# Calculate distances between defined non-breeding populations
nonbreedDist14 <- lapply(winteringPos14, distFromPos, surface = 'ellipsoid')

# Numbers of defined populations
nBreeding14 <- c(100, 10, 100, 10, 4, 100, 4, 10, 4, 10, 10, 10, 4, 4, 4)

nWintering14 <- c(100, 100, 10, 10, 100, 4, 10, 4, 4, 100, 10, 4, 100, 10, 4)

# Relative abundance by scenario and breeding population
breedingN14base <- rep(25, nBreeding14[1])
breedingN14 <- lapply(breedingSiteTrans14, rowsum, x=breedingN14base)
breedingRelN14 <- lapply(breedingN14, "/", sum(breedingN14base))

MC <- 0.25
o <- optimize(mlogitMC, MC.in = MC, origin.dist = breedDist,
              target.dist = nonbreedDist, origin.abund = breedingN14[[1]],
              interval = c(0, 10), tol = .Machine$double.eps^0.5)

slope <- o$minimum

psi <- mlogitMat(slope, breedDist)


nSample14 <- 100 # Total number sampled per simulation

# How many sampled from each natural population (sampling scenarios separate
# from definition scenarios)
sampleBreeding14 <- list(round(breedingRelN14[[1]]*nSample14),
                         c(rep(0, 22), round(breedingRelN14[[5]][1]*nSample14),
                           rep(0, 4), round(breedingRelN14[[5]][2]*nSample14),
                           rep(0, 44), round(breedingRelN14[[5]][3]*nSample14),
                           rep(0, 4), round(breedingRelN14[[5]][4]*nSample14),
                           rep(0, 22)),
                         rep(c(rep(0, 4),
                               rep(round(breedingRelN14[[2]][1]*nSample14/2), 2),
                               rep(0, 4)), 10))

lapply(sampleBreeding14, matrix, nrow=10, ncol=10)

# for the baseline use the simulation from above, sims
animalLoc14base <- sims$animalLoc

#transferring the simulated bird locations from the true populations to the researcher defined populations
changeLocations <- function(animalLoc, breedingSiteTrans, winteringSiteTrans) {
  animalLoc[,1,,] <- breedingSiteTrans[animalLoc[,1,,]]
  animalLoc[,2,,] <- winteringSiteTrans[animalLoc[,2,,]]
  return(animalLoc)
}

# Number of scenarios and number of simulations to run
nScenarios14 <- length(breedingSiteTrans14)
nSims14 <- 100
nSimsLarge14 <- 2500#0
nYears <- 1
nMonths <- 1

# Connections between scenarios and sampling regimes
scenarioToSampleMap14 <- c(rep(1, 9), rep(3, 3), rep(2, 3))

# Set up data structures for storing results
animalLoc14 <- vector("list", nScenarios14) #making an empty list to fill
sim14 <- vector("list", nSimsLarge14) #making an empty list to fill

compare14 <- data.frame(Scenario = c("True", nScenarios14),
                        MC = c(MC, rep(NA, nScenarios14)),
                        Mantel = c(MC, rep(NA, nScenarios14)))

compare14.array <- array(NA, c(nSimsLarge14, nScenarios14, 4), dimnames =
                           list(1:nSimsLarge14,
                                scenarios14,
                                c("MC", "MCA", "MCss", "Mantel")))

results14 <- vector("list", nScenarios14)

# Run simulations
set.seed(7)
system.time(for (sim in 101:nSimsLarge14) {
  cat("Simulation", sim, "of", nSimsLarge14, '\n')
  sim14[[sim]] <- lapply(sampleBreeding14, simMove, breedingDist = breedDist14[[1]],
                  winteringDist=nonbreedDist14[[1]], psi=psi, nYears=nYears,
                  nMonths=nMonths)
  for (i in 1:nScenarios14) {
    animalLoc14[[i]] <- changeLocations(sim14[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                        breedingSiteTrans14[[i]], winteringSiteTrans14[[i]])
    results14[[i]] <- calcPsiMC(breedDist14[[i]], nonbreedDist14[[i]],
                                animalLoc14[[i]],
                                originRelAbund = breedingRelN14[[i]],
                                verbose = F)
    compare14.array[sim, i, 'MC'] <- results14[[i]]$MC
    compare14.array[sim, i, 'MCA'] <- calcMC(breedDist14[[i]], nonbreedDist14[[i]],
                                             results14[[i]]$psi,
                                             originAbund = breedingN14[[i]])
    compare14.array[sim, i, 'MCss'] <- calcMC(breedDist14[[i]], nonbreedDist14[[i]],
                                             results14[[i]]$psi,
                                             originAbund = table(animalLoc14[[i]][,1,1,1]))
    compare14.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist14[[1]],
                                                         nonbreedDist14[[1]],
                                                         sim14[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                                         resamp=0)$correlation
  }
})

compare14.array[1:10, , ]
means14 <- apply(compare14.array, 2:3, mean, na.rm = T)
vars14 <- apply(compare14.array, 2:3, var, na.rm = T)
rmse14 <- apply(compare14.array[1:nSimsLarge14, , ], 2:3, function(x) sqrt(mean((x - MC)^2)))


# Run estimations
sim14.sub <- sim14[sample.int(nSimsLarge14, nSims14, T)]
for (sim in 1:nSims14) {
  cat("Simulation", sim, "of", nSims14, '\n')
  for (i in 1:nScenarios14) {
    cat("\tScenario", i, "\n")
    animalLoc14[[i]] <- changeLocations(sim14.sub[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                        breedingSiteTrans14[[i]], winteringSiteTrans14[[i]])
    results14[[i]] <- estMC(breedDist14[[i]], nonbreedDist14[[i]],
                            originRelAbund = breedingRelN14[[i]],
                            originPoints = capLocs14[[i]][animalLoc14[[i]][, 1, 1, 1], ],
                            targetPoints = targLocs14[[i]][animalLoc14[[i]][, 2, 1, 1], ],
                            originAssignment = animalLoc14[[i]][, 1, 1, 1],
                            targetAssignment = animalLoc14[[i]][, 2, 1, 1],
                            nSamples = 1000, verbose = 0, calcCorr = T)
    compare14.array[sim, i, 'MC'] <- results14[[i]]$MC
    compare14.array[sim, i, 'Mantel'] <- calcStrengthInd(breedDist14[[1]],
                                                         nonbreedDist14[[1]],
                                                         sim14[[scenarioToSampleMap14[i]]]$animalLoc,
                                                         resamp=0)$correlation
  }
}

# Compute means for each scenario
compare14$MC[1:nScenarios14 + 1] <- apply(compare14.array[,,'MC'], 2, mean)
compare14$Mantel[1:nScenarios14 + 1] <- apply(compare14.array[,,'Mantel'], 2, mean)

compare14 <- transform(compare14, MC.diff=MC - MC[1], Mantel.diff=Mantel - Mantel[1])
compare14
