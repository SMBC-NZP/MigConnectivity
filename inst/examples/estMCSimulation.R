calcPsiMC <- function(originDist, targetDist, originRelAbund, locations, verbose=F) {
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
          psiMat[originMat[bi], targetMat[wi]] <- psiMat[originMat[bi], targetMat[wi]] + 1
  }
  psiMat <- apply(psiMat, 2, "/", rowSums(psiMat))
  MC <- calcMC(originDist, targetDist, psi = psiMat,
               originRelAbund = originRelAbund, sampleSize = nInd)
  return(list(psi=psiMat, MC=MC))
}

set.seed(75)
\dontrun{
nBreeding <- 100 # Number of populations
nWintering <- 100 # Number of populations
breedingPos <- matrix(c(rep(seq(-99,-81,2), each=sqrt(nBreeding)),
                        rep(seq(49,31,-2), sqrt(nBreeding))), nBreeding, 2)
winteringPos <- matrix(c(rep(seq(-79,-61,2), each=sqrt(nWintering)),
                         rep(seq(9,-9,-2), sqrt(nWintering))), nWintering, 2)

capLocs<-sp::SpatialPoints(breedingPos, sp::CRS(WGS84))
targLocs<-sp::SpatialPoints(winteringPos, sp::CRS(WGS84))
#capLocsM<-spTransform(capLocs, CRS(Lambert))
#targLocsM<-spTransform(targLocs, CRS(Lambert))

breedDist <- MigConnectivity::distFromPos(breedingPos, 'ellipsoid') # calculate distance between populations
nonbreedDist <- MigConnectivity::distFromPos(winteringPos, 'ellipsoid') # calculate distance between populations


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
              target.dist = nonbreedDist, origin.abund = breedingRelN14[[1]],
              sample.size = sum(breedingN14[[1]]),
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
nSimsLarge14 <- 2500
nYears <- 1
nMonths <- 1

# Connections between scenarios and sampling regimes
scenarioToSampleMap14 <- c(rep(1, 9), rep(3, 3), rep(2, 3))

# Set up data structures for storing results
animalLoc14 <- vector("list", nScenarios14) #making an empty list to fill
sim14 <- vector("list", nSimsLarge14) #making an empty list to fill

compare14.array <- array(NA, c(nSimsLarge14, 2), dimnames =
                           list(1:nSimsLarge14,
                                c("MC", "Mantel")))

results14 <- vector("list", nScenarios14)

# Run simulations
set.seed(7)
system.time(for (sim in 1:nSimsLarge14) {
  cat("Simulation", sim, "of", nSimsLarge14, '\n')
  sim14[[sim]] <- lapply(sampleBreeding14, simMove, breedingDist = breedDist14[[1]],
                  winteringDist=nonbreedDist14[[1]], psi=psi, nYears=nYears,
                  nMonths=nMonths)
  for (i in 13) {
    animalLoc14[[i]] <- changeLocations(sim14[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                        breedingSiteTrans14[[i]], winteringSiteTrans14[[i]])
    results14[[i]] <- calcPsiMC(breedDist14[[i]], nonbreedDist14[[i]],
                                animalLoc14[[i]],
                                originRelAbund = breedingRelN14[[i]],
                                verbose = F)
    compare14.array[sim, 'MC'] <- results14[[i]]$MC
    compare14.array[sim, 'Mantel'] <- calcStrengthInd(breedDist14[[1]],
                                                      nonbreedDist14[[1]],
                                                      sim14[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                                      resamp=0)$correlation
  }
})

compare14.array[1:10, ]
means14 <- apply(compare14.array, 2, mean)
vars14 <- apply(compare14.array, 2, var)
rmse14 <- apply(compare14.array, 2, function(x) sqrt(mean((x - MC)^2)))

# Set up data structures for storing estimation results
est14.array <- array(NA, c(nSims14, 2), dimnames =
                           list(1:nSims14,
                                c("MC", "Mantel")))
var14.array <- array(NA, c(nSims14, 2), dimnames =
                           list(1:nSims14,
                                c("MC", "Mantel")))
ci14.array <- array(NA, c(nSims14, 2, 2), dimnames =
                           list(1:nSims14,
                                c("MC", "Mantel"),
                                c('lower', 'upper')))
animalLoc14 <- vector("list", nSims14) #making an empty list to fill
results14 <- vector("list", nSims14)

# Run estimations
set.seed(567)
sim14.sub <- sim14[sample.int(nSimsLarge14, nSims14, T)]
for (sim in 1:nSims14) {
  cat("Estimation", sim, "of", nSims14, '\n')
  for (i in 13) {#:nScenarios14) {
    animalLoc14[[sim]] <- changeLocations(sim14.sub[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                        breedingSiteTrans14[[i]], winteringSiteTrans14[[i]])
    results14[[sim]] <- estMC(breedDist14[[i]], nonbreedDist14[[i]],
                            originRelAbund = breedingRelN14[[i]],
                            originPoints = capLocs14[[i]][animalLoc14[[sim]][, 1, 1, 1], ],
                            targetPoints = targLocs14[[i]][animalLoc14[[sim]][, 2, 1, 1], ],
                            originAssignment = animalLoc14[[sim]][, 1, 1, 1],
                            targetAssignment = animalLoc14[[sim]][, 2, 1, 1],
                            nSamples = 1000, verbose = 0, calcCorr = T,
                            geoBias = c(0, 0), geoVCov = matrix(0, 2, 2))
    est14.array[sim, 'MC'] <- results14[[sim]]$meanMC
    est14.array[sim, 'Mantel'] <- results14[[sim]]$meanCorr
    var14.array[sim, 'MC'] <- results14[[sim]]$seMC ^ 2
    var14.array[sim, 'Mantel'] <- results14[[sim]]$seCorr ^ 2
    ci14.array[sim, 'MC', ] <- results14[[sim]]$bcCI
    ci14.array[sim, 'Mantel', ] <- results14[[sim]]$bcCICorr
    save(results14, est14.array, var14.array, ci14.array, file = 'results14c.gzip')
  }
}

summary(var14.array)
vars14
summary(est14.array)
means14
summary(ci14.array[, "MC", 'lower'] <= MC & ci14.array[, "MC", 'upper'] >= MC)
summary(ci14.array[, "Mantel", 'lower'] <= MC & ci14.array[, "Mantel", 'upper'] >= MC)

est.df <- data.frame(Parameter = rep(c("MC", 'rM'), 2, each = nSims14),
                     Quantity = rep(c('Mean', 'Variance'), each = 2 * nSims14),
                     sim = rep(1:nSims14, 4), Estimate = c(est14.array, var14.array))
trues.df <- data.frame(Parameter = rep(c("MC", 'rM'), 6),
                       Quantity = rep(c('Mean', 'Variance'), c(4, 2)),
                       Source = rep(c('"True" Value', 'Larger Simulation'), c(2, 4)),
                       Value = c(MC, MC, means14, vars14))
library(ggplot2)
g.est <- ggplot(est.df, aes(Estimate)) + geom_histogram(bins = 15) +
  facet_grid(Parameter ~ Quantity, scales = 'free_x') +
  geom_vline(aes(xintercept = Value, color = Source), data = trues.df) +
  theme_bw() + scale_x_continuous(breaks = c(0.002, 0.003, 0.1, 0.2, 0.3))
png('sim.est.MC.rM3.png', width = 6, height = 4, units = 'in', res=600)
g.est
dev.off()

qualities14 <- data.frame(Parameter = rep(c("MC", 'rM'), 8),
                       Quantity = rep(rep(c('Mean', 'Variance'), c(4, 2)), 3, 16),
                       Source = rep(rep(c('"True" Value', 'Larger Simulation'), c(2, 4)), 3, 16),
                       Measure = rep(c('Bias', 'RMSE', "Coverage"), c(6, 6, 4)),
                       Value = c(colMeans(est14.array - MC),
                                 colMeans(est14.array) - means14,
                                 colMeans(var14.array) - vars14,
                                 sqrt(colMeans((est14.array - MC)^2)),
                                 sqrt(mean((est14.array[,1] - means14[1])^2)),
                                 sqrt(mean((est14.array[,2] - means14[2])^2)),
                                 sqrt(mean((var14.array[,1] - vars14[1])^2)),
                                 sqrt(mean((var14.array[,2] - vars14[2])^2)),
                                 mean(ci14.array[, "MC", 'lower'] <= MC & ci14.array[, "MC", 'upper'] >= MC),
                                 mean(ci14.array[, "Mantel", 'lower'] <= MC & ci14.array[, "Mantel", 'upper'] >= MC),
                                 mean(ci14.array[, "MC", 'lower'] <= means14[1] & ci14.array[, "MC", 'upper'] >= means14[1]),
                                 mean(ci14.array[, "Mantel", 'lower'] <= means14[2] & ci14.array[, "Mantel", 'upper'] >= means14[2])))
format(qualities14, digits = 2, scientific = F)

# point.rM14 <- sapply(results14, function(x) x$pointCorr, simplify = 'array')
# mean(simplify2array(point.rM14))

save.image('estMCsims.c.RData')

# Try with biased scenario (5)
est14.array <- array(NA, c(nSims14, nScenarios14, 3), dimnames =
                           list(1:nSims14, scenarios14,
                                c("MC", "MCss", "Mantel")))
var14.array <- array(NA, c(nSims14, nScenarios14, 3), dimnames =
                           list(1:nSims14, scenarios14,
                                c("MC", "MCss", "Mantel")))
ci14.array <- array(NA, c(nSims14, nScenarios14, 3, 2), dimnames =
                           list(1:nSims14, scenarios14,
                                c("MC", "MCss", "Mantel"),
                                c('lower', 'upper')))
set.seed(567)
sim14.sub <- sim14[sample.int(nSimsLarge14, nSims14, T)]
for (sim in 1:nSims14) {
  cat("Estimation", sim, "of", nSims14, '\n')
  for (i in 5) {#:nScenarios14) {
    animalLoc14[[sim]] <- changeLocations(sim14.sub[[sim]][[scenarioToSampleMap14[i]]]$animalLoc,
                                        breedingSiteTrans14[[i]], winteringSiteTrans14[[i]])
    results14[[sim]] <- estMC(breedDist14[[i]], nonbreedDist14[[i]],
                            originRelAbund = breedingRelN14[[i]],
#                            originPoints = capLocs14[[i]][animalLoc14[[sim]][, 1, 1, 1], ],
#                           targetPoints = targLocs14[[i]][animalLoc14[[sim]][, 2, 1, 1], ],
                            originAssignment = animalLoc14[[sim]][, 1, 1, 1],
                            targetAssignment = animalLoc14[[sim]][, 2, 1, 1],
                            nSamples = 1000, verbose = 0, calcCorr = F)
    est14.array[sim, i, 'MCss'] <- results14[[sim]]$meanMC
    var14.array[sim, i, 'MCss'] <- results14[[sim]]$seMC ^ 2
    ci14.array[sim, i, 'MCss', ] <- results14[[sim]]$bcCI
    save(results14, est14.array, var14.array, ci14.array, file = 'results14d.gzip')
  }
}

rmse <- function(ests, truth) {
  return(sqrt(mean((ests - truth) ^ 2)))
}

bias <- function(ests, truth) {
  return(mean(ests - truth))
}

coverage <- function(ci, truth) {
  return(mean(ci[1,] <= truth & ci[2,] >= truth))
}

rmse(sapply(results14, function(x) x$pointMC), 0.25)
rmse(sapply(results14, function(x) x$meanMC), 0.25)
bias(sapply(results14, function(x) x$pointMC), 0.25)
bias(sapply(results14, function(x) x$meanMC), 0.25)
mean(ci14.array[, 5, "MCss", 'lower'] <= MC & ci14.array[, 5, "MCss", 'upper'] >= MC)
mean(ci14.array[, 5, "MCss", 'lower'] <= MC & ci14.array[, 5, "MCss", 'upper'] >= MC)
coverage(sapply(results14, function(x) x$simpleCI), 0.25)

###############################################################################
# Try with location uncertainty (GL data)
###############################################################################

# Assign geolocator bias / variance co-variance matrix
geoBias <- OVENdata$geo.bias
geoVCov <- OVENdata$geo.vcov

# Helper function to convert a string of XY coordinates of centroids to polygons #
toPoly <- function(siteCentroids,
                   projection.in = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                   projection.out = NA,
                   resolution = NA){
            # This automatically sets the resolution so that all polygons touch and cover the entire surface #
            # Alternatively the user can supply the resolution of the raster cells in the units of the
            # the input projection (projection.in) #
                   if(is.na(resolution)){
                     long <- unique(siteCentroids[,1])
                     lat <- unique(siteCentroids[,2])
                     long.res <- long[2]-long[1]
                     lat.res <- lat[2]-lat[1]
                     resolution <- c(long.res,long.res)
                   }
                             rast <- raster::rasterFromXYZ(cbind(cbind(siteCentroids),rep(1,nrow(siteCentroids))),
                             res = resolution,
                             crs = projection)
                      polys <- raster::rasterToPolygons(rast)
                      raster::crs(polys)<-projection.in

                      if(!is.na(projection.out)){
                      polys <- sp::spTransform(polys,sp::CRS(projection.out))
                      }

                      return(polys)}

# Define projections
WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Lambert <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# convert wintering locations to ploygons using the helper function
polys <- lapply(winteringPos14,toPoly,projection.in = WGS84, projection.out = NA)


simLocationError <- function(targetPoints, targetSites, geoBias, geoVCov) {
  nAnimals <- length(targetPoints)
  nSim <- 100
  geoBias2 <- matrix(rep(geoBias, nSim), nrow=nSim, byrow=T)
  target.sample <- rep(NA, nAnimals)
  target.point.sample <- matrix(NA, nAnimals, 2)
  for(i in 1:nAnimals){
    draws <- 0
    while (is.na(target.sample[i])) {
      draws <- draws + 1
      # Sample random point for each bird from parametric distribution of NB error
      point.sample <- sp::SpatialPoints(MASS::mvrnorm(n=nSim, mu=cbind(
        targetPoints@coords[animal.sample[i],1],
        targetPoints@coords[animal.sample[i],2]), Sigma=geoVCov)+
          geoBias2, sp::CRS(Lambert))
      # filtered to stay in NB areas (land)
      target.sample0 <- sp::over(point.sample, targetSites)
      target.sample[i]<-target.sample0[!is.na(target.sample0)][1]
    }
    target.point.sample[i, ]<-point.sample[!is.na(target.sample0)][1]@coords
  }
  return(targetSample = target.sample, targetPointSample = target.point.sample)
}
}

