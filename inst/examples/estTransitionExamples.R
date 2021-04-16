# Test new code for sampling different types of data
###############################################################################
# Example 1 (all released origin; some telemetry, some GL, some prob
# some both GL and prob)
###############################################################################
library(VGAM)
library(MigConnectivity)
nAnimals <- 40
isGL <- c(OVENdata$isGL, FALSE)
#isGL[5] <- FALSE
isTelemetry <- c(!OVENdata$isGL, FALSE)
isRaster <- rep(FALSE, nAnimals)
isProb <- rep(FALSE, nAnimals)
targetPoints <- rbind(OVENdata$targetPoints, OVENdata$targetPoints[1,])
targetSites <- OVENdata$targetSites
originSites <- OVENdata$originSites
resampleProjection <- sf::st_crs(OVENdata$targetPoints)
targetPoints <- sf::st_transform(targetPoints, crs = resampleProjection)
targetSites <- sf::st_transform(targetSites, crs = resampleProjection)
targetNames <- OVENdata$targetNames
originNames <- OVENdata$originNames
targetAssignment <- array(0, dim = c(nAnimals, 3), dimnames = list(NULL, targetNames))
assignment0 <- unclass(sf::st_intersects(x = targetPoints, y = targetSites,
                                         sparse = TRUE))
assignment0[sapply(assignment0, function(x) length(x)==0)] <- 0
assignment0 <- array(unlist(assignment0), nAnimals)
for (ani in 1:nAnimals) {
  if (assignment0[ani]>0)
    targetAssignment[ani, assignment0[ani]] <- 1
  else{
    targetAssignment[ani, ] <- rdiric(1, c(15, 1, 1))
    isProb[ani] <- TRUE
  }
}
targetAssignment
isProb
system.time(test <- MigConnectivity:::locSample(isGL, isRaster, isProb, isTelemetry,
                  geoBias = OVENdata$geo.bias, #[, 2:1, drop = FALSE]
                  geoVCov = OVENdata$geo.vcov,#*1.5,#[2:1,2:1]
                  points = targetPoints, assignment = targetAssignment,
                  sites = targetSites, resampleProjection = resampleProjection,
                  nSim = 20000, maxTries = 300))
test
originPoints <- rbind(OVENdata$originPoints,
                      OVENdata$originPoints[39,])
originPoints <- sf::st_transform(originPoints, crs = resampleProjection)
originSites <- sf::st_transform(OVENdata$originSites, crs = resampleProjection)
system.time(test2 <-
              MigConnectivity:::estTransition(isGL = isGL, isRaster = isRaster,
                                              isProb = isProb,
                                              isTelemetry = isTelemetry,
                                              geoBias = OVENdata$geo.bias, #[, 2:1, drop = FALSE]
                                              geoVCov = OVENdata$geo.vcov,#*1.5,#[2:1,2:1]
                                              targetPoints = targetPoints,
                                              targetAssignment = targetAssignment,
                                              targetSites = targetSites,
                                              resampleProjection = resampleProjection,
                                              nSim = 5000, maxTries = 300,
                                              originSites = originSites,
                                              originPoints = originPoints,
                                              captured = "origin",
                                              originNames = OVENdata$originNames,
                                              targetNames = OVENdata$targetNames,
                                              verbose = 3,
                                              nSamples = 1000))
test2
plot(test2, legend = "top",
     main = paste("OVENlike w/", sum(isGL & !isProb), "GL,",
                  sum(!isGL & isProb), "probs,",
                  sum(isGL & isProb), "both, and", sum(isTelemetry), "GPS"),
     col = RColorBrewer::brewer.pal(nNonBreeding, "Dark2"))
plot(targetPoints)
plot(targetSites, add = T)
plot(targetSites)
plot(targetPoints[isProb & isGL, ], add = T, col = 1:sum(isProb & isGL))

################################################################
# Example 2 (add prob animals released on other end)
################################################################
nAnimals <- 45
captured <- rep(c("origin", "target"), c(40, 5))
isGL <- c(OVENdata$isGL, rep(FALSE, 6))
#isGL[5] <- FALSE
isTelemetry <- c(!OVENdata$isGL, rep(FALSE, 6))
isRaster <- rep(FALSE, nAnimals)
isProb <- rep(FALSE, nAnimals)
targetPoints <- rbind(OVENdata$targetPoints, OVENdata$targetPoints[c(1:4,6:7),])
targetSites <- OVENdata$targetSites
originSites <- OVENdata$originSites
resampleProjection <- sf::st_crs(OVENdata$targetPoints)
targetPoints <- sf::st_transform(targetPoints, crs = resampleProjection)
targetSites <- sf::st_transform(targetSites, crs = resampleProjection)
targetNames <- OVENdata$targetNames
originNames <- OVENdata$originNames
targetAssignment <- array(0, dim = c(nAnimals, 3), dimnames = list(NULL, targetNames))
assignment0 <- unclass(sf::st_intersects(x = targetPoints, y = targetSites,
                                         sparse = TRUE))
assignment0[sapply(assignment0, function(x) length(x)==0)] <- 0
assignment0 <- array(unlist(assignment0), nAnimals)
assignment0[nAnimals] <- 3
for (ani in 1:nAnimals) {
  if (assignment0[ani]>0)
    targetAssignment[ani, assignment0[ani]] <- 1
  else{
    targetAssignment[ani, ] <- rdiric(1, c(15, 1, 1))
    isProb[ani] <- TRUE
  }
}
targetAssignment
isProb
originPoints <- rbind(OVENdata$originPoints,
                      OVENdata$originPoints[34:39,])
originPoints <- sf::st_transform(originPoints, crs = resampleProjection)
originSites <- sf::st_transform(OVENdata$originSites, crs = resampleProjection)
assignment1 <- unclass(sf::st_intersects(x = originPoints, y = originSites,
                                         sparse = TRUE))
assignment1[sapply(assignment1, function(x) length(x)==0)] <- 0
assignment1 <- array(unlist(assignment1), nAnimals)
originAssignment <- array(0, dim = c(nAnimals, nOriginSites),
                          dimnames = list(NULL, originNames))
for (ani in 1:40) {
  originAssignment[ani, assignment1[ani]] <- 1
}
for (ani in 41:nAnimals) {
  originAssignment[ani, ] <- rdiric(1, c(1, 1))
  isProb[ani] <- TRUE
}
originAssignment
isProb
system.time(test3 <-
              MigConnectivity:::estTransition(isGL = isGL, isRaster = isRaster,
                                              isProb = isProb,
                                              isTelemetry = isTelemetry,
                                              geoBias = OVENdata$geo.bias, #[, 2:1, drop = FALSE]
                                              geoVCov = OVENdata$geo.vcov,#*1.5,#[2:1,2:1]
                                              targetPoints = targetPoints,
                                              targetAssignment = targetAssignment,
                                              targetSites = targetSites,
                                              resampleProjection = resampleProjection,
                                              nSim = 5000, maxTries = 300,
                                              originSites = originSites,
                                              #originPoints = originPoints,
                                              originAssignment = originAssignment,
                                              captured = captured,
                                              originNames = OVENdata$originNames,
                                              targetNames = OVENdata$targetNames,
                                              verbose = 3,
                                              nSamples = 100))
test3
plot(test3, legend = "top",
     main = paste(sum(isGL & !isProb), "GL,",
                  sum(!isGL & isProb & captured == "origin"), "prob,",
                  sum(isGL & isProb), "both,",
                  sum(isTelemetry), "GPS (all\ncaptured origin), and",
                  sum(isProb & captured == "target"), "probs (captured target)"),
     col = RColorBrewer::brewer.pal(nNonBreeding, "Dark2"))
MC3 <- MigConnectivity:::estStrength(OVENdata$originDist, OVENdata$targetDist,
                                     OVENdata$originRelAbund, test3,
                                     sampleSize = nAnimals)
MC3

###############################################################################
# Example 3
###############################################################################
