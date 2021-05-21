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
              estTransition(isGL = isGL, isRaster = isRaster,
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

nOriginSites <- nrow(originSites)

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
              estTransition(isGL = isGL, isRaster = isRaster,
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

nNonBreeding <- nrow(OVENdata$targetSites)

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
# Example 3 (all raster, from our OVEN example)
###############################################################################
getCSV <- function(filename) {
  tmp <- tempdir()
  url1 <- paste0('https://github.com/SMBC-NZP/MigConnectivity/blob/isodev/data-raw/',
                 filename, '?raw=true')
  temp <- paste(tmp, filename, sep = '/')
  utils::download.file(url1, temp, mode = 'wb')
  csv <- read.csv(temp)
  unlink(temp)
  return(csv)

}
getRDS <- function(speciesDist) {
  tmp <- tempdir()
  extension <- '.rds'
  filename <- paste0(speciesDist, extension)
  url1 <- paste0('https://github.com/SMBC-NZP/MigConnectivity/blob/isodev/data-raw/Spatial_Layers/',
                 filename, '?raw=true')
  temp <- paste(tmp, filename, sep = '/')
  utils::download.file(url1, temp, mode = 'wb')
  shp <- readRDS(temp)
  unlink(temp)
  return(shp)
}
OVENdist <- getRDS("OVENdist")

raster::crs(OVENdist) <- MigConnectivity::projections$WGS84

OVENvals <- getCSV("deltaDvalues.csv")

OVENvals <- OVENvals[grep(x=OVENvals$Sample,"NH", invert = TRUE),]

originSites <- getRDS("originSites")
originDist <- distFromPos(rgeos::gCentroid(originSites,byid = TRUE)@coords)
originSites <- sf::st_as_sf(originSites)

EVER <- length(grep(x=OVENvals$Sample,"EVER"))
JAM <- length(grep(x=OVENvals$Sample,"JAM"))

originRelAbund <- matrix(c(EVER,JAM),nrow = 1,byrow = TRUE)
originRelAbund <- prop.table(originRelAbund,1)

op <- sf::st_centroid(originSites)

originPoints <- array(NA,c(EVER+JAM,2), list(NULL, c("x","y")))
originPoints[grep(x = OVENvals$Sample,"JAM"),1] <- sf::st_coordinates(op)[1,1]
originPoints[grep(x = OVENvals$Sample,"JAM"),2] <- sf::st_coordinates(op)[1,2]
originPoints[grep(x = OVENvals$Sample,"EVER"),1] <- sf::st_coordinates(op)[2,1]
originPoints[grep(x = OVENvals$Sample,"EVER"),2] <- sf::st_coordinates(op)[2,2]

originPoints <- sf::st_as_sf(data.frame(originPoints),
                             coords = c("x", "y"),
                             crs = 4326)

iso <- isoAssign(isovalues = OVENvals[,2],
                 isoSTD = 12,       # this value is for demonstration only
                 intercept = -10,   # this value is for demonstration only
                 slope = 0.8,       # this value is for demonstration only
                 odds = NULL,
                 restrict2Likely = TRUE,
                 nSamples = 1000,
                 sppShapefile = OVENdist,
                 assignExtent = c(-179,-60,15,89),
                 element = "Hydrogen",
                 surface = FALSE,
                 period = "Annual",
                 seed = 12345,
                 verbose=1)

nAnimals <- dim(iso$probassign)[3]
isGL <-rep(F, nAnimals); isRaster <- rep(T, nAnimals)
isProb <- rep(F, nAnimals); isTelemetry <- rep(F, nAnimals)
targetSites <- rgeos::gUnaryUnion(iso$targetSites, id = iso$targetSites$targetSite)
targetSites <- sf::st_as_sf(targetSites)


system.time(test4 <-
              estTransition(isGL = isGL,
                                              isRaster = isRaster,
                                              isProb = isProb,
                                              isTelemetry = isTelemetry,
                                              #geoBias = OVENdata$geo.bias, #[, 2:1, drop = FALSE]
                                              #geoVCov = OVENdata$geo.vcov,#*1.5,#[2:1,2:1]
                                              #targetPoints = targetPoints,
                                              #targetAssignment = targetAssignment,
                                              targetSites = targetSites,
                                              resampleProjection = resampleProjection,
                                              targetRaster = iso,
                                              #nSim = 5000, maxTries = 300,
                                              originSites = originSites,
                                              originPoints = originPoints,
                                              #originAssignment = originAssignment,
                                              captured = rep("origin", nAnimals),
                                              #originNames = OVENdata$originNames,
                                              #targetNames = OVENdata$targetNames,
                                              verbose = 3,
                                              nSamples = 1000))
test4

ovenMC <- estMC(originRelAbund = originRelAbund,
                targetIntrinsic = iso,
                originPoints = originPoints,
                originSites = originSites,
                targetSites = targetSites,
                originDist = originDist,
                nSamples = 200,
                verbose = 1,
                calcCorr = TRUE,
                alpha = 0.05,
                approxSigTest = FALSE,
                isIntrinsic = TRUE)
ovenMC

test4$psi$mean - ovenMC$psi$mean
test4$psi$se - ovenMC$psi$se

iso2 <- iso
iso2$SingleCell <- NULL
something <- MigConnectivity:::locSample(isGL = isGL,
                                             isRaster = isRaster,
                                             isProb = isProb,
                                             isTelemetry = isTelemetry,
                                             #geoBias = OVENdata$geo.bias, #[, 2:1, drop = FALSE]
                                             #geoVCov = OVENdata$geo.vcov,#*1.5,#[2:1,2:1]
                                             #targetPoints = targetPoints,
                                             #targetAssignment = targetAssignment,
                                             sites = targetSites,
                                             resampleProjection = resampleProjection,
                                             matvals = raster::rasterToPoints(iso2$probassign),
                                            nSim = 10, maxTries = 300,
                                         assignment = targetAssignment)

system.time(test5 <-
              estTransition(isGL = isGL,
                                              isRaster = isRaster,
                                              isProb = isProb,
                                              isTelemetry = isTelemetry,
                                              #geoBias = OVENdata$geo.bias, #[, 2:1, drop = FALSE]
                                              #geoVCov = OVENdata$geo.vcov,#*1.5,#[2:1,2:1]
                                              #targetPoints = targetPoints,
                                              #targetAssignment = targetAssignment,
                                              targetSites = targetSites,
                                              resampleProjection = resampleProjection,
                                              targetRaster = iso2,
                                              #nSim = 5000, maxTries = 300,
                                              originSites = originSites,
                                              originPoints = originPoints,
                                              #originAssignment = originAssignment,
                                              captured = rep("origin", nAnimals),
                                              #originNames = OVENdata$originNames,
                                              #targetNames = OVENdata$targetNames,
                                              verbose = 3,
                                              nSamples = 10))
test5$psi$mean
test4$psi$mean - test5$psi$mean
test4$psi$se - test5$psi$se

ovenMC5 <- suppressMessages(estMC(originRelAbund = originRelAbund,
                targetIntrinsic = iso2,
                originPoints = originPoints,
                originSites = originSites,
                targetSites = targetSites,
                originDist = originDist,
                nSamples = 200,
                verbose = 1,
                calcCorr = F,
                alpha = 0.05,
                approxSigTest = FALSE,
                isIntrinsic = TRUE))
ovenMC5



# Don't run, only works if you've run a bunch of stuff from within functions
ts <- sf::st_transform(targetSites, resampleProjection)
ps2 <- lapply(point.sample2, sf::st_transform, crs = resampleProjection)
tar <- sf::st_transform(targetAssignRast, resampleProjection)

plot(ts)
for(i in 2) {
  plot(ps2[[i]], col = RColorBrewer::brewer.pal(7, "Dark2")[i], add = T)
}
plot(tar[1:7, ], col = RColorBrewer::brewer.pal(7, "Dark2"), add = T, pch = 19)
sc1 <- sf::st_as_sf(as.data.frame(iso$SingleCell[,,1]), coords = c("Longitude", "Latitude"),
                    crs = MigConnectivity::projections$WGS84)
sc1 <- sf::st_transform(sc1, resampleProjection)
plot(sc1, col = RColorBrewer::brewer.pal(7, "Dark2")[1], add = T, pch = 5)
sc2 <- sf::st_as_sf(as.data.frame(iso$SingleCell[,,2]), coords = c("Longitude", "Latitude"),
                    crs = MigConnectivity::projections$WGS84)
sc2 <- sf::st_transform(sc2, resampleProjection)
plot(sc2, col = RColorBrewer::brewer.pal(7, "Dark2")[2], add = T, pch = 5, alpha = 0.5)

###############################################################################
# Example 4 (raster plus prob tables) (haven't gotten far setting this one up)
###############################################################################
isProb <- rep(c(F, T), c(20, nAnimals - 20))
targetAssignment <- array(0, dim = c(nAnimals, 3), dimnames = list(NULL, targetNames))
xyTargetRast <- apply(targetRasterXYZ[,3:ncol(targetRasterXYZ)],
                      MARGIN = 2,
                      FUN = function(x){
                        #xy <-cbind(Hmisc::wtd.quantile(targetRasterXYZ[,1],probs = 0.5, weight = x, na.rm = TRUE),
                        #           Hmisc::wtd.quantile(targetRasterXYZ[,2],probs = 0.5, weight = x, na.rm = TRUE))
                        # xy <- cbind(weighted.mean(targetRasterXYZ[,1], w = x, na.rm = TRUE),
                        #            weighted.mean(targetRasterXYZ[,2], w = x, na.rm = TRUE))
                        # Select cell with the maximum posterior probability #
                        xy <- cbind(targetRasterXYZ[which.max(x)[1],1],
                                    targetRasterXYZ[which.max(x)[1],2])
                        return(xy)})
# returns a point estimate for each bird - turn it into a sf object
xyTargetRast <- t(xyTargetRast)
colnames(xyTargetRast) <- c("x","y")
# right now the assignment CRS is WGS84 - should be the same as the origin raster
targetAssignRast <- sf::st_as_sf(data.frame(xyTargetRast), coords = c("x","y"), crs = 4326)
# transform to match originSites
targetSites_wgs <- sf::st_transform(targetSites, 4326)
#targetAssignRast <- sf::st_transform(targetAssignRast, sf::st_crs(targetSites))
#targetAssignment <- suppressMessages(array(unlist(unclass(sf::st_intersects(x = targetAssignRast,
#                                                           y = targetSites_wgs,
#                                                           sparse = TRUE)))))
# Check which points are in target sites
ta_assign <- suppressMessages(sf::st_intersects(x = targetAssignRast,
                                                y = targetSites_wgs,
                                                sparse = TRUE))

# Give values without intersection an NA
ta_intersect <- lengths(ta_assign)
# quick check to ensure that all the points fall exactly in one targetSite #
if(any(ta_intersect)>1){stop("Overlapping targetSites not allowed \n")}

ta_bool_intersect <- lengths(ta_intersect)>0

targetAssignment <- unlist(as.numeric(ta_assign))


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
