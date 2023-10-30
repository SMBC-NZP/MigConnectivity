\dontrun{
###############################################################################
# Examples 1 (banding data: first example is based on common tern banding data;
#   the second is made up data to demonstrate data with two ages)
###############################################################################
COTE_banded <- c(10360, 1787, 2495, 336)
COTE_reencountered <- matrix(c(12, 0, 38, 15,
                               111, 7, 6, 2,
                               5, 0, 19, 4,
                               1123, 40, 41, 7),
                             4, 4,
                             dimnames = list(LETTERS[1:4], 1:4))
COTE_psi <- estTransition(originNames = LETTERS[1:4],
                          targetNames = 1:4,
                          banded = COTE_banded,
                          reencountered = COTE_reencountered,
                          verbose = 1,
                          nSamples = 60000, nBurnin = 20000,
                          method = "MCMC")
COTE_psi

COTE_banded2 <- matrix(rep(COTE_banded, 2), 4, 2)
COTE_reencountered2 <- array(c(12, 0, 38, 15, 6, 0, 17, 7,
                               111, 7, 6, 2, 55, 3, 3, 1,
                               5, 0, 19, 4, 2, 0, 10, 2,
                               1123, 40, 41, 7, 660, 20, 20, 3),
                             c(4, 2, 4),
                             dimnames = list(LETTERS[1:4], c("J", "A"), 1:4))
COTE_psi2 <- estTransition(originNames = LETTERS[1:4],
                          targetNames = 1:4,
                          banded = COTE_banded2,
                          reencountered = COTE_reencountered2,
                          verbose = 0,
                          nSamples = 60000, nBurnin = 20000,
                          method = "MCMC")
COTE_psi2

###############################################################################
# Example 2 (geolocator and telemetry ovenbirds captured on origin sites)
###############################################################################
data(OVENdata) # Ovenbird

nSamplesGLGPS <- 100 # Number of bootstrap iterations
#\dontrun{
#  nSamplesGLGPS <- 10000 # Number of bootstrap iterations
#}

# Estimate MC only, treat all data as geolocator
GL_psi <- estTransition(isGL=TRUE,
                        geoBias = OVENdata$geo.bias,
                        geoVCov = OVENdata$geo.vcov,
                        targetSites = OVENdata$targetSites,
                        originSites = OVENdata$originSites,
                        originPoints = OVENdata$originPoints,
                        targetPoints = OVENdata$targetPoints,
                        verbose = 2,
                        nSamples = nSamplesGLGPS,
                        resampleProjection = sf::st_crs(OVENdata$targetPoints))

# Estimate MC and rM, treat all data as is
Combined.psi <- estTransition(isGL=OVENdata$isGL,
                        isTelemetry = !OVENdata$isGL,
                geoBias = OVENdata$geo.bias, # Light-level GL location bias
                geoVCov = OVENdata$geo.vcov, # Location covariance matrix
                targetSites = OVENdata$targetSites, # Non-breeding target sites
                originSites = OVENdata$originSites, # Breeding origin sites
                originPoints = OVENdata$originPoints, # Capture Locations
                targetPoints = OVENdata$targetPoints, # Device target locations
                verbose = 2,   # output options
                nSamples = nSamplesGLGPS, # This is set low for example
                resampleProjection = sf::st_crs(OVENdata$targetPoints))

print(Combined.psi)

# For treating all data as GPS,
# Move the latitude of birds with locations that fall offshore
int <- sf::st_intersects(OVENdata$targetPoints, OVENdata$targetSites)
any(lengths(int)<1)
plot(OVENdata$targetPoints)
plot(OVENdata$targetSites,add=TRUE)
tp<-sf::st_coordinates(OVENdata$targetPoints)
text(tp[,1], tp[,2], label=c(1:39))

tp[5,2] <- 2450000
tp[10,2]<- 2240496
tp[1,2]<- 2240496
tp[11,2]<- 2026511
tp[15,2]<- 2031268
tp[16,2]<- 2031268

oven_targetPoints<-sf::st_as_sf(as.data.frame(tp),
                                coords = c("X","Y"),
                                crs = sf::st_crs(OVENdata$targetPoints))
inter <- sf::st_intersects(oven_targetPoints, OVENdata$targetSites)
any(lengths(inter)<1)
plot(oven_targetPoints,add=TRUE, col = "green")
plot(oven_targetPoints[lengths(inter)<1,],add=TRUE, col = "darkblue")

# Estimate MC only, treat all data as GPS
GPS_psi <- estTransition(isTelemetry = TRUE,
              targetSites = OVENdata$targetSites, # Non-breeding target sites
              originSites = OVENdata$originSites, # Breeding origin sites
              originPoints = OVENdata$originPoints, # Capture Locations
              targetPoints = oven_targetPoints, # Device target locations
              verbose = 2,   # output options
              nSamples = nSamplesGLGPS) # This is set low for example



###############################################################################
# Example 3 (all released origin; some telemetry, some GL, some probability
# tables, some both GL and probability tables; data modified from ovenbird
# example)
###############################################################################
#\dontrun{
#  install.packages(c('VGAM'))
#}
library(VGAM)
nAnimals <- 40
isGL <- c(OVENdata$isGL, FALSE)
isTelemetry <- c(!OVENdata$isGL, FALSE)
isRaster <- rep(FALSE, nAnimals)
isProb <- rep(FALSE, nAnimals)
targetPoints <- rbind(OVENdata$targetPoints, OVENdata$targetPoints[1,])
targetSites <- OVENdata$targetSites
originSites <- OVENdata$originSites
resampleProjection <- sf::st_crs(OVENdata$targetPoints)
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
nSamplesTry <- 100 # Number of bootstrap iterations
originPoints <- rbind(OVENdata$originPoints,
                      OVENdata$originPoints[39,])
system.time(psi3 <-
              estTransition(isGL = isGL, isRaster = isRaster,
                            isProb = isProb,
                            isTelemetry = isTelemetry,
                            geoBias = OVENdata$geo.bias,
                            geoVCov = OVENdata$geo.vcov,
                            targetPoints = targetPoints,
                            targetAssignment = targetAssignment,
                            targetSites = targetSites,
                            resampleProjection = resampleProjection,
                            nSim = 20000, maxTries = 300,
                            originSites = originSites,
                            originPoints = originPoints,
                            captured = "origin",
                            originNames = OVENdata$originNames,
                            targetNames = OVENdata$targetNames,
                            verbose = 3,
                            nSamples = nSamplesTry))
psi3

nNonBreeding <- nrow(OVENdata$targetSites)

plot(psi3, legend = "top",
     main = paste("OVENlike w/", sum(isGL & !isProb), "GL,",
                  sum(!isGL & isProb), "probs,",
                  sum(isGL & isProb), "both, and", sum(isTelemetry), "GPS"),
     col = RColorBrewer::brewer.pal(nNonBreeding, "Dark2"))

################################################################
# Example 4 (add prob animals released on other end)
################################################################
nAnimals <- 45
captured <- rep(c("origin", "target"), c(40, 5))
isGL <- c(OVENdata$isGL, rep(FALSE, 6))
isTelemetry <- c(!OVENdata$isGL, rep(FALSE, 6))
isRaster <- rep(FALSE, nAnimals)
isProb <- rep(FALSE, nAnimals)
targetPoints <- rbind(OVENdata$targetPoints, OVENdata$targetPoints[c(1:3,19,23,31),])
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
system.time(psi4 <-
              estTransition(isGL = isGL, isRaster = isRaster,
                            isProb = isProb,
                            isTelemetry = isTelemetry,
                            geoBias = OVENdata$geo.bias,
                            geoVCov = OVENdata$geo.vcov,
                            targetPoints = targetPoints,
                            targetAssignment = targetAssignment,
                            targetSites = targetSites,
                            resampleProjection = resampleProjection,
                            nSim = 15000, maxTries = 300,
                            originSites = originSites,
                            originAssignment = originAssignment,
                            captured = captured,
                            originNames = OVENdata$originNames,
                            targetNames = OVENdata$targetNames,
                            verbose = 2,
                            nSamples = nSamplesTry,
                            targetRelAbund = c(0.1432, 0.3577, 0.4991)))
psi4

plot(psi4, legend = "top",
     main = paste(sum(isGL & !isProb), "GL,",
                  sum(!isGL & isProb & captured == "origin"), "prob,",
                  sum(isGL & isProb), "both,",
                  sum(isTelemetry), "GPS (all\ncaptured origin), and",
                  sum(isProb & captured == "target"), "probs (captured target)"),
     col = RColorBrewer::brewer.pal(nNonBreeding, "Dark2"))
MC4 <- estStrength(OVENdata$originDist, OVENdata$targetDist,
                                     OVENdata$originRelAbund, psi4,
                                     sampleSize = nAnimals)
MC4

###############################################################################
# Example 5 (all raster, from our OVEN example)
###############################################################################
getCSV <- function(filename) {
  tmp <- tempdir()
  url1 <- paste0('https://github.com/SMBC-NZP/MigConnectivity/blob/master/data-raw/',
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
  url1 <- paste0('https://github.com/SMBC-NZP/MigConnectivity/blob/master/data-raw/Spatial_Layers/',
                 filename, '?raw=true')
  temp <- paste(tmp, filename, sep = '/')
  utils::download.file(url1, temp, mode = 'wb')
  shp <- readRDS(temp)
  unlink(temp)
  return(shp)
}
OVENdist <- getRDS("OVENdist")

OVENdist <- sf::st_as_sf(OVENdist)

OVENdist <- sf::st_transform(OVENdist, 4326)

OVENvals <- getCSV("deltaDvalues.csv")

OVENvals <- OVENvals[grep(x=OVENvals$Sample,"NH", invert = TRUE),]

originSites <- getRDS("originSites")
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
                             crs = sf::st_crs(originSites))

iso <- isoAssign(isovalues = OVENvals[,2],
                 isoSTD = 12,       # this value is for demonstration only
                 intercept = -10,   # this value is for demonstration only
                 slope = 0.8,       # this value is for demonstration only
                 odds = NULL,
                 restrict2Likely = FALSE,
                 nSamples = 1000,
                 sppShapefile = OVENdist,
                 assignExtent = c(-179,-60,15,89),
                 element = "Hydrogen",
                 surface = FALSE,
                 period = "Annual",
                 seed = 12345,
                 verbose=1)

nAnimals <- dim(iso$probassign)[3]
isGL <-rep(FALSE, nAnimals); isRaster <- rep(TRUE, nAnimals)
isProb <- rep(FALSE, nAnimals); isTelemetry <- rep(FALSE, nAnimals)
targetSites <- sf::st_as_sf(iso$targetSites)
targetSites <- sf::st_make_valid(targetSites)
targetSites <- sf::st_union(targetSites, by_feature = TRUE)


system.time(psi5 <-
              estTransition(isGL = isGL,
                            isRaster = isRaster,
                            isProb = isProb,
                            isTelemetry = isTelemetry,
                            targetSites = targetSites,
                            resampleProjection = resampleProjection,
                            targetRaster = iso,
                            originSites = originSites,
                            originPoints = originPoints,
                            captured = rep("origin", nAnimals),
                            verbose = 2,
                            nSamples = nSamplesTry))
psi5
}
