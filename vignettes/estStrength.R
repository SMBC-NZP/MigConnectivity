## ----setup, include = FALSE------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE----------------------------------------------------------------------------
#  install.packages('devtools')
#  
#  devtools::install_github("SMBC-NZP/MigConnectivity", build_vignettes = TRUE)
#  
#  library(MigConnectivity)

## ----warning=FALSE,message=FALSE,echo=FALSE--------------------------------------------------
library(MigConnectivity)
library(sf)
set.seed(123)

## --------------------------------------------------------------------------------------------
nBreeding <- 4 #number of breeding regions
nNonBreeding <- 4 #number of non-breeding regions

## --------------------------------------------------------------------------------------------
originPos <- matrix(c(seq(-99, -81, 6),
                    rep(35, nBreeding)), nBreeding, 2)
targetPos <- matrix(c(seq(-93, -87, 2),
                    rep(9, nNonBreeding)), nNonBreeding, 2)

#distances between centroids of four breeding regions
breedDist <- distFromPos(originPos, surface = "ellipsoid")

#distances between centroids of four nonbreeding regions
nonBreedDist <- distFromPos(targetPos, surface = "ellipsoid")

## ----echo=FALSE------------------------------------------------------------------------------
op <- graphics::par(no.readonly = TRUE)
on.exit(graphics::par(op))
par(bty="n")
plot(originPos[,1],originPos[,2],
     ylim=c(8.5,36),
     xlim=c(-100, -80),
     pch=19,
     axes=FALSE,
     yaxt="n",
     xaxt="n",
     cex=2,
     ylab="",
     xlab="Nonbreeding",
     main="Breeding")
points(targetPos[,1],targetPos[,2],pch=15,cex=2)

## --------------------------------------------------------------------------------------------
#transition probabilities form each breeding to each nonbreeding region
psiTrue <- samplePsis[["Low"]]

psiEst <- getCMRexample(1)

# Define total sample size for psi data 
# for small sample corrected version of MC 

n <- length(grep('[2-5]', psiEst$data$data$ch))


## ----echo=FALSE------------------------------------------------------------------------------
plot(originPos[,1],originPos[,2],
     ylim=c(8.5,36),
     xlim=c(-100, -80),
     pch=19,
     axes=FALSE,
     yaxt="n",
     xaxt="n",
     cex=2,
     ylab="",
     xlab="Nonbreeding",
     main="Breeding")
points(targetPos[,1],targetPos[,2],pch=15,cex=2)

for (b in 1:nBreeding) {
  for (nb in 1:nNonBreeding) {
    segments(originPos[b, 1], originPos[b, 2], 
             targetPos[nb, 1], targetPos[nb, 2],
             lwd = psiTrue[b, nb]*10, lty = b)
  }
}

## --------------------------------------------------------------------------------------------
#equal relative abundance among the four breeding regions, must sum to 1                      
relNTrue <- rep(1/nBreeding, nBreeding) 

relNEst <- abundExamples[[1]]


## --------------------------------------------------------------------------------------------
# Calculate the strength of migratory connectivity 
MCtrue <- calcMC(originDist = breedDist,
                 targetDist = nonBreedDist,
                 originRelAbund = relNTrue,
                 psi = psiTrue)


## --------------------------------------------------------------------------------------------
MCtrue

## --------------------------------------------------------------------------------------------
MCest <- estStrength(originDist = breedDist,
             targetDist = nonBreedDist,
             originRelAbund = relNEst, 
             psi = psiEst,
             sampleSize = n,
             originSites = 5:8, 
             targetSites = c(3,2,1,4),
             nSamples = 1000)

## --------------------------------------------------------------------------------------------
MCest
MCest$MC$mean - MCtrue

## --------------------------------------------------------------------------------------------
# Read in the processed Yellow Warbler data 
# see the Worked Examples Vignette for details on how data structured/created
newDir <- tempdir()
baseURL <- 'https://github.com/SMBC-NZP/MigConnectivity/blob/devpsi2/data-raw/YEWA/'
file.name <- "yewa_estTrans_data.rds"
url1 <- paste0(baseURL, file.name, '?raw=true')
temp <- paste(newDir, file.name, sep = '/')
utils::download.file(url1, temp, mode = 'wb')
YEWA_target_sites <- readRDS(temp)
YEWAdata <- readRDS(temp)
unlink(temp)

# take a quick look at the data # 
str(YEWAdata,1)

## ----eval = FALSE----------------------------------------------------------------------------
#  ## Run analysis for psi
#  psiYEWA <- estTransition(originSites = YEWAdata$originSites,
#                           targetSites = YEWAdata$targetSites,
#                           originPoints = YEWAdata$originPoints,
#                           targetPoints = YEWAdata$targetPoints,
#                           originAssignment = YEWAdata$originAssignment,
#                           originNames = YEWAdata$originNames,
#                           targetNames = YEWAdata$targetNames,
#                           nSamples = 10,  # Set low for illustration
#                           isGL = YEWAdata$isGL,
#                           isTelemetry = YEWAdata$isTelemetry,
#                           isRaster = YEWAdata$isRaster,
#                           isProb = YEWAdata$isProb,
#                           captured = YEWAdata$captured,
#                           geoBias = YEWAdata$geoBias,
#                           geoVCov = YEWAdata$geoVCov,
#                           originRaster = YEWAdata$originRaster,
#                           verbose = 0,
#                           maxTries = 400,
#                           resampleProjection = YEWAdata$resampleProjection,
#                           nSim = 30,
#                           dataOverlapSetting = "none",
#                           targetRelAbund = YEWAdata$targetRelAbund,
#                           returnAllInput = FALSE)

## --------------------------------------------------------------------------------------------
# calculate distance between originSites
originCentersYEWA <- st_centroid(YEWAdata$originSites)
originCentersYEWA <- st_transform(originCentersYEWA, 4326)
originDistYEWA <- distFromPos(st_coordinates(originCentersYEWA$geometry))

# calculate distance between targetSites 
targetCentersYEWA <- st_centroid(YEWAdata$targetSites)
targetCentersYEWA <- st_transform(targetCentersYEWA, 4326)
targetDistYEWA <- distFromPos(st_coordinates(targetCentersYEWA$geometry))

## ----echo = FALSE----------------------------------------------------------------------------
newDir <- tempdir()
baseURL <- 'https://github.com/SMBC-NZP/MigConnectivity/blob/devpsi2/data-raw/YEWA/'
file.name <- "psiYEWA.rds"
url1 <- paste0(baseURL, file.name, '?raw=true')
temp <- paste(newDir, file.name, sep = '/')
utils::download.file(url1, temp, mode = 'wb')
psiYEWA <- readRDS(temp)
unlink(temp)


## --------------------------------------------------------------------------------------------
# Estimate the strength of migratory connectivity
YEWAmc <- estStrength(originDist = originDistYEWA,
                      targetDist = targetDistYEWA,
                      originRelAbund = YEWAdata$originRelAbund,
                      psi = psiYEWA,
                      nSamples = 5000,
                      verbose = 0)
YEWAmc

## ----echo = FALSE----------------------------------------------------------------------------
#saveRDS(YEWAmc, "YEWAmc.rds")

