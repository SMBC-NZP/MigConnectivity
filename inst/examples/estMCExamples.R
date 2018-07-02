\dontrun{
set.seed(101)
# Uncertainty in detection with equal abundances
# Number of resampling iterations for generating confidence intervals

nSamplesCMR <- 100 #10000
nSimulationsCMR <- 10 #100
#\dontrun{
#  nSamplesCMR <- 10000
#  nSimulationsCMR <- 100
#}
originPos13 <- matrix(c(rep(seq(-99, -81, 2), each = 10),
                        rep(seq(49, 31, -2), 10)), 100, 2)
targetPos13 <- matrix(c(rep(seq(-79, -61, 2), each = 10),
                        rep(seq(9, -9, -2), 10)), 100, 2)
originPosCMR <- rowsum(originPos13, c(rep(1:2, 5, each = 5),
                                      rep(3:4, 5, each = 5))) / 25
originPosCMR
targetPosCMR <- rowsum(targetPos13, c(rep(1:2, 5, each = 5),
                                      rep(3:4, 5, each = 5))) / 25
targetPosCMR

originDist <- distFromPos(originPosCMR, 'ellipsoid')
targetDist <- distFromPos(targetPosCMR, 'ellipsoid')
originRelAbundTrue <- rep(0.25, 4)
# the second intermediate psi scenario, the "low" level
psiTrue <- samplePsis[["Low"]]
trueMC <- calcMC(originDist, targetDist, originRelAbundTrue, psiTrue)
trueMC

# Storage matrix for samples
cmrMCSample <- matrix(NA, nSamplesCMR, nSimulationsCMR)
summaryCMR <- data.frame(Simulation = 1:nSimulationsCMR, True=trueMC,
                         estimate=NA, mean=NA, median=NA, se=NA, lcl.simple=NA,
                         ucl.simple=NA, lcl.BC=NA, ucl.BC=NA)
for (r in 1:nSimulationsCMR) {
  cat("Simulation",r,"of",nSimulationsCMR,"\n")
  # Note: getCMRexample requires a valid internet connection and that GitHub is
  # accessible
  fm <- getCMRexample(r)
  results <- estMC(originRelAbund = originRelAbundTrue, psi = fm,
                   originDist = originDist, targetDist = targetDist,
                   originSites = 5:8, targetSites = c(3,2,1,4),
                   nSamples = nSamplesCMR, verbose = 0,
                   sampleSize = length(grep('[2-5]', fm$data$data$ch)))
  #sampleSize argument not really needed (big sample sizes)
  cmrMCSample[ , r] <- results$sampleMC
  summaryCMR$estimate[r] <- results$pointMC
  summaryCMR$mean[r] <- results$meanMC
  summaryCMR$median[r] <- results$medianMC
  summaryCMR$se[r] <- results$seMC
  # Calculate confidence intervals using quantiles of sampled MC
  summaryCMR[r, c('lcl.simple', 'ucl.simple')] <- results$simpleCI
  summaryCMR[r, c('lcl.BC', 'ucl.BC')] <- results$bcCI
}

summaryCMR <- transform(summaryCMR, coverage.simple = (True>=lcl.simple &
                                                         True<=ucl.simple),
                        coverage.BC=(True>=lcl.BC & True<=ucl.BC))
summaryCMR
summary(summaryCMR)
biasCMR <- mean(summaryCMR$estimate) - trueMC
biasCMR
mseCMR <- mean((summaryCMR$estimate - trueMC)^2)
mseCMR
rmseCMR <- sqrt(mseCMR)
rmseCMR


# Simulation of BBS data to quantify uncertainty in relative abundance

nSamplesAbund <- 700 #1700 are stored
nSimulationsAbund <- 10 #length(abundExamples) is 100
#\dontrun{
#  nSamplesAbund <- 1700
#  nSimulationsAbund <- length(abundExamples)
#}
# Storage matrix for samples
abundMCSample <- matrix(NA, nSamplesAbund, nSimulationsAbund)
summaryAbund <- data.frame(Simulation = 1:nSimulationsAbund, True = trueMC,
                           estimate = NA, mean = NA, median = NA, se = NA,
                           lcl.simple = NA, ucl.simple = NA,
                           lcl.BC = NA, ucl.BC = NA, lclHPD = NA, uclHPD = NA)
for (r in 1:nSimulationsAbund) {
  cat("Simulation",r,"of",nSimulationsAbund,"\n")
  row0 <- nrow(abundExamples[[r]]) - nSamplesAbund
  results <- estMC(originRelAbund = abundExamples[[r]], psi = psiTrue,
                   originDist = originDist, targetDist = targetDist,
                   row0 = row0, nSamples = nSamplesAbund, verbose = 1)
  abundMCSample[ , r] <- results$sampleMC
  summaryAbund$estimate[r] <- results$pointMC
  summaryAbund$mean[r] <- results$meanMC
  summaryAbund$median[r] <- results$medianMC
  summaryAbund$se[r] <- results$seMC
  # Calculate confidence intervals using quantiles of sampled MC
  summaryAbund[r, c('lcl.simple', 'ucl.simple')] <- results$simpleCI
  summaryAbund[r, c('lcl.BC', 'ucl.BC')] <- results$bcCI
  summaryAbund[r, c('lclHPD', 'uclHPD')] <- results$hpdCI
}

summaryAbund <- transform(summaryAbund, coverage.simple = (True >= lcl.simple &
                                                           True <= ucl.simple),
                          coverage.BC=(True>=lcl.BC & True<=ucl.BC),
                          coverage.HPD=(True>=lclHPD & True<=uclHPD))
summaryAbund
summary(summaryAbund)
biasAbund <- mean(summaryAbund$estimate) - trueMC
biasAbund
mseAbund <- mean((summaryAbund$estimate - trueMC)^2)
mseAbund
rmseAbund <- sqrt(mseAbund)
rmseAbund

# Ovenbird example with GL and GPS data

nSamplesGLGPS <- 100 # Number of bootstrap iterations
#\dontrun{
#  nSamplesGLGPS <- 10000 # Number of bootstrap iterations
#  install.packages(c('raster', 'maptools', 'rgdal', 'rgeos', 'Rcpp'))
#}

# Estimate MC only, treat all data as geolocator
GL_mc<-estMC(isGL=TRUE, # Logical vector: light-level geolocator(T)/GPS(F)
             geoBias = OVENdata$geo.bias, #Light-level geolocator location bias
             geoVCov = OVENdata$geo.vcov, # Location covariance matrix
             targetDist = OVENdata$targetDist, # targetSites distance matrix
             originDist = OVENdata$originDist, # originSites distance matrix
             targetSites = OVENdata$targetSites, # Non-breeding target sites
             originSites = OVENdata$originSites, # Breeding origin sites
             originPoints = OVENdata$originPoints, # Capture Locations
             targetPoints = OVENdata$targetPoints, # Device target locations
             originRelAbund = OVENdata$originRelAbund,#Origin relative abund.
             verbose = 1,   # output options
             nSamples = nSamplesGLGPS,# This is set low for example
             resampleProjection = raster::projection(OVENdata$targetSites))

str(GL_mc)

# Estimate MC and rM, treat all data as is
Combined<-estMC(isGL=OVENdata$isGL, #Logical vector:light-level GL(T)/GPS(F)
                geoBias = OVENdata$geo.bias, # Light-level GL location bias
                geoVCov = OVENdata$geo.vcov, # Location covariance matrix
                targetDist = OVENdata$targetDist, # targetSites distance matrix
                originDist = OVENdata$originDist, # originSites distance matrix
                targetSites = OVENdata$targetSites, # Non-breeding target sites
                originSites = OVENdata$originSites, # Breeding origin sites
                originPoints = OVENdata$originPoints, # Capture Locations
                targetPoints = OVENdata$targetPoints, # Device target locations
                originRelAbund = OVENdata$originRelAbund,#Origin relative abund
                verbose = 1,   # output options
                calcCorr = TRUE, # estimate rM as well
                nSamples = nSamplesGLGPS, # This is set low for example
                approxSigTest = TRUE,
                resampleProjection = raster::projection(OVENdata$targetSites))

# For treating all data as GPS,
# Move the latitude of birds with locations that fall off shore - only change
# Latitude Estimate #
tp<-OVENdata$targetPoints@coords
sp::plot(OVENdata$targetPoints)
sp::plot(OVENdata$targetSites,add=TRUE)
text(OVENdata$targetPoints@coords[,1], OVENdata$targetPoints@coords[,2],
     label=c(1:39))

tp[5,2]<- -1899469
tp[10,2]<- -2007848
tp[1,2]<- -2017930
tp[11,2]<- -2136511
tp[15,2]<- -2121268
tp[16,2]<- -2096063

oven_targetPoints<-sp::SpatialPoints(cbind(tp[,1],tp[,2]))
raster::crs(oven_targetPoints)<-raster::crs(OVENdata$targetPoints)

# Estimate MC only, treat all data as GPS
GPS_mc<-estMC(isGL=FALSE, # Logical vector: light-level geolocator(T)/GPS(F)
              targetDist = OVENdata$targetDist, # targetSites distance matrix
              originDist = OVENdata$originDist, # originSites distance matrix
              targetSites = OVENdata$targetSites, # Non-breeding target sites
              originSites = OVENdata$originSites, # Breeding origin sites
              originPoints = OVENdata$originPoints, # Capture Locations
              targetPoints = oven_targetPoints, # Device target locations
              originRelAbund = OVENdata$originRelAbund,#Origin relative abund.
              verbose = 1,   # output options
              nSamples = nSamplesGLGPS) # This is set low for example

str(GPS_mc)
str(Combined)
str(GL_mc)


# Generate probabilistic assignments using intrinsic markers (stable-hydrogen isotopes)

OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
raster::crs(OVENdist) <- MigConnectivity::projections$WGS84

OVENvals <- read.csv("data-raw/deltaDvalues.csv")

world <- raster::shapefile("data-raw/Spatial_Layers/TM_WORLD_BORDERS-0.3.shp")
states <- raster::shapefile("data-raw/Spatial_Layers/st99_d00.shp")

JAMrock <- world[world$NAME=="Jamaica","NAME"]
Florida <- states[states$NAME=="Florida","NAME"]
Florida <- rgeos::gUnaryUnion(Florida)
Florida$NAME <- "Florida"

originSites <- rbind(JAMrock,Florida)
originDist <- originDist <- MigConnectivity::distFromPos(rgeos::gCentroid(originSites,byid = TRUE)@coords)

EVER <- length(grep(x=OVENvals$Sample,"EVER"))
JAM <- length(grep(x=OVENvals$Sample,"JAM"))

originRelAbund <- matrix(c(EVER,JAM),nrow = 1,byrow = TRUE)
originRelAbund <- prop.table(originRelAbund,1)

op <- rgeos::gCentroid(originSites,byid = TRUE)

originPoints <- array(NA,c(EVER+JAM,2))
originPoints[grep(x = OVENvals$Sample,"JAM"),1] <- sp::coordinates(op[1])[,1]
originPoints[grep(x = OVENvals$Sample,"JAM"),2] <- sp::coordinates(op[1])[,2]
originPoints[grep(x = OVENvals$Sample,"EVER"),1] <- sp::coordinates(op[2])[,1]
originPoints[grep(x = OVENvals$Sample,"EVER"),2] <- sp::coordinates(op[2])[,2]

originPoints <- sp::SpatialPoints(originPoints)
raster::crs(originPoints)<- MigConnectivity::projections$WGS84

b <- isoAssign(isovalues = OVENvals[,2],
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

simEst <- estMC(originRelAbund = originRelAbund,
                targetIntrinsic = b,
                originPoints = originPoints,
                originSites = originSites,
                originDist = originDist,
                nSamples = 5,
                verbose=1,
                calcCorr=TRUE,
                alpha = 0.05,
                approxSigTest = F,
                sigConst = 0,
                isIntrinsic = TRUE)


# Identify weights to use for abundance and isotope values when making assignments

OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
raster::crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

OVENvals <- read.csv("data-raw/deltaDvalues.csv")

HBEFbirds <- OVENvals[grep("NH",OVENvals[,1]),]
knownLocs <- cbind(rep(-73,nrow(HBEFbirds)),rep(43,nrow(HBEFbirds)))

utils::download.file("https://www.mbr-pwrc.usgs.gov/bbs/ra15/ra06740.zip", destfile = "oven.zip")
utils::unzip("oven.zip")
oven_dist <- raster::shapefile("ra06740.shp")

# Empty raster with the same dimensions as isoscape and Ovenbird distribution
r <- raster::raster(nrow = 83, ncol = 217, res = c(0.333333, 0.333333),
                     xmn = -125.0001, xmx = -52.66679, ymn = 33.33321, ymx = 60.99985,
                     crs = MigConnectivity::projections$WGS84)

relativeAbund <- raster::rasterize(sp::spTransform(oven_dist, sp::CRS(r@crs@projargs)),r)
relativeAbund <- relativeAbund /raster::cellStats(relativeAbund ,sum)

BE <- weightAssign(knownLocs = knownLocs,
                  isovalues = HBEFbirds[,2],
                  isoSTD = 12,     # this value is for demonstration only
                  intercept = -10, # this value is for demonstration only
                  slope = 0.8,     # this value is for demonstration only
                  odds = 0.67,
                  relAbund = relativeAbund,
                  weightRange = c(-1,1),
                  sppShapefile = OVENdist,
                  assignExtent = c(-179,-60,15,89),
                  element = "Hydrogen",
                  surface = FALSE,
                  period = "Annual")}
