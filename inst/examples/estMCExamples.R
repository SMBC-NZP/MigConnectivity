\dontrun{
set.seed(101)
# Uncertainty in detection (RMark estimates) with equal abundances
# Number of resampling iterations for generating confidence intervals

nSamplesCMR <- 100
nSimulationsCMR <- 10
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
                         mean=NA, se=NA, lcl=NA, ucl=NA)
# Get RMark psi estimates and estimate MC from each
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
  cmrMCSample[ , r] <- results$MC$sample
  summaryCMR$mean[r] <- results$MC$mean
  summaryCMR$se[r] <- results$MC$se
  # Calculate confidence intervals using quantiles of sampled MC
  summaryCMR[r, c('lcl', 'ucl')] <- results$MC$simpleCI
}

summaryCMR <- transform(summaryCMR, coverage = (True>=lcl & True<=ucl))
summaryCMR
summary(summaryCMR)
biasCMR <- mean(summaryCMR$mean) - trueMC
biasCMR
mseCMR <- mean((summaryCMR$mean - trueMC)^2)
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
                           mean = NA, se = NA, lcl = NA, ucl = NA)
for (r in 1:nSimulationsAbund) {
  cat("Simulation",r,"of",nSimulationsAbund,"\n")
  row0 <- nrow(abundExamples[[r]]) - nSamplesAbund
  results <- estMC(originRelAbund = abundExamples[[r]], psi = psiTrue,
                   originDist = originDist, targetDist = targetDist,
                   row0 = row0, nSamples = nSamplesAbund, verbose = 1)
  abundMCSample[ , r] <- results$MC$sample
  summaryAbund$mean[r] <- results$MC$mean
  summaryAbund$se[r] <- results$MC$se
  # Calculate confidence intervals using quantiles of sampled MC
  summaryAbund[r, c('lcl', 'ucl')] <- results$MC$simpleCI
}

summaryAbund <- transform(summaryAbund, coverage = (True >= lcl & True <= ucl))
summaryAbund
summary(summaryAbund)
biasAbund <- mean(summaryAbund$mean) - trueMC
biasAbund
mseAbund <- mean((summaryAbund$mean - trueMC)^2)
mseAbund
rmseAbund <- sqrt(mseAbund)
rmseAbund

# Ovenbird example with GL and GPS data
data(OVENdata) # Ovenbird

nSamplesGLGPS <- 100 # Number of bootstrap iterations

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
             resampleProjection = terra::crs(OVENdata$targetSites))

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
                resampleProjection = terra::crs(OVENdata$targetSites),
                originNames = OVENdata$originNames,
                targetNames = OVENdata$targetNames)

print(Combined)

# For treating all data as GPS,
# Move the latitude of birds with locations that fall offshore - only change
# Latitude
int <- sf::st_intersects(OVENdata$targetPoints, OVENdata$targetSites)
any(lengths(int)<1)
plot(OVENdata$targetPoints)
plot(OVENdata$targetSites,add=TRUE)
tp<-sf::st_coordinates(OVENdata$targetPoints)
text(tp[,1], tp[,2], label=c(1:39))

tp[5,2]<- -1899469
tp[10,2]<- -1927848
tp[1,2]<- -1927930
tp[11,2]<- -2026511
tp[15,2]<- -2021268
tp[16,2]<- -1976063

oven_targetPoints<-sf::st_as_sf(as.data.frame(tp),
                                coords = c("X","Y"),
                                crs = sf::st_crs(OVENdata$targetPoints))
inter <- sf::st_intersects(oven_targetPoints, OVENdata$targetSites)
any(lengths(inter)<1)
plot(oven_targetPoints,add=TRUE, col = "green")

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

str(GPS_mc, max.level = 2)
str(Combined, max.level = 2)
str(GL_mc, max.level = 2)
if (length(find.package("RColorBrewer", quiet = TRUE))==0)
  install.packages(c('RColorBrewer'))
plot(Combined, col = RColorBrewer::brewer.pal(3, "Dark2"), legend = "top",
     main = "Ovenbird GL and GPS")
text(1.1, 0.98, cex = 1,
     labels = paste("MC = ", round(Combined$MC$mean, 2), "+/-",
                    round(Combined$MC$se, 2)))


# Generate probabilistic assignments using intrinsic markers (stable-hydrogen isotopes)
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
  url1 <- paste0(
    'https://github.com/SMBC-NZP/MigConnectivity/blob/noSpWeightMantel/data-raw/Spatial_Layers/',
                 filename, '?raw=true')
  temp <- paste(tmp, filename, sep = '/')
  utils::download.file(url1, temp, mode = 'wb')
  shp <- readRDS(temp)
  unlink(temp)
  return(shp)
}
OVENdist <- getRDS("OVENdist")

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
originDist <- distFromPos(sf::st_coordinates(op))


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

ovenMC <- estMC(originRelAbund = originRelAbund,
                targetIntrinsic = iso,
                originPoints = originPoints,
                originSites = originSites,
                originDist = originDist,
                nSamples = 200,
                verbose = 1,
                calcCorr = TRUE,
                alpha = 0.05,
                approxSigTest = FALSE,
                isIntrinsic = TRUE)

ovenMC
}
