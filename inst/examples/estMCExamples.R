set.seed(101)
# Uncertainty in detection with equal abundances
# Number of resampling iterations for generating confidence intervals
nSamplesCMR <- 100 #10000
nSimulationsCMR <- 10 #100
\dontrun{
  nSamplesCMR <- 10000
  nSimulationsCMR <- 100
}
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
trueMC <- calcMC(originDist, targetDist, psiTrue, originRelAbundTrue)
trueMC

# Storage matrix for samples
cmrMCSample <- matrix(NA, nSamplesCMR, nSimulationsCMR)
summaryCMR <- data.frame(Simulation = 1:nSimulationsCMR, True=trueMC, estimate=NA,
                        mean=NA, median=NA, se=NA, lcl.simple=NA, ucl.simple=NA,
                        lcl.BC=NA, ucl.BC=NA)
for (r in 1:nSimulationsCMR) {
  cat("Simulation",r,"of",nSimulationsCMR,"\n")
  # Note: getCMRexample requires a valid internet connection and that GitHub is accessible
  fm <- getCMRexample(r)
  results <- estMC(originRelAbund = originRelAbundTrue, psi = fm,
                   originDist = originDist, targetDist = targetDist,
                   originSites = 5:8, targetSites = c(3,2,1,4),
                   nSamples = nSamplesCMR, verbose = 0)
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
\dontrun{
  nSamplesAbund <- 1700
  nSimulationsAbund <- length(abundExamples)
}
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
nSamplesGLGPS <- 500 # Number of bootstrap iterations
\dontrun{
  nSamplesGLGPS <- 10000 # Number of bootstrap iterations
}

# Estimate MC only, treat all data as geolocator
GL_mc<-estMC(isGL=TRUE, # Logical vector indicating light-level geolocator (TRUE) or GPS (F)
             geoBias = OVENdata$geo.bias, # Light-level geolocator location bias
             geoVCov = OVENdata$geo.vcov, # Light-level geolocator co-variance matrix
             targetDist = OVENdata$targetDist, # Non-breeind target distribution distance matrix
             originDist = OVENdata$originDist, # Breeding / origin distribution distance matrix
             targetSites = OVENdata$targetSites, # Non-breeding target sites
             originSites = OVENdata$originSites, # Breeding origin sites
             originPoints = OVENdata$originPoints, # Capture Locations
             targetPoints = OVENdata$targetPoints, # Non-breeding Locations derived from devices
             originRelAbund = OVENdata$originRelAbund, # Relative abundance within OriginSites
             verbose = 1,   # output options
             nSamples = nSamplesGLGPS) # This is set low for example

# Estimate MC and rM, treat all data as is
Combined<-estMC(isGL=OVENdata$isGL, # Logical vector for light-level geolocator (TRUE) or GPS (F)
                geoBias = OVENdata$geo.bias, # Light-level geolocator location bias
                geoVCov = OVENdata$geo.vcov, # Light-level geolocator co-variance matrix
                targetDist = OVENdata$targetDist, # Non-breeind target distribution distance matrix
                originDist = OVENdata$originDist, # Breeding / origin distribution distance matrix
                targetSites = OVENdata$targetSites, # Non-breeding target sites
                originSites = OVENdata$originSites, # Breeding origin sites
                originPoints = OVENdata$originPoints, # Capture Locations
                targetPoints = OVENdata$targetPoints, # Non-breeding Locations derived from devices
                originRelAbund = OVENdata$originRelAbund, # Relative abundance within OriginSites
                verbose = 1,   # output options
                calcCorr = TRUE, # estimate rM as well
                nSamples = nSamplesGLGPS) # This is set low for example

# For treating all data as GPS,
# Move the latitude of birds with locations that fall off shore - only change Latitude Estimate #
tp<-OVENdata$targetPoints@coords
sp::plot(OVENdata$targetPoints)
sp::plot(OVENdata$targetSites,add=TRUE)
text(OVENdata$targetPoints@coords[,1],OVENdata$targetPoints@coords[,2],label=c(1:39))

tp[5,2]<- -1899469
tp[10,2]<- -2007848
tp[1,2]<- -2017930
tp[11,2]<- -2136511
tp[15,2]<- -2121268
tp[16,2]<- -2096063

oven_targetPoints<-sp::SpatialPoints(cbind(tp[,1],tp[,2]))
raster::crs(oven_targetPoints)<-raster::crs(OVENdata$targetPoints)

# Estimate MC only, treat all data as GPS
GPS_mc<-estMC(isGL=FALSE, # Logical vector indicating light-level geolocator (TRUE) or GPS (F)
              geoBias = OVENdata$geo.bias, # Light-level geolocator location bias
              geoVCov = OVENdata$geo.vcov, # Light-level geolocator co-variance matrix
              targetDist = OVENdata$targetDist, # Non-breeind target distribution distance matrix
              originDist = OVENdata$originDist, # Breeding / origin distribution distance matrix
              targetSites = OVENdata$targetSites, # Non-breeding target sites
              originSites = OVENdata$originSites, # Breeding origin sites
              originPoints = OVENdata$originPoints, # Capture Locations
              targetPoints = oven_targetPoints, # Non-breeding Locations derived from devices
              originRelAbund = OVENdata$originRelAbund, # Relative abundance within OriginSites
              verbose = 1,   # output options
              nSamples = nSamplesGLGPS) # This is set low for example

str(GPS_mc)
str(Combined)
str(GL_mc)
