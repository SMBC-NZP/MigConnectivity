set.seed(101)
# Uncertainty in detection with equal abundances
nSamplesCMR <- 100 #10000 # Number of resampling iterations for generating confidence intervals
nSimulationsCMR <- 10 #length(cmrExamples)
\dontrun{
  nSamplesCMR <- 10000 # Number of resampling iterations for generating confidence intervals
  nSimulationsCMR <- length(cmrExamples)
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
  fm <- cmrExamples[[r]]
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
nSamplesGLGPS <- 1000 # Number of bootstrap iterations
\dontrun{
  nSamplesGLGPS <- 10000 # Number of bootstrap iterations
}

# Estimate MC only, treat all data as is
glGPSMC <- estMC(isGL = OVENdata$isGL,
                 geoBias = OVENdata$geo.bias,
                 geoVCov = OVENdata$geo.vcov,
                 targetSites = OVENdata$targetSites,
                 targetPoints = OVENdata$targetPoints,
                 targetDist = OVENdata$targetDist,
                 originPoints = OVENdata$originPoints,
                 originSites = OVENdata$originSites,
                 originDist = OVENdata$originDist,
                 originRelAbund = OVENdata$originRelAbund,
                 nSamples = nSamplesGLGPS,
                 verbose=1,
                 calcCorr=FALSE)

# Estimate MC and correlation, treat all data as is
glGPSMCCorr <- estMC(isGL = OVENdata$isGL,
                     geoBias = OVENdata$geo.bias,
                     geoVCov = OVENdata$geo.vcov,
                     targetSites = OVENdata$targetSites,
                     targetPoints = OVENdata$targetPoints,
                     targetDist = OVENdata$targetDist,
                     originPoints = OVENdata$originPoints,
                     originSites = OVENdata$originSites,
                     originDist = OVENdata$originDist,
                     originRelAbund = OVENdata$originRelAbund,
                     nSamples = nSamplesGLGPS,
                     verbose=1,
                     calcCorr=TRUE)

# Estimate MC and correlation, treat GPS data as GL
glMCCorr <- estMC(isGL = TRUE,
                  geoBias = OVENdata$geo.bias,
                  geoVCov = OVENdata$geo.vcov,
                  targetSites = OVENdata$targetSites,
                  targetPoints = OVENdata$targetPoints,
                  targetDist = OVENdata$targetDist,
                  originPoints = OVENdata$originPoints,
                  originSites = OVENdata$originSites,
                  originDist = OVENdata$originDist,
                  originRelAbund = OVENdata$originRelAbund,
                  nSamples = nSamplesGLGPS,
                  verbose=1,
                  calcCorr=TRUE)
str(glGPSMC)
str(glGPSMCCorr)
str(glMCCorr)
# Estimate MC and correlation, treat GPS data as GL
gpsMCCorr <- estMC(isGL = FALSE,
                  geoBias = OVENdata$geo.bias,
                  geoVCov = OVENdata$geo.vcov,
                  targetSites = OVENdata$targetSites,
                  targetPoints = OVENdata$targetPoints,
                  targetDist = OVENdata$targetDist,
                  originPoints = OVENdata$originPoints,
                  originSites = OVENdata$originSites,
                  originDist = OVENdata$originDist,
                  originRelAbund = OVENdata$originRelAbund,
                  nSamples = nSamplesGLGPS,
                  verbose=1,
                  calcCorr=TRUE)

glGPSMCCorr0 <- estMC(isGL = OVENdata$isGL,
                     geoBias = rep(0, 2),
                     geoVCov = OVENdata$geo.vcov,
                     targetSites = OVENdata$targetSites,
                     targetPoints = OVENdata$targetPoints,
                     targetDist = OVENdata$targetDist,
                     originPoints = OVENdata$originPoints,
                     originSites = OVENdata$originSites,
                     originDist = OVENdata$originDist,
                     originRelAbund = OVENdata$originRelAbund,
                     nSamples = nSamplesGLGPS,
                     verbose=1,
                     calcCorr=TRUE)

# Estimate MC and correlation, treat GPS data as GL
glMCCorr0 <- estMC(isGL = TRUE,
                  geoBias = rep(0, 2),
                  geoVCov = OVENdata$geo.vcov,
                  targetSites = OVENdata$targetSites,
                  targetPoints = OVENdata$targetPoints,
                  targetDist = OVENdata$targetDist,
                  originPoints = OVENdata$originPoints,
                  originSites = OVENdata$originSites,
                  originDist = OVENdata$originDist,
                  originRelAbund = OVENdata$originRelAbund,
                  nSamples = nSamplesGLGPS,
                  verbose=1,
                  calcCorr=TRUE)
str(glGPSMCCorr0)
str(glMCCorr0)
# Estimate MC and correlation, treat GPS data as GL
gpsMCCorr0 <- estMC(isGL = FALSE,
                   geoBias = rep(0, 2),
                   geoVCov = OVENdata$geo.vcov,
                   targetSites = OVENdata$targetSites,
                   targetPoints = OVENdata$targetPoints,
                   targetDist = OVENdata$targetDist,
                   originPoints = OVENdata$originPoints,
                   originSites = OVENdata$originSites,
                   originDist = OVENdata$originDist,
                   originRelAbund = OVENdata$originRelAbund,
                   nSamples = nSamplesGLGPS,
                   verbose=1,
                   calcCorr=TRUE)
str(gpsMCCorr)
