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
    results <- estStrength(originRelAbund = originRelAbundTrue, psi = fm,
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
  nSimulationsAbund <- 10
  #\dontrun{
  #  nSamplesAbund <- 1700
  #}
  # Storage matrix for samples
  abundMCSample <- matrix(NA, nSamplesAbund, nSimulationsAbund)
  summaryAbund <- data.frame(Simulation = 1:nSimulationsAbund, True = trueMC,
                             mean = NA, se = NA, lcl = NA, ucl = NA)
  for (r in 1:nSimulationsAbund) {
    cat("Simulation",r,"of",nSimulationsAbund,"\n")
    row0 <- nrow(abundExamples[[r]]) - nSamplesAbund
    results <- estStrength(originRelAbund = abundExamples[[r]], psi = psiTrue,
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

  nSamplesGLGPS <- 100 # Number of bootstrap iterations, set low for example

  # Estimate transition probabilities
  Combined.psi<-estTransition(isGL=OVENdata$isGL, #Logical vector:light-level GL(T)/GPS(F)
                              isTelemetry = !OVENdata$isGL,
                  geoBias = OVENdata$geo.bias, # Light-level GL location bias
                  geoVCov = OVENdata$geo.vcov, # Location covariance matrix
                  targetSites = OVENdata$targetSites, # Non-breeding target sites
                  originSites = OVENdata$originSites, # Breeding origin sites
                  originPoints = OVENdata$originPoints, # Capture Locations
                  targetPoints = OVENdata$targetPoints, # Device target locations
                  verbose = 3,   # output options
                  nSamples = nSamplesGLGPS, # This is set low for example
                  resampleProjection = sf::st_crs(OVENdata$targetPoints),
                  #targetAssignment = targetAssignment[1:39,],
                  nSim = 1000)

  # Can estimate MC from previous psi estimate
  Combo.MC1 <- estStrength(targetDist = OVENdata$targetDist, # targetSites distance matrix
                           originDist = OVENdata$originDist, # originSites distance matrix
                           targetSites = OVENdata$targetSites, # Non-breeding target sites
                           originSites = OVENdata$originSites, # Breeding origin sites
                           psi = Combined.psi,
                           originRelAbund = OVENdata$originRelAbund,
                           nSamples = nSamplesGLGPS,
                           sampleSize = nrow(OVENdata$targetPoints))
  Combo.MC1

  # Doesn't have to be an estPsi object - can simply be array of psi samples
  Combo.MC2 <- estStrength(targetDist = OVENdata$targetDist, # targetSites distance matrix
                           originDist = OVENdata$originDist, # originSites distance matrix
                           targetSites = OVENdata$targetSites, # Non-breeding target sites
                           originSites = OVENdata$originSites, # Breeding origin sites
                           psi = Combined.psi$psi$sample,
                           originRelAbund = OVENdata$originRelAbund,
                           nSamples = nSamplesGLGPS,
                           sampleSize = nrow(OVENdata$targetPoints))
  Combo.MC2
}
