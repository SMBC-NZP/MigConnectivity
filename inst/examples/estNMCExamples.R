\donttest{
  set.seed(101)
  # Uncertainty in detection (RMark estimates)
  # Number of resampling iterations for generating confidence intervals
  nSamplesCMR <- 100
  nSimulationsCMR <- 10

  # the second intermediate psi scenario, the "low" level
  psiTrue <- samplePsis[["Low"]]
  trueNMC <- calcNMC(psiTrue)
  trueNMC

  # Storage matrix for samples
  cmrNMCSample <- matrix(NA, nSamplesCMR, nSimulationsCMR)
  summaryCMR <- data.frame(Simulation = 1:nSimulationsCMR, True=trueNMC$NMC,
                           mean=NA, se=NA, lcl=NA, ucl=NA)
  # Get 'RMark' psi estimates and estimate MC from each
  for (r in 1:nSimulationsCMR) {
    cat("Simulation",r,"of",nSimulationsCMR,"\n")
    # Note: getCMRexample() requires a valid internet connection and that GitHub
    # is accessible
    fm <- getCMRexample(r)
    results <- estNMC(psi = fm,
                      originSites = 5:8, targetSites = c(3,2,1,4),
                      nSamples = nSamplesCMR, verbose = 0)
    cmrNMCSample[ , r] <- results$NMC$sample
    summaryCMR$mean[r] <- results$NMC$mean
    summaryCMR$se[r] <- results$NMC$se
    # Calculate confidence intervals using quantiles of sampled MC
    summaryCMR[r, c('lcl', 'ucl')] <- results$NMC$simpleCI
  }

  summaryCMR <- transform(summaryCMR, coverage = (True>=lcl & True<=ucl))
  summaryCMR
  summary(summaryCMR)
  biasCMR <- mean(summaryCMR$mean) - trueNMC$NMC
  biasCMR
  mseCMR <- mean((summaryCMR$mean - trueNMC$NMC)^2)
  mseCMR
  rmseCMR <- sqrt(mseCMR)
  rmseCMR

  # Ovenbird example with GL and GPS data
  data(OVENdata) # Ovenbird

  nSamplesGLGPS <- 100 # Number of bootstrap iterations, set low for example

  # Estimate transition probabilities
  Combined.psi<-estTransition(isGL=OVENdata$isGL, #Light-level geolocator (T/F)
                              isTelemetry = !OVENdata$isGL,
                  geoBias = OVENdata$geo.bias, # Light-level GL location bias
                  geoVCov = OVENdata$geo.vcov, # Location covariance matrix
                  targetSites = OVENdata$targetSites, # Nonbreeding/target sites
                  originSites = OVENdata$originSites, # Breeding/origin sites
                  originPoints = OVENdata$originPoints, # Capture Locations
                  targetPoints = OVENdata$targetPoints, #Device target locations
                  verbose = 3,   # output options
                  nSamples = nSamplesGLGPS, # This is set low for example
                  resampleProjection = sf::st_crs(OVENdata$targetPoints),
                  nSim = 1000)

  # Can estimate NMC from previous psi estimate
  Combo.NMC1 <- estNMC(psi = Combined.psi,
                       nSamples = nSamplesGLGPS)
  Combo.NMC1

  # Doesn't have to be an estPsi object - can simply be array of psi samples
  Combo.MC2 <- estNMC(psi = Combined.psi$psi$sample, # Array of samples
                      originNames = Combined.psi$input$originNames,
                      nSamples = nSamplesGLGPS)
  Combo.MC2
}
