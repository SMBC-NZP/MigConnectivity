if (interactive()){
  ### Uncertainty in detection with equal abundances
  # Number of resampling iterations for generating confidence intervals (set
  # low for example)
  nSamplesCMR <- 100
  nSimulationsCMR <- 10 # Set low for example-try length(cmrExamples) to run all
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

  ### Simulation of BBS data to quantify uncertainty in relative abundance
  nSamplesAbund <- 700 #1700 are stored
  nSimulationsAbund <- 10
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
}
