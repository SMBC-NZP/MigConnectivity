## Example 1: sample psis and relative abundances from Cohen et al. (2018)
## (no uncertainty in psi or relative abundance)
for (i in 1:length(samplePsis)) {
 for (j in 1:length(sampleOriginRelN)){
  cat("For psi:\n")
  print(samplePsis[[i]])
  cat("and origin relative abundance:", sampleOriginRelN[[j]], "\n")
  print(reverseTransition(samplePsis[[i]], sampleOriginRelN[[j]]))
 }
}

## Example 2: Common tern banding example (uncertainty in psi, not relative
## abundance)
# Number of MCMC iterations
ni. <- 1000 # reduced from 70000 for example speed
# Number of iterations to thin from posterior
nt. <- 1
# Number of iterations to discard as burn-in
nb. <- 500 # reduced from 20000 for example speed
# Number of MCMC chains
nc. <- 1 # reduced from 3 for example speed
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
                          nSamples = (ni. - nb.) / nt. * nc., nBurnin = nb.,
                          nThin = nt., nChains = nc.,
                          method = "MCMC")
COTE_psi
COTE_rev <- reverseTransition(COTE_psi, sampleOriginRelN[[1]],
                               nSamples = 2000)
COTE_rev

\donttest{
## Example 3: Uncertainty in both psi and relative abundance
# Number of populations
nOriginSites <- 3; originNames <- LETTERS[1:nOriginSites]
nTargetSites <- 4; targetNames <- 1:nTargetSites

originRelAbund <- c(1/3, 1/3, 1/3)

psiTrue <- array(0, c(nOriginSites, nTargetSites),
                 list(originNames, targetNames))
psiTrue[1,] <- c(0.22, 0.52, 0.16, 0.10)
psiTrue[2,] <- c(0.41, 0.31, 0.17, 0.11)
psiTrue[3,] <- c(0.10, 0.15, 0.42, 0.33)
rowSums(psiTrue)

rev <- reverseTransition(psiTrue, originRelAbund)

# Simulate abundance data on origin sites
# Number of routes w/i each population (assumed to be balanced)
routePerPop. <- 30 # reduced for example speed
# Number of years
nYears. <- 5 # reduced for example speed
# log(Expected number of birds counted at each route)
alphaPop. <- 1.95
# standard deviation of normal distribution assumed for route/observer random
# effects
sdRoute. <- 0.6
# standard deviation of normal distribution assumed for year random effects
sdYear. <- 0.18
# Number of MCMC iterations
ni. <- 1000 # reduced from 70000 for example speed
# Number of iterations to thin from posterior
nt. <- 1
# Number of iterations to discard as burn-in
nb. <- 500 # reduced from 20000 for example speed
# Number of MCMC chains
nc. <- 1 # reduced from 3 for example speed

sim_data <- simCountData(nStrata = nOriginSites, routesPerStrata = routePerPop.,
                         nYears = nYears., alphaStrat = alphaPop.,
                         sdRoute = sdRoute., sdYear = sdYear.)
# Estimate population-level abundance
out_mcmc <- modelCountDataJAGS(count_data = sim_data, ni = ni., nt = nt.,
                               nb = nb., nc = nc.)

# Simulate movement data
sampleSize <- list(rep(20, nOriginSites), NULL)
captured <- rep("origin", sum(sampleSize[[1]]))
isTelemetry <- rep(TRUE:FALSE, c(sum(sampleSize[[1]]), sum(sampleSize[[2]])))
isProb <- rep(FALSE:TRUE, c(sum(sampleSize[[1]]), sum(sampleSize[[2]])))

# Telemetry data (released origin)
data1 <- simTelemetryData(psi = psiTrue,
                          sampleSize = sampleSize[[1]],
                          captured = "origin")
tt <- data1$targetAssignment
oa <- data1$originAssignment

# Estimate transition probabilities (psi)
est1 <- estTransition(targetAssignment = tt,
                      originAssignment = oa,
                      originNames = originNames,
                      targetNames = targetNames,
                      nSamples = 500, isGL = FALSE,
                      isTelemetry = isTelemetry,
                      isRaster = FALSE,
                      isProb = isProb,
                      captured = captured,
                      nSim = 10, verbose = 0)
# Reverse estimates
rev1 <- reverseTransition(psi = est1, originRelAbund = out_mcmc)
# Compare estimates of gamma, target relative abundance, and pi with calculation
# from true values
rev
rev1
}
