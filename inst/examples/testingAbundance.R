library(sf)
library(MigConnectivity)
library(methods)
library(VGAM)
# Number of populations
nOriginSites <- 3
nTargetSites <- 4
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
ni. <- 1000 # reduced for example speed
# Number of iterations to thin from posterior
nt. <- 1
# Number of iterations to discard as burn-in
nb. <- 500 # reduced for example speed
# Number of MCMC chains
nc. <- 1 # reduced for example speed

# Simulation ---
# Simulate data
sim_data <- simCountData(nPops = nOriginSites, routePerPop = routePerPop.,
                         nYears = nYears., alphaPop = alphaPop.,
                         sdRoute = sdRoute., sdYear = sdYear.)

originRelAbund <- c(1/3, 1/3, 1/3)

psiTrue <- array(0, c(nOriginSites, nTargetSites))
psiTrue[1,] <- c(0.10608 + 0.10856, 0.52536, 0.16, 0.10)
psiTrue[2,] <- c(0.04441 + 0.36664, 0.30895, 0.17, 0.11)
psiTrue[3,] <- c(0.10, 0.15, 0.42463, 0.32537)
rowSums(psiTrue)


rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

sim_data2 <- simCountData(nPops = nTargetSites, routePerPop = routePerPop.,
                          nYears = nYears., alphaPop = rev$targetRelAbund * 2,
                          sdRoute = sdRoute., sdYear = sdYear.)
# Estimate population-level abundance
out_mcmc <- modelCountDataJAGS(count_data = sim_data, ni = ni., nt = nt.,
                               nb = nb., nc = nc.)
# Estimate winter abundance
out_mcmc2 <- modelCountDataJAGS(count_data = sim_data2, ni = ni., nt = nt.,
                                nb = nb., nc = nc.)

originDist <- matrix(c(0, 624.3587, 1497.2287,
                       624.3587, 0, 942.2186,
                       1497.2287, 942.2186, 0),
                     nOriginSites, nOriginSites,
                     dimnames = list(originNames, originNames))
targetDist <- matrix(c(0, 1505.0301, 2256.9304, 2299.9251,
                       1505.0301, 0, 1795.6660, 1261.6037,
                       2256.9304, 1795.6660, 0, 860.5803,
                       2299.9251, 1261.6037, 860.5803, 0),
                     nTargetSites, nTargetSites,
                     dimnames = list(targetNames, targetNames))
sampleSize <- list(rep(20, nOriginSites), rep(75, nTargetSites))

MCtrue <- calcMC(originDist, targetDist, originRelAbund, psiTrue)

shapesO <- matrix(c(5.5, 0.25, 0.01,
                    0.25, 5.5, 0.01,
                    0.004, 0.006, 5.75),
                  nOriginSites, nOriginSites, byrow = T,
                  dimnames = list(originNames, originNames))
prop.table(shapesO, 1)

shapesT <- matrix(c(15, 0.24, 0.0045, 0.0045,
                    0.24, 15, 0.0045, 0.0045,
                    0.0045, 0.0068, 18, 0.15,
                    0.0045, 0.0068, 0.15, 18), 4, 4, byrow = T,
                  dimnames = list(targetNames, targetNames))
prop.table(shapesT, 1)

data1 <- simGenData(psi = psiTrue,
                    originRelAbund = originRelAbund,
                    sampleSize = sampleSize[[1]],
                    shapes = shapesT,
                    captured = "origin",
                    verbose = 0)
tt <- data1$genProbs
oa <- data1$originAssignment
ot2 <- matrix(0, length(oa), nOriginSites)
for (i in 1:length(oa))
  ot2[i, oa[i]] <- 1
data2 <- simGenData(psi = psiTrue,
                    originRelAbund = originRelAbund,
                    sampleSize = sampleSize[[2]],
                    shapes = shapesO,
                    captured = "target",
                    verbose = 0)
ot <- data2$genProbs
ta <- data2$targetAssignment
tt2 <- matrix(0, length(ta), nTargetSites)
for (i in 1:length(ta))
  tt2[i, ta[i]] <- 1
tt <- rbind(tt, tt2)
ot <- rbind(ot2, ot)
captured <- rep("origin", sum(sampleSize[[1]]))
captured <- c(captured, rep("target", sum(sampleSize[[2]])))

est1 <- estTransition(targetAssignment = tt,
                      originAssignment = ot, maxTries = 1000,
                      #originRaster = or, #
                      originNames = originNames,
                      targetNames = targetNames,
                      nSamples = 500, isGL = F,
                      isTelemetry = F,
                      isRaster = F,
                      isProb = T,
                      captured = captured,
                      nSim = 100, verbose = 0,
                      fixedZero = fixedZero,
                      targetRelAbund = out_mcmc2)
est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                    originRelAbund = out_mcmc,
                    est1, nSamples = 500)

est1 <- estTransition(targetAssignment = tt,
                      originAssignment = ot, maxTries = 1000,
                      #originRaster = or, #
                      originNames = originNames,
                      targetNames = targetNames,
                      nSamples = 500, isGL = F,
                      isTelemetry = F,
                      isRaster = F,
                      isProb = T,
                      captured = captured,
                      nSim = 100, verbose = 0,
                      fixedZero = fixedZero,
                      targetRelAbund = as.matrix(out_mcmc2))
tr <- as.matrix(out_mcmc2)
tr <- tr[, 12:15]
names(tr) <- NULL
est1 <- estTransition(targetAssignment = tt,
                      originAssignment = ot, maxTries = 1000,
                      #originRaster = or, #
                      originNames = originNames,
                      targetNames = targetNames,
                      nSamples = 500, isGL = F,
                      isTelemetry = F,
                      isRaster = F,
                      isProb = T,
                      captured = captured,
                      nSim = 100, verbose = 0,
                      fixedZero = fixedZero,
                      targetRelAbund = tr)

out_mc <- estStrength(originDist = originDist, targetDist = targetDist,
                      originRelAbund = coda::as.mcmc(out_mcmc),
                      est1, nSamples = 500)
out_mc <- estStrength(originDist = originDist, targetDist = targetDist,
                      originRelAbund = as.matrix(out_mcmc),
                      est1, nSamples = 500)
or <- as.matrix(out_mcmc)
or <- or[, 10:12]
names(or) <- NULL
out_mc <- estStrength(originDist = originDist, targetDist = targetDist,
                      originRelAbund = or,
                      est1, nSamples = 500)
