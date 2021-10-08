################################################################################
# Simulation to test new estTransition/estMantel/estStrength options, including
# data with location uncertainty on origin side
# This version simulates GL data released on origin and genoscape data gathered
# on the target side
################################################################################
library(raster)
library(sf)
library(MigConnectivity)
library(methods)
args <- commandArgs(TRUE)
instance <- as.integer(args[1])
set.seed(instance)

nSims <- 20
scenarios <- c('Proportional Released Origin',
               'Proportional Released Target',
               'Proportional Released Both',
               'Equal Released Origin',
               'Equal Released Target',
               'Origin Release East; Target Release West',
               'Small Origin Release East; Target Release West')

nScenarios <- length(scenarios)
nOriginSites <- 3
nTargetSites <- 5

# psiTrue <- psiPABU$psi$mean
# Transition probability (psi) estimates (mean):
#            Pacific and Interior Mexico Atlantic Lowland Mexico Central America
# Central                        0.10608               1.086e-01       0.7853552
# Louisiana                      0.04441               4.866e-01       0.4689479
# East Coast                     0.00000               6.017e-05       0.0002171
#            Southeastern US Caribbean
# Central             0.0000    0.0000
# Louisiana           0.0000    0.0000
# East Coast          0.5444    0.4554


psiTrue <- array(0, c(nOriginSites, nTargetSites))
psiTrue[1,] <- c(0.10608, 0.10856, 0.78536, 0, 0)
psiTrue[2,] <- c(0.04441, 0.48664, 0.46895, 0, 0)
psiTrue[3,] <- c(0, 0.00006, 0.00022, 0.54435, 0.45537)
rowSums(psiTrue)

fixedZero <- which(psiTrue==0, arr.ind = TRUE)

originNames <- LETTERS[1:nOriginSites]
targetNames <- as.character(1:nTargetSites)
dimnames(psiTrue) <- list(originNames, targetNames)
geoBias <- c(22250.92, -32407.14)
geoVCov <- matrix(c(5637340653, -2692682155, -2692682155, 6012336962),
                  nrow = 2, ncol = 2)

# Changed these slightly from PABU estimates so we'd have at least two GL per
# population in first and third scenario
originRelAbund <- c(45/60, 9/60, 6/60)

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

originSites <- sf::st_read("data-raw/PABU_breeding_regions.shp")
targetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")
originSites$originSite <- originNames
targetSites$targetSite <- targetNames


# looks like originSites needs to be projected
originSites <- sf::st_transform(originSites, "ESRI:102010")
targetSites <- sf::st_transform(targetSites, "ESRI:102010")

# Load in genPops
load(file = "data-raw/genPopsSim.RData")

S <- vector("list", nScenarios)
S[[1]] <- S[[2]] <- S[[3]] <- S[[4]] <- S[[5]] <- S[[6]] <- S[[7]] <- 1
p <- vector("list", nScenarios)
for (i in 1:nScenarios)
  p[[i]] <- list(1, 1)

sampleSize <- vector("list", nScenarios)
sampleSize[[1]] <- list(60 * originRelAbund, NULL)
sampleSize[[2]] <- list(NULL, round(rev$targetRelAbund * 300))
sampleSize[[3]] <- list(60 * originRelAbund, round(rev$targetRelAbund * 300))
sampleSize[[4]] <- list(c(20, 20, 20), NULL)
sampleSize[[5]] <- list(NULL, rep(60, nTargetSites))
sampleSize[[6]] <- list(c(0, 0, 60), c(100, 100, 100, 0, 0))
sampleSize[[7]] <- list(c(0, 0, 6), c(100, 100, 100, 0, 0))
sampleSizeGL <- lapply(sampleSize, function(x) list(x[[1]], NULL))
sampleSizeGeno <- lapply(sampleSize, function(x) x[[2]])

captured <- vector("list", nScenarios)
isGL <- vector("list", nScenarios); isTelemetry <- vector("list", nScenarios)
isRaster <- vector("list", nScenarios); isProb <- vector("list", nScenarios)
for (i in 1:nScenarios) {
  if (!is.null(sampleSize[[i]][[1]])) {
    captured[[i]] <- rep("origin", sum(sampleSize[[i]][[1]]))
    isGL[[i]] <- rep(TRUE, sum(sampleSize[[i]][[1]]))
    isRaster[[i]] <- rep(FALSE, sum(sampleSize[[i]][[1]]))
    isTelemetry[[i]] <- rep(FALSE, sum(sampleSize[[i]][[1]]))
    isProb[[i]] <- rep(FALSE, sum(sampleSize[[i]][[1]]))
  }
  if (!is.null(sampleSize[[i]][[2]])){
    captured[[i]] <- c(captured[[i]], rep("target", sum(sampleSize[[i]][[2]])))
    isGL[[i]] <- c(isGL[[i]], rep(FALSE, sum(sampleSize[[i]][[2]])))
    isRaster[[i]] <- c(isRaster[[i]], rep(F, sum(sampleSize[[i]][[2]])))
    #isRaster[[i]] <- c(isRaster[[i]], rep(T, sum(sampleSize[[i]][[2]])))
    isTelemetry[[i]] <- c(isTelemetry[[i]], rep(FALSE, sum(sampleSize[[i]][[2]])))
    isProb[[i]] <- c(isProb[[i]], rep(T, sum(sampleSize[[i]][[2]])))
    #isProb[[i]] <- c(isProb[[i]], rep(F, sum(sampleSize[[i]][[2]])))
  }
}

targetCenters <- sf::st_transform(targetSites, 4326) %>% st_centroid()

targetDist <- distFromPos(st_coordinates(targetCenters$geometry))
originCenters <- st_centroid(originSites)
originCenters <- st_transform(originCenters, 4326)
originDist <- distFromPos(st_coordinates(originCenters$geometry))

MCtrue <- calcMC(originDist, targetDist, originRelAbund, psiTrue)

psiEst <- array(NA, c(nOriginSites, nTargetSites, nScenarios, nSims),
                list(originNames, targetNames, scenarios, NULL))
psiCI <- array(NA, c(2, nOriginSites, nTargetSites, nScenarios, nSims),
               list(c("lower", "upper"), originNames, targetNames, scenarios,
                    NULL))
MCest <- array(NA, c(nScenarios, nSims),
               list(scenarios, NULL))
MCCI <- array(NA, c(2, nScenarios, nSims),
              list(c("lower", "upper"), scenarios, NULL))
# MantelEst <- array(NA, c(nScenarios, nSims),
#                list(scenarios, NULL))
# MantelCI <- array(NA, c(2, nScenarios, nSims),
#               list(c("lower", "upper"), scenarios, NULL))
# sampleSizes <- array(NA, c(nOriginSites, nTargetSites, nScenarios, nSims),
#                 list(originNames, targetNames, scenarios, NULL))

dataStore <- vector("list", nSims)


for (sim in 1:nSims) {
  cat("Simulation", sim, "of", nSims, "at", date(), " ")
  dataStore[[sim]] <- vector("list", nScenarios)
  for (sc in 1:nScenarios) {
    cat(sc, format(Sys.time(), "%H:%M:%S"), " ")
    or <- NULL
    if (!is.null(sampleSizeGL[[sc]][[1]])){
      data1 <- simGLData(psi = psiTrue, originRelAbund = originRelAbund,
                     sampleSize = sampleSizeGL[[sc]],
                     originSites = originSites, targetSites = targetSites,
                     #captured = "origin",
                     geoBias, geoVCov,
                     S = S[[sc]], p = p[[sc]],
                     requireEveryOrigin = is.null(sampleSizeGeno[[sc]]),
                     verbose = 1)
      op <- data1$originPointsObs
      tp <- data1$targetPointsObs
    }else{
      data1 <- NULL
      op <- NULL
      tp <- NULL
    }
    if (!is.null(sampleSizeGeno[[sc]])){
      data2 <- simGeneticData(genPops = genPops, psi = psiTrue,
                              originRelAbund = originRelAbund,
                              sampleSize = sampleSizeGeno[[sc]],
                              originSites = originSites,
                              targetSites = targetSites,
                              captured = "target",
                              verbose = 0)
      tp <- rbind(tp, data2$targetPointsTrue)
      or <- data2$genRaster
      ot <- data2$genProbs
      #originSites <- sf::st_transform(originSites, crs(or, TRUE))
      #crs(or) <- sf::st_crs(originSites)
    }else{
      data2 <- NULL
      or <- NULL
      ot <- NULL
    }
    est1 <- estTransition(originSites,
                          targetSites,
                          op,
                          tp,
                          originAssignment = ot, maxTries = 1000,
                          #originRaster = or, #
                          originNames = originNames,
                          targetNames = targetNames,
                          nSamples = 1000, isGL = isGL[[sc]],
                          isTelemetry = isTelemetry[[sc]],
                          isRaster = isRaster[[sc]],
                          isProb = isProb[[sc]],
                          captured = captured[[sc]],
                          geoBias = geoBias, geoVCov = geoVCov,
                          resampleProjection = sf::st_crs(targetSites),
                          nSim = 100, verbose = 3,
                          dataOverlapSetting = "none",
                          fixedZero = fixedZero)
    est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                        originRelAbund = originRelAbund,
                        est1)
    # est3 <- estMantel(tp, op, isGL[[sc]],
    #                   geoBias, geoVCov, targetSites, 10, 80, 0,
    #                   resampleProjection = sf::st_crs(targetSites),
    #                   originSites = originSites,
    #                   captured = captured[[sc]],
    #                   isTelemetry = isTelemetry[[sc]],
    #                   isRaster = isRaster[[sc]],
    #                   originRaster = or,
    #                   dataOverlapSetting = "none")
    psiEst[,,sc,sim] <- est1$psi$mean
    psiCI[,,,sc,sim] <- est1$psi$simpleCI
    MCest[sc,sim] <- est2$MC$mean
    MCCI[,sc,sim] <- est2$MC$simpleCI
    # MantelEst[sc,sim] <- est3$corr$mean
    # MantelCI[,sc,sim] <- est3$corr$simpleCI
    dataStore[[sim]][[sc]] <- list(data1 = data1, data2 = data2)
    save(psiEst, psiCI, MCest, MCCI, sim, sc, dataStore,
         file = paste('testConnectivity5', instance, 'RData', sep = '.'))
  }
  cat("\n")
}
