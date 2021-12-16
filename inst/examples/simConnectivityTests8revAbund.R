library(raster)
library(sf)
library(MigConnectivity)
library(methods)
library(VGAM)
instance <- 205# as.integer(args[1])
set.seed(instance)

nSims <- 100

simGenData <- function(psi,
                       originRelAbund,
                       sampleSize,
                       shapes,
                       captured = "target",
                       requireEveryOrigin = FALSE,
                       resampleNAs = TRUE,
                       verbose = 0) {
  if (verbose > 0)
    cat("Setting up...\n")
  nOriginSites <- nrow(psi)
  nTargetSites <- ncol(psi)
  rev <- MigConnectivity:::reversePsiRelAbund(psi, originRelAbund)
  gamma <- rev$gamma
  targetRelAbund <- rev$targetRelAbund
  nAnimals <- sum(sampleSize)
  if (length(captured)>1)
    captured <- captured[1]
  targetAssignment <- rep(NA, nAnimals)
  originAssignment <- rep(NA, nAnimals)
  if (captured=="origin" && length(sampleSize)>1){
    originAssignment <- rep(1:nOriginSites, sampleSize)
    if (requireEveryOrigin){
      if (any(sampleSize < 1))
        stop("Some origin site or sites have no samples ", sampleSize)
    }
  }
  if (captured=="target" && length(sampleSize)>1){
    targetAssignment <- rep(1:nTargetSites, sampleSize)
  }

  if (verbose > 0)
    cat("Assigning sites to animals...\n")
  runWell <- FALSE
  while (!runWell){
    for (a in 1:nAnimals) {
      if (captured=="origin") {
        if (is.na(originAssignment[a]))
          originAssignment[a] <- sample.int(nOriginSites, 1, prob = originRelAbund)
        targetAssignment[a] <- sample.int(nTargetSites, 1, prob = psi[originAssignment[a], ])
      }
      else {
        if (is.na(targetAssignment[a]))
          targetAssignment[a] <- sample.int(nTargetSites, 1, prob = targetRelAbund)
        originAssignment[a] <- sample.int(nOriginSites, 1, prob = gamma[targetAssignment[a], ])
      }
    }
    runWell <- !requireEveryOrigin || captured=="target" ||
      all.equal(unique(originAssignment), 1:nOriginSites)
  }

  if (verbose > 0)
    cat("Assigning genProbs...\n")

  if (captured == "target") {

    genProbs <- array(NA, c(nAnimals, nOriginSites))
    for (a in 1:nAnimals) {
      genProbs[a, ] <- rdiric(1, shapes[originAssignment[a], ])
    }
  }
  else {
    genProbs <- array(NA, c(nAnimals, nTargetSites))
    for (a in 1:nAnimals) {
      genProbs[a, ] <- rdiric(1, shapes[targetAssignment[a], ])
    }
  }


  return(list(originAssignment = originAssignment,
              targetAssignment = targetAssignment,
              genProbs = genProbs,
              input = list(psi = psi,
                           originRelAbund = originRelAbund,
                           captured = captured,
                           shapes = shapes)))
}

scenarios <- c('Proportional Release Origin',
               'Proportional Release Target',
               'Proportional Release Both',
               'Equal Release Origin',
               'Equal Release Target',
               'Equal Release Both',
               'Origin Release East; Target Release West',
               'Small Origin Release East; Target Release West')

nScenarios <- length(scenarios)
nOriginSites <- 3
nTargetSites <- 5
fastScenarios <- c(2, 5)
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

# Changed these slightly from PABU estimates so we'd have at least six GL per
# population in first and third scenario
# UPDATE: also reversed them
originRelAbund <- c(6/60, 9/60, 45/60)

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

originSites <- sf::st_read("data-raw/PABU_breeding_regions.shp")
targetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")
originSites$originSite <- originNames
targetSites$targetSite <- targetNames


# looks like originSites needs to be projected
originSites <- sf::st_transform(originSites, "ESRI:102010")
targetSites <- sf::st_transform(targetSites, "ESRI:102010")


S <- vector("list", nScenarios)
p <- vector("list", nScenarios)
for (i in 1:nScenarios){
  p[[i]] <- list(1, 1)
  S[[i]] <- 1
}

sampleSize <- vector("list", nScenarios)
sampleSize[[1]] <- list(60 * originRelAbund, NULL)
sampleSize[[2]] <- list(NULL, round(rev$targetRelAbund * 300))
sampleSize[[3]] <- list(60 * originRelAbund, round(rev$targetRelAbund * 300))
sampleSize[[4]] <- list(c(20, 20, 20), NULL)
sampleSize[[5]] <- list(NULL, rep(60, nTargetSites))
sampleSize[[6]] <- list(rep(20, nOriginSites), rep(60, nTargetSites))
sampleSize[[7]] <- list(c(0, 0, 60), c(100, 100, 100, 0, 0))
sampleSize[[8]] <- list(c(0, 0, 6), c(100, 100, 100, 0, 0))
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


shapes <- matrix(c(5.5, 0.25, 0.01,
                   0.25, 5.5, 0.01,
                   0.004, 0.006, 5.75), 3, 3, byrow = T,
                 dimnames = list(originNames, originNames))
prop.table(shapes, 1)

for (sim in 1:nSims) {
  cat("Simulation", sim, "of", nSims, "at", date(), " ")
  dataStore[[sim]] <- vector("list", nScenarios)
  for (sc in fastScenarios) { #1:nScenarios
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
                         verbose = 0)
      op <- data1$originPointsObs
      tp <- data1$targetPointsObs
    }else{
      data1 <- NULL
      op <- NULL
      tp <- NULL
    }
    if (!is.null(sampleSizeGeno[[sc]])){
      data2 <- simGenData(psi = psiTrue,
                          originRelAbund = originRelAbund,
                          sampleSize = sampleSizeGeno[[sc]],
                          shapes = shapes,
                          captured = "target",
                          verbose = 0)
      ot <- data2$genProbs
      ta <- data2$targetAssignment
      #originSites <- sf::st_transform(originSites, crs(or, TRUE))
      #crs(or) <- sf::st_crs(originSites)
    }else{
      data2 <- NULL
      or <- NULL
      ot <- NULL
    }
    est1 <- estTransition(originSites,
                          targetSites,
                          targetAssignment = ta,
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
                          nSim = 100, verbose = 0)
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
         file = paste('testConnectivity8', instance, 'RData', sep = '.'))
  }
  cat("\n")
}
