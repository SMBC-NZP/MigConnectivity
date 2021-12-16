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
library(jagsUI)
#args <- commandArgs(TRUE)
#instance <- as.integer(args[1])
instance <- 6444
set.seed(instance)

nSims <- 20
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
originRelAbund <- c(45/60, 9/60, 6/60)

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

originSites <- sf::st_read("data-raw/PABU_breeding_regions.shp")
targetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")
originSites$originSite <- factor(originNames, levels = originNames)
targetSites$targetSite <- factor(targetNames, levels = targetNames)


# looks like originSites needs to be projected
originSites <- sf::st_transform(originSites, "ESRI:102010")
targetSites <- sf::st_transform(targetSites, "ESRI:102010")
bound <- sf::st_bbox(targetSites)
rast <- raster(xmn = bound[1], ymn = bound[2], xmx = bound[3], ymx = bound[4],
               res = 50000, crs = CRS("ESRI:102010"))
s_rast <- raster::rasterize(as(targetSites, "Spatial"), rast,
                            field = targetSites$targetSite)
plot(s_rast)
s_mat <- values(s_rast, format = "matrix")
s_vec <- values(s_rast)
rasterX <- xFromCol(s_rast)
rasterY <- yFromRow(s_rast)
locPrior <- matrix(0, length(s_vec), nTargetSites)
for (i in 1:nTargetSites) {
  cells <- which(s_vec==i)
  locPrior[cells, i] <- 1/length(cells)
}
# Load in genPops
load(file = "data-raw/genPopsSim.RData")

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

psiEst7 <- array(NA, c(nOriginSites, nTargetSites, nScenarios, nSims),
                list(originNames, targetNames, scenarios, NULL))
psiCI7 <- array(NA, c(2, nOriginSites, nTargetSites, nScenarios, nSims),
               list(c("lower", "upper"), originNames, targetNames, scenarios,
                    NULL))
MCest7 <- array(NA, c(nScenarios, nSims),
               list(scenarios, NULL))
MCCI7 <- array(NA, c(2, nScenarios, nSims),
              list(c("lower", "upper"), scenarios, NULL))
# MantelEst <- array(NA, c(nScenarios, nSims),
#                list(scenarios, NULL))
# MantelCI <- array(NA, c(2, nScenarios, nSims),
#               list(c("lower", "upper"), scenarios, NULL))
# sampleSizes <- array(NA, c(nOriginSites, nTargetSites, nScenarios, nSims),
#                 list(originNames, targetNames, scenarios, NULL))

#dataStore <- vector("list", nSims)
instances <- paste0(instance, '.', 1:nSims)
filenames <- paste0("~/MC work/simConnectivity/testConnectivity5.", instances, ".RData")

if (!is.null(fixedZero)) {
  psiFixed <- matrix(NA, nOriginSites, nTargetSites)
  for (i in 1:nrow(fixedZero)) {
    psiFixed[fixedZero[i, 1], fixedZero[i, 2]] <- 0
  }
}

ncols <- ncol(s_mat)
nrows <- nrow(s_mat)
vecLoc <- 1:length(s_vec)
trueCell <- matrix(NA, length(s_vec), 2)
trueCell[,2] <- round(vecLoc / ncols + 0.499)
trueCell[,1] <- vecLoc - (trueCell[,2] - 1) * ncols
summary(trueCell)
head(trueCell)
any(duplicated(as.data.frame(trueCell)))
trueCell[340:350, ]

for (sim7 in 1:nSims) {
  cat("Simulation", sim7, "of", nSims, "at", date(), " ")
  load(filenames[sim7])
  for (sc in 1:nScenarios) { #1:nScenarios
    cat(sc, format(Sys.time(), "%H:%M:%S"), " ")
    or <- NULL
    data1 <- dataStore[[sc]]$data1
    data2 <- dataStore[[sc]]$data2
    jags.data <- list(npop = nOriginSites, ndest = nTargetSites,
                      m0 = psiFixed)
    if (!is.null(sampleSizeGL[[sc]][[1]])){
      op <- data1$originPointsObs
      tp <- st_coordinates(data1$targetPointsObs)
      oa <- data1$originAssignment
      jags.data$pop_gl <- oa
      jags.data$bias <- geoBias
      jags.data$Sigma <- geoVCov
      jags.data$nGL <- nrow(tp)
      jags.data$tpObs <- tp
      jags.data$locPrior <- locPrior
      jags.data$ncells <- length(s_vec)
      jags.data$ncols <- ncol(s_mat)
      jags.data$rasterY <- rasterY
      jags.data$rasterX <- rasterX
    }else{
      op <- NULL
      tp <- NULL
      oa <- NULL
    }
    if (!is.null(sampleSizeGeno[[sc]])){
      #tp <- rbind(tp, data2$targetPointsTrue)
      or <- data2$genRaster
      ot <- data2$genProbs
      ta <- data2$targetAssignment
      jags.data$nProb <- nrow(ot)
      jags.data$dest <- ta
      jags.data$prob <- ot
    } else{
      or <- NULL
      ot <- matrix()
    }
    jags.inits <- function()list()
    params <- 'psi'
    file <- paste0(find.package('MigConnectivity'),
                       ifelse(is.null(sampleSizeGL[[sc]][[1]]),
                              "/JAGS/multinomial_prob_origin.txt",
                              ifelse(is.null(sampleSizeGeno[[sc]]),
                                     "/JAGS/multinomial_geolocator_target.txt",
                              "/JAGS/multinomial_geolocator_target_prob_origin.txt")))
    est1 <- autojags(jags.data, jags.inits, params, file, n.chains = 3,
                     iter.increment = 1000, parallel = TRUE, DIC = FALSE,
                     verbose = FALSE)
    # maxRhat <- max(est1$Rhat$psi, na.rm = TRUE)
    # if (maxRhat < 1.1)
    #   cat("Successful convergence based on Rhat values (all < 1.1).  ")
    # else
    #   cat("**WARNING** Rhat values indicate convergence failure.  ")

    # est1 <- estTransition(originSites,
    #                       targetSites,
    #                       op,
    #                       tp,
    #                       originAssignment = ot, maxTries = 1000,
    #                       #originRaster = or, #
    #                       originNames = originNames,
    #                       targetNames = targetNames,
    #                       nSamples = 1000, isGL = isGL[[sc]],
    #                       isTelemetry = isTelemetry[[sc]],
    #                       isRaster = isRaster[[sc]],
    #                       isProb = isProb[[sc]],
    #                       captured = captured[[sc]],
    #                       geoBias = geoBias, geoVCov = geoVCov,
    #                       resampleProjection = sf::st_crs(targetSites),
    #                       nSim = 100, verbose = 0,
    #                       dataOverlapSetting = "none",
    #                       fixedZero = fixedZero)
    est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                        originRelAbund = originRelAbund,
                        est1$sims.list$psi)
    # est3 <- estMantel(tp, op, isGL[[sc]],
    #                   geoBias, geoVCov, targetSites, 10, 80, 0,
    #                   resampleProjection = sf::st_crs(targetSites),
    #                   originSites = originSites,
    #                   captured = captured[[sc]],
    #                   isTelemetry = isTelemetry[[sc]],
    #                   isRaster = isRaster[[sc]],
    #                   originRaster = or,
    #                   dataOverlapSetting = "none")
    psiEst7[,,sc,sim7] <- est1$mean$psi
    psiCI7[1,,,sc,sim7] <- est1$q2.5$psi
    psiCI7[2,,,sc,sim7] <- est1$q97.5$psi
    MCest7[sc,sim7] <- est2$MC$mean
    MCCI7[,sc,sim7] <- est2$MC$simpleCI
    # MantelEst[sc,sim] <- est3$corr$mean
    # MantelCI[,sc,sim] <- est3$corr$simpleCI
    #dataStore[[sim]][[sc]] <- list(data1 = data1, data2 = data2)
    save(psiEst7, psiCI7, MCest7, MCCI7, sim7, sc, #dataStore,
         file = paste('testConnectivity7b', instance, 'RData', sep = '.'))
  }
  cat("\n")
}

