# SIMULATIONS FOR ESTIMATING MIGRATORY CONNECTIVITY USING GENETIC DATA

################################################################################
# Load required packages
library(raster)
library(sf)
library(ks)
library(MigConnectivity)

# Set the  number of simulations
nSims <- 1000

# Establish the scenarios #
scenarios <- c('Proportional Released Origin',
               'Proportional Released Target',
               'Proportional Released Both',
               'Equal Released Origin',
               'Equal Released Target',
               'Origin Release East; Target Release West',
               'Small Origin Release East; Target Release West')

# sample sizes, etc #
nScenarios <- length(scenarios)
nOriginSites <- 3
nTargetSites <- 5

# True transition matrix

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

originNames <- LETTERS[1:nOriginSites]
targetNames <- as.character(1:nTargetSites)

dimnames(psiTrue) <- list(originNames, targetNames)

# Geolocator bias and covariance structure (units = meters)
geoBias <- c(22250.92, -32407.14)
geoVCov <- matrix(c(5637340653, -2692682155,
                    -2692682155, 6012336962),
                  nrow = 2, ncol = 2)

# Changed these from PABU estimates so we'd have at least one GL per
# population in first and third scenario
originRelAbund <- c(27/30, 2/30, 1/30)

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

originSites <- sf::st_read("data-raw/PABU_breeding_regions.shp")
targetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")

originSites$originSite <- originNames
targetSites$targetSite <- targetNames

# looks like originSites needs to be projected
originSites <- sf::st_transform(originSites, "ESRI:102010")
targetSites <- sf::st_transform(targetSites, "ESRI:102010")



# nSims <- 100

# set up simulation where equal number of birds sampled in target
sampleSize <- list(NULL, c(30, 39, 221, 5, 5))

captured <- rep("target", sum(sampleSize[[2]]))
isGL <- rep(FALSE, sum(sampleSize[[2]]))
isRaster <- rep(TRUE, sum(sampleSize[[2]]))
isTelemetry <- rep(FALSE, sum(sampleSize[[2]]))
isProb <- rep(F, sum(sampleSize[[2]]))

targetCenters <- sf::st_transform(targetSites, 4326) %>% st_centroid()

targetDist <- distFromPos(st_coordinates(targetCenters$geometry))
originCenters <- st_centroid(originSites)
originCenters <- st_transform(originCenters, 4326)
originDist <- distFromPos(st_coordinates(originCenters$geometry))

MCtrue <- calcMC(originDist, targetDist, originRelAbund, psiTrue)

psiEstRaster <- array(NA, c(nOriginSites, nTargetSites,nSims),
                      list(originNames, targetNames,NULL))
psiCIRaster <- array(NA, c(2, nOriginSites, nTargetSites,nSims),
                     list(c("lower", "upper"), originNames, targetNames,NULL))
MCestRaster <- array(NA, c(nSims),
                     list(NULL))
MCCIRaster <- array(NA, c(2, nSims),
                    list(c("lower", "upper"), NULL))
psiEstProb <- array(NA, c(nOriginSites, nTargetSites,nSims),
                      list(originNames, targetNames,NULL))
psiCIProb <- array(NA, c(2, nOriginSites, nTargetSites,nSims),
                     list(c("lower", "upper"), originNames, targetNames,NULL))
MCestProb <- array(NA, c(nSims),
                     list(NULL))
MCCIProb <- array(NA, c(2, nSims),
                    list(c("lower", "upper"), NULL))

buffered <- c(rep(TRUE,nSims/2),rep(FALSE,nSims/2))

for(sim in 1:nSims){
  cat("Simulation", sim, "of", nSims, "at", date(), "\n")

  # simulate data with a 50km buffer (low genetic uncertainty)

  genPops <- simGeneticPops(popBoundaries = list(originSites[1, ],
                                                 originSites[2, ],
                                                 originSites[3, ]),
                            popNames = originNames,
                            res = c(50000, 50000), # 50x50km
                            bufferRegions = buffered[sim],
                            bufferDist = 50000, # 50km
                            npts = 100,
                            verbose = 0)

  data2 <- simGeneticData(genPops = genPops, psi = psiTrue,
                          originRelAbund = originRelAbund,
                          sampleSize = sampleSize[[2]],
                          originSites = originSites,
                          targetSites = targetSites,
                          captured = "target",
                          verbose = 0)

  # simple plot to see the simulated 'true' locations for both seasons #
  plot(targetSites$geometry, col = "grey50", border = "white")
  plot(originSites$geometry, add = TRUE, col = "gray", border = "white")
  plot(data2$originPointsTrue$geometry, add = TRUE, col = "black", pch = 19, cex = 0.5)
  plot(data2$targetPointsTrue$geometry, add = TRUE, col = "black", pch = 19, cex = 0.5)

  segments(x0 = st_coordinates(data2$originPointsTrue)[,1],
           y0 = st_coordinates(data2$originPointsTrue)[,2],
           x1 = st_coordinates(data2$targetPointsTrue)[,1],
           y1 = st_coordinates(data2$targetPointsTrue)[,2],
           col = rgb(150,150,150,150, m = 255))


  tp <- data2$targetPointsTrue
  or <- data2$genRaster
  ot <- data2$genProbs

  est1 <- estTransition(originSites,
                        targetSites,
                        op,
                        tp,
                        #originAssignment = ot,
                        originRaster = or, #
                        originNames = originNames,
                        targetNames = targetNames,
                        nSamples = 200,
                        isGL = isGL,
                        isTelemetry = isTelemetry,
                        isRaster = isRaster,
                        isProb = isProb,
                        captured = captured,
                        geoBias = geoBias,
                        geoVCov = geoVCov,
                        resampleProjection = sf::st_crs(targetSites),
                        nSim = 80,
                        verbose = 0,
                        dataOverlapSetting = "none")
  est1a <- estTransition(originSites,
                        targetSites,
                        op,
                        tp,
                        originAssignment = ot,
                        #originRaster = or, #
                        originNames = originNames,
                        targetNames = targetNames,
                        nSamples = 200,
                        isGL = isGL,
                        isTelemetry = isTelemetry,
                        isRaster = isProb,
                        isProb = isRaster,
                        captured = captured,
                        geoBias = geoBias,
                        geoVCov = geoVCov,
                        resampleProjection = sf::st_crs(targetSites),
                        nSim = 80,
                        verbose = 0,
                        dataOverlapSetting = "none")

  est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                      originRelAbund = originRelAbund,
                      est1)
  est2a <- estStrength(originDist = originDist, targetDist = targetDist,
                      originRelAbund = originRelAbund,
                      est1a)

  psiEstRaster[,,sim] <- est1$psi$mean
  psiCIRaster[,,,sim] <- est1$psi$simpleCI
  MCestRaster[sim] <- est2$MC$mean
  MCCIRaster[,sim] <- est2$MC$simpleCI
  psiEstProb[,,sim] <- est1a$psi$mean
  psiCIProb[,,,sim] <- est1a$psi$simpleCI
  MCestProb[sim] <- est2a$MC$mean
  MCCIProb[,sim] <- est2a$MC$simpleCI
}

png("RasterMC2a.png", width = 5, height = 3.5, units = "in", res = 600)
par(bty = "l")
plot(MCestRaster, ylim = c(0,0.60), pch = 19, las = 1, ylab = "MC",
     xlab = "Simulation", main = "Raster")
polygon(x = c(0,50,50,0), y=c(0,0,0.5,0.5),
        col = rgb(150,150,150,200, m = 255),
        border = rgb(150,150,150,200, m = 255))
segments(y0 = MCCIRaster[1,],
         y1 = MCCIRaster[2,],
         x0 = 1:100,
         x1 = 1:100)
text(20,0.55, "50km buffer")
abline(h = MCtrue, lty = 3)
dev.off()

png("ProbMC2a.png", width = 5, height = 3.5, units = "in", res = 600)
par(bty = "l")
plot(MCestProb, ylim = c(0,0.6), pch = 19, las = 1, ylab = "MC",
     xlab = "Simulation", main = "Probability Table")
polygon(x = c(0,50,50,0), y=c(0,0,0.5,0.5),
        col = rgb(150,150,150,200, m = 255),
        border = rgb(150,150,150,200, m = 255))
segments(y0 = MCCIProb[1,],
         y1 = MCCIProb[2,],
         x0 = 1:100,
         x1 = 1:100)
text(20,0.55, "50km buffer")
abline(h = MCtrue, lty = 3)
dev.off()

loc_name <- c("A","B","C")
png("psiRasterProb2a.png", width = 6, height = 6, units = "in", res = 600)
par(bty = "l")
par(mfrow = c(2,2), mar = c(4,3.5,2,2))
for(loc in 1:3){
  plot(psiTrue[loc,], ylim = c(0,1),las = 1, pch = 19,
       ylab = "psi", xlab = "target")
  text(1,0.9,loc_name[loc])
  for(i in 1:50){
    points(psiEstRaster[loc,,i]~c(0.90:4.90),pch = 19, col = "lightblue")
    points(psiEstProb[loc,,i]~c(0.80:4.80),pch = 19, col = "darkblue")
  }
  for(i in 51:100){
    points(psiEstRaster[loc,,i]~c(1.10:5.10), pch = 19, col = "lightgreen")
    points(psiEstProb[loc,,i]~c(1.20:5.20),pch = 19, col = "darkgreen")
  }
}
plot(0:2, 0:2, axes = F, col = "white", xlab = "", ylab = "")
legend(1, 1,
       c("Buffered Raster", "Buffered Prob", "No Buff Raster", "No Buff Prob"),
       pch = 19, col = c("lightblue", "darkblue", "lightgreen", "darkgreen"),
       xjust = 0.5, yjust = 0.5)
dev.off()

MAE <- function (x) mean(abs(x))
psiTrueRep <- array(psiTrue, c(nOriginSites, nTargetSites, nSims/2),
                    list(originNames, targetNames, NULL))
psiErrorRasterBuffer <- psiEstRaster[,,buffered] - psiTrueRep
psiErrorRasterNoBuffer <- psiEstRaster[,,!buffered] - psiTrueRep
psiErrorProbBuffer <- psiEstProb[,,buffered] - psiTrueRep
psiErrorProbNoBuffer <- psiEstProb[,,!buffered] - psiTrueRep
(psiBiasRasterBuffer <- apply(psiErrorRasterBuffer, 1:2, mean))
(psiBiasRasterNoBuffer <- apply(psiErrorRasterNoBuffer, 1:2, mean))
(psiBiasProbBuffer <- apply(psiErrorProbBuffer, 1:2, mean))
(psiBiasProbNoBuffer <- apply(psiErrorProbNoBuffer, 1:2, mean))
(psiMAERasterBuffer <- apply(psiErrorRasterBuffer, 1:2, MAE))
(psiMAERasterNoBuffer <- apply(psiErrorRasterNoBuffer, 1:2, MAE))
(psiMAEProbBuffer <- apply(psiErrorProbBuffer, 1:2, MAE))
(psiMAEProbNoBuffer <- apply(psiErrorProbNoBuffer, 1:2, MAE))
psiEstAll <- array(c(psiEstRaster, psiEstProb),
                   c(nOriginSites, nTargetSites, nSims/2, 2, 2),
                   list(originNames, targetNames, NULL,
                        c("Buffered", "Unbuffered"), c("Raster", "Prob")))

# Variation between simulations
(psiEstAllSD <- apply(apply(psiEstAll, c(1:2, 4:5), sd), 1:3, mean))
# Variation between probability tables and rasters
(psiEstSDAll <- apply(apply(psiEstAll, 1:4, sd), c(1:2, 4), mean))


save.image("simGeneticDataTest2a.RData")

################################################################################
#
#
# Simulate genetic data with 200km buffer - moderate uncertainty
#
#
################################################################################
for(sim in 1:nSims){
  cat("Simulation", sim, "of", nSims, "at", date(), "\n")

  # simulate data with a 50km buffer (low genetic uncertainty)

  genPops <- simGeneticPops(popBoundaries = list(originSites[1, ],
                                                 originSites[2, ],
                                                 originSites[3, ]),
                            popNames = originNames,
                            res = c(200000, 200000), # 50x50km
                            bufferRegions = buffered[sim],
                            bufferDist = 50000, # 50km
                            npts = 100,
                            verbose = 0)

  data2 <- simGeneticData(genPops = genPops, psi = psiTrue,
                          originRelAbund = originRelAbund,
                          sampleSize = sampleSize[[2]],
                          originSites = originSites,
                          targetSites = targetSites,
                          captured = "target",
                          verbose = 0)

  # simple plot to see the simulated 'true' locations for both seasons #
  plot(targetSites$geometry, col = "grey50", border = "white")
  plot(originSites$geometry, add = TRUE, col = "gray", border = "white")
  plot(data2$originPointsTrue$geometry, add = TRUE, col = "black", pch = 19, cex = 0.5)
  plot(data2$targetPointsTrue$geometry, add = TRUE, col = "black", pch = 19, cex = 0.5)

  segments(x0 = st_coordinates(data2$originPointsTrue)[,1],
           y0 = st_coordinates(data2$originPointsTrue)[,2],
           x1 = st_coordinates(data2$targetPointsTrue)[,1],
           y1 = st_coordinates(data2$targetPointsTrue)[,2],
           col = rgb(150,150,150,150, m = 255))


  tp <- data2$targetPointsTrue
  or <- data2$genRaster
  ot <- data2$genProbs

  est1 <- estTransition(originSites,
                        targetSites,
                        op,
                        tp,
                        #originAssignment = ot,
                        originRaster = or, #
                        originNames = originNames,
                        targetNames = targetNames,
                        nSamples = 200,
                        isGL = isGL,
                        isTelemetry = isTelemetry,
                        isRaster = isRaster,
                        isProb = isProb,
                        captured = captured,
                        geoBias = geoBias,
                        geoVCov = geoVCov,
                        resampleProjection = sf::st_crs(targetSites),
                        nSim = 80,
                        verbose = 0,
                        dataOverlapSetting = "none")
  est1a <- estTransition(originSites,
                         targetSites,
                         op,
                         tp,
                         originAssignment = ot,
                         #originRaster = or, #
                         originNames = originNames,
                         targetNames = targetNames,
                         nSamples = 200,
                         isGL = isGL,
                         isTelemetry = isTelemetry,
                         isRaster = isProb,
                         isProb = isRaster,
                         captured = captured,
                         geoBias = geoBias,
                         geoVCov = geoVCov,
                         resampleProjection = sf::st_crs(targetSites),
                         nSim = 80,
                         verbose = 0,
                         dataOverlapSetting = "none")

  est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                      originRelAbund = originRelAbund,
                      est1)
  est2a <- estStrength(originDist = originDist, targetDist = targetDist,
                       originRelAbund = originRelAbund,
                       est1a)

  psiEstRaster[,,sim] <- est1$psi$mean
  psiCIRaster[,,,sim] <- est1$psi$simpleCI
  MCestRaster[sim] <- est2$MC$mean
  MCCIRaster[,sim] <- est2$MC$simpleCI
  psiEstProb[,,sim] <- est1a$psi$mean
  psiCIProb[,,,sim] <- est1a$psi$simpleCI
  MCestProb[sim] <- est2a$MC$mean
  MCCIProb[,sim] <- est2a$MC$simpleCI
}


png("RasterMC2_200km.png", width = 5, height = 3.5, units = "in", res = 600)
par(bty = "l")
plot(MCestRaster, ylim = c(0,0.60), pch = 19, las = 1, ylab = "MC",
     xlab = "Simulation", main = "Raster")
polygon(x = c(0,50,50,0), y=c(0,0,0.5,0.5),
        col = rgb(150,150,150,200, m = 255),
        border = rgb(150,150,150,200, m = 255))
segments(y0 = MCCIRaster[1,],
         y1 = MCCIRaster[2,],
         x0 = 1:100,
         x1 = 1:100)
text(20,0.55, "50km buffer")
abline(h = MCtrue, lty = 3)
dev.off()

png("ProbMC2_200km.png", width = 5, height = 3.5, units = "in", res = 600)
par(bty = "l")
plot(MCestProb, ylim = c(0,0.6), pch = 19, las = 1, ylab = "MC",
     xlab = "Simulation", main = "Probability Table")
polygon(x = c(0,50,50,0), y=c(0,0,0.5,0.5),
        col = rgb(150,150,150,200, m = 255),
        border = rgb(150,150,150,200, m = 255))
segments(y0 = MCCIProb[1,],
         y1 = MCCIProb[2,],
         x0 = 1:100,
         x1 = 1:100)
text(20,0.55, "50km buffer")
abline(h = MCtrue, lty = 3)
dev.off()

loc_name <- c("A","B","C")
png("psiRasterProb2_200km.png", width = 6, height = 6, units = "in", res = 600)
par(bty = "l")
par(mfrow = c(2,2), mar = c(4,3.5,2,2))
for(loc in 1:3){
  plot(psiTrue[loc,], ylim = c(0,1),las = 1, pch = 19,
       ylab = "psi", xlab = "target")
  text(1,0.9,loc_name[loc])
  for(i in 1:50){
    points(psiEstRaster[loc,,i]~c(0.90:4.90),pch = 19, col = "lightblue")
    points(psiEstProb[loc,,i]~c(0.80:4.80),pch = 19, col = "darkblue")
  }
  for(i in 51:100){
    points(psiEstRaster[loc,,i]~c(1.10:5.10), pch = 19, col = "lightgreen")
    points(psiEstProb[loc,,i]~c(1.20:5.20),pch = 19, col = "darkgreen")
  }
}
plot(0:2, 0:2, axes = F, col = "white", xlab = "", ylab = "")
legend(1, 1,
       c("Buffered Raster", "Buffered Prob", "No Buff Raster", "No Buff Prob"),
       pch = 19, col = c("lightblue", "darkblue", "lightgreen", "darkgreen"),
       xjust = 0.5, yjust = 0.5)
dev.off()

MAE <- function (x) mean(abs(x))
psiTrueRep <- array(psiTrue, c(nOriginSites, nTargetSites, nSims/2),
                    list(originNames, targetNames, NULL))
psiErrorRasterBuffer <- psiEstRaster[,,buffered] - psiTrueRep
psiErrorRasterNoBuffer <- psiEstRaster[,,!buffered] - psiTrueRep
psiErrorProbBuffer <- psiEstProb[,,buffered] - psiTrueRep
psiErrorProbNoBuffer <- psiEstProb[,,!buffered] - psiTrueRep
(psiBiasRasterBuffer <- apply(psiErrorRasterBuffer, 1:2, mean))
(psiBiasRasterNoBuffer <- apply(psiErrorRasterNoBuffer, 1:2, mean))
(psiBiasProbBuffer <- apply(psiErrorProbBuffer, 1:2, mean))
(psiBiasProbNoBuffer <- apply(psiErrorProbNoBuffer, 1:2, mean))
(psiMAERasterBuffer <- apply(psiErrorRasterBuffer, 1:2, MAE))
(psiMAERasterNoBuffer <- apply(psiErrorRasterNoBuffer, 1:2, MAE))
(psiMAEProbBuffer <- apply(psiErrorProbBuffer, 1:2, MAE))
(psiMAEProbNoBuffer <- apply(psiErrorProbNoBuffer, 1:2, MAE))
psiEstAll <- array(c(psiEstRaster, psiEstProb),
                   c(nOriginSites, nTargetSites, nSims/2, 2, 2),
                   list(originNames, targetNames, NULL,
                        c("Buffered", "Unbuffered"), c("Raster", "Prob")))

# Variation between simulations
(psiEstAllSD <- apply(apply(psiEstAll, c(1:2, 4:5), sd), 1:3, mean))
# Variation between probability tables and rasters
(psiEstSDAll <- apply(apply(psiEstAll, 1:4, sd), c(1:2, 4), mean))


save.image("simGeneticDataTest2_200km.RData")
