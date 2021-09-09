################################################################################
# Simulation to test new estTransition/estMantel/estStrength options, including
# data with location uncertainty on origin side
# This version simulates GL data released on origin and genoscape data gathered
# on the target side
################################################################################
library(raster)
library(sf)
library(ks)
library(MigConnectivity)
# library(spdep)
# library(mapview)

nSims <- 1000
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

# Transition probability (psi) estimates (mean):
#            Pacific and Interior Mexico Atlantic Lowland Mexico Central America
# Central                        0.10529                  0.1068       0.7878945
# Louisiana                      0.05526                  0.4737       0.4710340
# East Coast                     0.00000                  0.0000       0.0001181
#            Southeastern US Caribbean
# Central             0.0000    0.0000
# Louisiana           0.0000    0.0000
# East Coast          0.5433    0.4566

psiTrue <- array(0, c(nOriginSites, nTargetSites))
psiTrue[1,] <- c(0.10529, 0.10682, 0.78789, 0, 0)
psiTrue[2,] <- c(0.05526, 0.47371, 0.47103, 0, 0)
psiTrue[3,] <- c(0, 0, 0.00012, 0.54331, 0.45657)
rowSums(psiTrue)
# psiTrue <- psiPABU$psi$mean

originNames <- LETTERS[1:nOriginSites]
targetNames <- as.character(1:nTargetSites)
dimnames(psiTrue) <- list(originNames, targetNames)
geoBias <- c(22250.92, -32407.14)
geoVCov <- matrix(c(5637340653, -2692682155, -2692682155, 6012336962),
                  nrow = 2, ncol = 2)
# Changed these from PABU estimates so we'd average at least one GL per
# population in first and third scenario

originRelAbund <- c(27/30, 2/30, 1/30)

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

originSites <- sf::st_read("data-raw/PABU_breeding_regions.shp")
targetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")
originSites$originSite <- originNames
targetSites$targetSite <- targetNames

# sf::st_crs(originSitesPABU)

# looks like originSites needs to be projected
originSites <- sf::st_transform(originSites, "ESRI:102010")
targetSites <- sf::st_transform(targetSites, "ESRI:102010")

genPops <- simGeneticPops(popBoundaries = list(originSites[1, ],
                                               originSites[2, ],
                                               originSites[3, ]),
                          popNames = originNames, res = c(50000, 50000),
                          bufferRegions = F, bufferDist = 200000,
                          npts = 1000,
                          verbose = 1)
calibSampleSize <- table(factor(unlist(overlapT2)[-breeders],
                                levels = 1:nTargetSites))
data2 <- simGeneticData(genPops = genPops, psi = psiTrue,
                        originRelAbund = originRelAbund,
                        sampleSize = calibSampleSize,
                        originSites = originSites,
                        targetSites = targetSites,
                        captured = "target",
                        verbose = 1)
mean(apply(data2$genProbs, 1, var))
mean(apply(originAssignmentPABU[-breeders, ], 1, var))
mean(apply(data2$genProbs, 1, max)>0.99999)
mean(apply(originAssignmentPABU[-breeders, ], 1, max)>0.99999)
mean(apply(data2$genProbs, 1, max)>0.99)
mean(apply(originAssignmentPABU[-breeders, ], 1, max)>0.99)
S <- vector("list", nScenarios)
S[[1]] <- S[[2]] <- S[[3]] <- S[[4]] <- S[[5]] <- S[[6]] <- S[[7]] <- 1
p <- vector("list", nScenarios)
for (i in 1:nScenarios)
  p[[i]] <- list(1, 1)

sampleSize <- vector("list", nScenarios)
sampleSize[[1]] <- list(30 * originRelAbund, NULL)
sampleSize[[2]] <- list(NULL, round(rev$targetRelAbund * 300))
sampleSize[[3]] <- list(30 * originRelAbund, round(rev$targetRelAbund * 300))
sampleSize[[4]] <- list(c(10, 10, 10), NULL)
sampleSize[[5]] <- list(NULL, rep(60, nTargetSites))
sampleSize[[6]] <- list(c(0, 0, 30), c(100, 100, 100, 0, 0))
sampleSize[[7]] <- list(c(0, 0, 5), c(100, 100, 100, 0, 0))
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
    # isRaster[[i]] <- c(isRaster[[i]], rep(F, sum(sampleSize[[i]][[2]])))
    isRaster[[i]] <- c(isRaster[[i]], rep(T, sum(sampleSize[[i]][[2]])))
    isTelemetry[[i]] <- c(isTelemetry[[i]], rep(FALSE, sum(sampleSize[[i]][[2]])))
    #isProb[[i]] <- c(isProb[[i]], rep(T, sum(sampleSize[[i]][[2]])))
    isProb[[i]] <- c(isProb[[i]], rep(F, sum(sampleSize[[i]][[2]])))
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
MantelEst <- array(NA, c(nScenarios, nSims),
               list(scenarios, NULL))
MantelCI <- array(NA, c(2, nScenarios, nSims),
              list(c("lower", "upper"), scenarios, NULL))
# sampleSizes <- array(NA, c(nOriginSites, nTargetSites, nScenarios, nSims),
#                 list(originNames, targetNames, scenarios, NULL))

dataStore <- vector("list", nSims)


for (sim in 1:nSims) {
  cat("Simulation", sim, "of", nSims, "at", date(), " ")
  dataStore[[sim]] <- vector("list", nScenarios)
  for (sc in 1:nScenarios) {
    or <- NULL
    if (!is.null(sampleSizeGL[[sc]][[1]])){
      data1 <- simGL(psi = psiTrue, originRelAbund = originRelAbund,
                     sampleSize = sampleSizeGL[[sc]],
                     originSites = originSites, targetSites = targetSites,
                     captured = "origin",
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
                              verbose = 1)
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
                          #originAssignment = ot,
                          originRaster = or, #
                          originNames = originNames,
                          targetNames = targetNames,
                          nSamples = 1000, isGL = isGL[[sc]],
                          isTelemetry = isTelemetry[[sc]],
                          isRaster = isRaster[[sc]],
                          isProb = isProb[[sc]],
                          captured = captured[[sc]],
                          geoBias = geoBias, geoVCov = geoVCov,
                          resampleProjection = sf::st_crs(targetSites),
                          nSim = 80, verbose = 3,
                          dataOverlapSetting = "none")
    est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                        originRelAbund = originRelAbund,
                        est1)
    est3 <- estMantel(tp, op, isGL[[sc]],
                      geoBias, geoVCov, targetSites, 10, 80, 3,
                      resampleProjection = sf::st_crs(targetSites),
                      originSites = originSites,
                      captured = captured[[sc]],
                      isTelemetry = isTelemetry[[sc]],
                      isRaster = isRaster[[sc]],
                      originRaster = or,
                      dataOverlapSetting = "none")
    psiEst[,,sc,sim] <- est1$psi$mean
    psiCI[,,,sc,sim] <- est1$psi$simpleCI
    MCest[sc,sim] <- est2$MC$mean
    MCCI[,sc,sim] <- est2$MC$simpleCI
    MantelEst[sc,sim] <- est3$corr$mean
    MantelCI[,sc,sim] <- est3$corr$simpleCI
    # sampleSizes[,,sc,sim] <- table(data1$originAssignment[which(data1$recaptured==1)],
    #                                data1$targetAssignment[which(data1$recaptured==1)])
    dataStore[[sim]][[sc]] <- list(data1 = data1, data2 = data2)
    cat(sc)
    save(psiEst, psiCI, MCest, MCCI, MantelEst, MantelCI, sim, sc, #sampleSizes,
         dataStore, file = 'testConnectivity3a.RData')
  }
  cat("\n")
}

# nSims <- 20
# psiEstRaster <- array(NA, c(nOriginSites, nTargetSites, nSims),
#                       list(originNames, targetNames, NULL))
# psiSDRaster <- array(NA, c(nOriginSites, nTargetSites, nSims),
#                       list(originNames, targetNames, NULL))
# psiCIRaster <- array(NA, c(2, nOriginSites, nTargetSites, nSims),
#                      list(c("lower", "upper"), originNames, targetNames, NULL))
# psiEstProb <- array(NA, c(nOriginSites, nTargetSites, nSims),
#                     list(originNames, targetNames, NULL))
# psiSDProb <- array(NA, c(nOriginSites, nTargetSites, nSims),
#                     list(originNames, targetNames, NULL))
# psiCIProb <- array(NA, c(2, nOriginSites, nTargetSites, nSims),
#                    list(c("lower", "upper"), originNames, targetNames, NULL))
# for (sim in 1:nSims) {
#   cat("Simulation", sim, "of", nSims, "at", date(), " ")
#   for (sc in 5) {
#     or <- NULL
#     if (!is.null(sampleSizeGL[[sc]][[1]])){
#       data1 <- simGL(psi = psiTrue, originRelAbund = originRelAbund,
#                      sampleSize = sampleSizeGL[[sc]],
#                      originSites = originSites, targetSites = targetSites,
#                      captured = "origin",
#                      geoBias, geoVCov,
#                      S = S[[sc]], p = p[[sc]],
#                      requireEveryOrigin = is.null(sampleSizeGeno[[sc]]),
#                      verbose = 0)
#       op <- data1$originPointsObs
#       tp <- data1$targetPointsObs
#     }else{
#       data1 <- NULL
#       op <- NULL
#       tp <- NULL
#     }
#     if (!is.null(sampleSizeGeno[[sc]])){
#       data2 <- simGeneticData(genPops = genPops, psi = psiTrue,
#                               originRelAbund = originRelAbund,
#                               sampleSize = sampleSizeGeno[[sc]],
#                               originSites = originSites,
#                               targetSites = targetSites,
#                               captured = "target",
#                               verbose = 0)
#       tp <- rbind(tp, data2$targetPointsTrue)
#       or <- data2$genRaster
#       ot <- data2$genProbs
#       #originSites <- sf::st_transform(originSites, crs(or, TRUE))
#       #crs(or) <- sf::st_crs(originSites)
#     }else{
#       data2 <- NULL
#       or <- NULL
#       ot <- NULL
#     }
#     est1 <- estTransition(originSites,
#                           targetSites,
#                           op,
#                           tp,
#                           #originAssignment = ot,
#                           originRaster = or, #
#                           originNames = originNames,
#                           targetNames = targetNames,
#                           nSamples = 1000, isGL = isGL[[sc]],
#                           isTelemetry = isTelemetry[[sc]],
#                           isRaster = isRaster[[sc]],
#                           isProb = isProb[[sc]],
#                           captured = captured[[sc]],
#                           geoBias = geoBias, geoVCov = geoVCov,
#                           resampleProjection = sf::st_crs(targetSites),
#                           nSim = 80, verbose = 1,
#                           dataOverlapSetting = "none")
#     est1a <- estTransition(originSites,
#                            targetSites,
#                            op,
#                            tp,
#                            originAssignment = ot,
#                            #originRaster = or, #
#                            originNames = originNames,
#                            targetNames = targetNames,
#                            nSamples = 1000, isGL = isGL[[sc]],
#                            isTelemetry = isTelemetry[[sc]],
#                            isRaster = isProb[[sc]],
#                            isProb = isRaster[[sc]],
#                            captured = captured[[sc]],
#                            geoBias = geoBias, geoVCov = geoVCov,
#                            resampleProjection = sf::st_crs(targetSites),
#                            nSim = 80, verbose = 0,
#                            dataOverlapSetting = "none")
#     psiEstRaster[,,sim] <- est1$psi$mean
#     psiCIRaster[,,,sim] <- est1$psi$simpleCI
#     psiSDRaster[,,sim] <- est1$psi$se
#     psiEstProb[,,sim] <- est1a$psi$mean
#     psiSDProb[,,sim] <- est1a$psi$se
#     psiCIProb[,,,sim] <- est1a$psi$simpleCI
#     dataStore[[sim]] <- list(data1 = data1, data2 = data2)
#     save(psiEstRaster, psiCIRaster, psiSDRaster, psiEstProb, psiCIProb, psiSDProb,#sampleSizes,
#          dataStore, file = 'testRasterProb3a.RData')
#   }
#   cat("\n")
# }
#
# nSimsDone <- sim - 1
# psiEstRaster <- psiEstRaster[,,1:nSimsDone]
# psiCIRaster <- psiCIRaster[,,,1:nSimsDone]
# psiSDRaster <- psiSDRaster[,,1:nSimsDone]
# psiEstProb <- psiEstProb[,,1:nSimsDone]
# psiCIProb <- psiCIProb[,,,1:nSimsDone]
# psiSDProb <- psiSDProb[,,1:nSimsDone]
# psiEstAll <- array(c(psiEstRaster, psiEstProb),
#                    c(nOriginSites, nTargetSites, nSimsDone, 2),
#                    list(originNames, targetNames, NULL, c("Raster", "Prob")))
#
# psiEstRasterSD <- apply(psiEstRaster, 1:2, sd)
# psiEstProbSD <- apply(psiEstProb, 1:2, sd)
# psiEstAllSD <- (psiEstRasterSD + psiEstProbSD) / 2
# psiEstSDAll <- apply(apply(psiEstAll, c(1,2,3), sd), 1:2, mean)
# psiEstAllSD
# # Average sd between simulations
# #              1           2          3           4           5
# # A 0.0007509665 0.000000000 0.01306195 0.008530197 0.018982294
# # B 0.0454965864 0.005774737 0.02549208 0.038386530 0.003482010
# # C 0.0304489928 0.008777679 0.01751065 0.025716835 0.004284463
# psiEstSDAll
# # Average sd between data types
# #              1            2            3            4            5
# # A 0.0000433691 0.0000000000 0.0005234190 0.0002079959 3.692237e-04
# # B 0.0018138104 0.0001174937 0.0012673675 0.0013989689 1.313988e-04
# # C 0.0009039616 0.0010716161 0.0007846056 0.0007073056 7.510227e-05
# # More variation between simulations than between data types (good!)

nSimsDone <- sim
psiEst <- psiEst[,,,1:nSimsDone]
psiCI <- psiCI[,,,,1:nSimsDone]
psiTrue2 <- array(c(psiTrue), c(nOriginSites, nTargetSites, nScenarios, nSimsDone),
                  list(originNames, targetNames, scenarios, NULL))
all(psiTrue2[,,2,2]==psiTrue)
psiEstError <- psiEst - psiTrue2
(psiEstBias <- apply(psiEstError, 1:3, mean))
(psiEstMAE <- apply(psiEstError, 1:3, function(x) mean(abs(x))))
(psiSimpleCover <- apply(psiCI[1,,,,] <= psiTrue2 & psiCI[2,,,,] >= psiTrue2,
                        1:3, mean))
(psiEstMAE2 <- apply(psiEstMAE, 3, mean))

MCest <- MCest[,1:nSimsDone]
MCCI <- MCCI[,,1:nSimsDone]
MCestError <- MCest - MCtrue
(MCestBias <- apply(MCestError, 1, mean))
(MCestMAE <- apply(MCestError, 1, function(x) mean(abs(x))))
(MCsimpleCover <- apply(MCCI[1,,] <= MCtrue & MCCI[2,,] >= MCtrue,
                         1, mean))

library(ggplot2)
psiEst.df <- data.frame(origin = factor(rep(originNames, nTargetSites * nScenarios * nSimsDone),
                                        levels = originNames),
                        target = factor(rep(targetNames, nScenarios*nSimsDone, each = nOriginSites),
                                        levels = targetNames),
                        scenario = factor(rep(scenarios, nSimsDone, each=nOriginSites*nTargetSites),
                                          levels = scenarios),
                        psi = c(psiEst))
summary(psiEst.df)
psiTrue.df <- data.frame(origin = factor(rep(originNames, nTargetSites),
                                         levels = originNames),
                         target = factor(rep(targetNames, each = nOriginSites),
                                         levels = targetNames),
                         psi = c(psiTrue))
bias.df <- data.frame(origin = factor(rep(originNames, nTargetSites * nScenarios),
                                      levels = originNames),
                      target = factor(rep(targetNames, nScenarios, each = nOriginSites),
                                      levels = targetNames),
                      scenario = factor(rep(scenarios, each=nOriginSites*nTargetSites),
                                        levels = scenarios),
                      bias = paste("Bias =", format(psiEstBias, digits = 1)),
                      x = 0.5, y = 50)
mae.df <- data.frame(origin = factor(rep(originNames, nTargetSites * nScenarios),
                                      levels = originNames),
                      target = factor(rep(targetNames, nScenarios, each = nOriginSites),
                                      levels = targetNames),
                      scenario = factor(rep(scenarios, each=nOriginSites*nTargetSites),
                                        levels = scenarios),
                      x = 0.5, y = 42,
                      mae = paste("MAE =", format(psiEstMAE, digits = 2)))
cov1.df <- data.frame(origin = factor(rep(originNames, nTargetSites * nScenarios),
                                      levels = originNames),
                      target = factor(rep(targetNames, nScenarios, each = nOriginSites),
                                      levels = targetNames),
                      scenario = factor(rep(scenarios, each=nOriginSites*nTargetSites),
                                        levels = scenarios),
                      x = 0.5, y = 34,
                      coverSimple = paste("C (q) =", format(psiSimpleCover, digits = 2)),
                      cover = paste("C =", format(psiSimpleCover, digits = 2)))
cov2.df <- data.frame(origin = factor(rep(originNames, nTargetSites * nScenarios),
                                      levels = originNames),
                      target = factor(rep(targetNames, nScenarios, each = nOriginSites),
                                      levels = targetNames),
                      scenario = factor(rep(scenarios, each=nOriginSites*nTargetSites),
                                        levels = scenarios),
                      x = 0.5, y = 34,
                      coverBC = paste("C (bc) =", format(psibcCover, digits = 2)),
                      cover = paste("C =", format(psibcCover, digits = 2)))
plot1 <- ggplot(subset(psiEst.df), aes(psi)) +#, scenario == "released origin"
  geom_histogram() + theme_bw() + facet_grid(origin + scenario ~ target) +
  geom_vline(aes(xintercept = psi), psiTrue.df) +
  geom_text(aes(x = x, y = y, label = bias), data = bias.df, size = 3) +
  geom_text(aes(x = x, y = y, label = mae), data = mae.df, size = 3) +
#  geom_text(aes(x = x, y = y, label = coverSimple), data = cov1.df, size = 3) +
  geom_text(aes(x = x, y = y, label = cover), data = cov2.df, size = 3)
plot1
png('simGL_psi_hists2.png', height = 8, width = 6, units = "in", res = 1200)
plot1
dev.off()
MCest.df <- data.frame(scenario = factor(rep(scenarios, nSimsDone),
                                         levels = scenarios),
                       MC = c(MCest))
MCbias.df <- data.frame(scenario = factor(scenarios, levels = scenarios),
                        bias = paste("Bias =", format(MCestBias, digits = 1)),
                        x = 0.5, y = 30)
MCmae.df <- data.frame(scenario = factor(scenarios, levels = scenarios),
                       x = 0.5, y = 28,
                       mae = paste("MAE =", format(MCestMAE, digits = 2)))
MCcov1.df <- data.frame(scenario = factor(scenarios, levels = scenarios),
                        x = 0.5, y = 26,
                        coverSimple = paste("C (q) =", format(MCsimpleCover, digits = 2)),
                        cover = paste("C =", format(MCsimpleCover, digits = 2)))
MCcov2.df <- data.frame(scenario = factor(scenarios, levels = scenarios),
                        x = 0.5, y = 26,
                        coverSimple = paste("C (q) =", format(MCbcCover, digits = 2)),
                        cover = paste("C =", format(MCbcCover, digits = 2)))
plot2 <- ggplot(MCest.df, aes(MC)) +#, scenario == "released origin"
  geom_histogram() + theme_bw() + facet_grid(. ~ scenario) +
  geom_vline(xintercept = MCtrue) +
  geom_text(aes(x = x, y = y, label = bias), data = MCbias.df, size = 4) +
  geom_text(aes(x = x, y = y, label = mae), data = MCmae.df, size = 4) +
  #  geom_text(aes(x = x, y = y, label = coverSimple), data = cov1.df, size = 3) +
  geom_text(aes(x = x, y = y, label = cover), data = MCcov2.df, size = 4)
plot2
png('simGL_MC_hists2.png', height = 4, width = 6, units = "in", res = 1200)
plot2
dev.off()

oa <- factor(data2$originAssignment, labels = originNames)
ggplot() +
  geom_sf(data = originSites, aes(fill = originSite)) +
  geom_sf(data = data2$originPointsTrue, size = 2,
          aes(shape = oa)) +
  ggtitle("Breeding Sites") +
  coord_sf()
ta <- factor(data2$targetAssignment, labels = targetNames)
ggplot() +
  geom_sf(data = targetSites, aes(fill = targetSite)) +
  geom_sf(data = data2$targetPointsTrue, size = 2,
          aes(shape = ta)) +
  ggtitle("Wintering Sites") +
  coord_sf()
oa1 <- factor(data1$originAssignment, labels = originNames)
ggplot() +
  geom_sf(data = originSites, aes(fill = originSite)) +
  geom_sf(data = data1$originPointsTrue, size = 2,
          aes(shape = oa1)) +
  ggtitle("Breeding Sites") +
  coord_sf()
test <- MigConnectivity:::randomPoints(sites = originSites,
                                       assignment = data1$originAssignment)
ggplot() +
  geom_sf(data = originSites, aes(fill = originSite)) +
  geom_sf(data = test, size = 2,
          aes(shape = oa1)) +
  ggtitle("Breeding Sites") +
  coord_sf()
