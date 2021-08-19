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
nTargetSites <- 6
psiTrue <- Combined$psi$mean
psiTrue[2, 3] <- psiTrue[2, 1]
psiTrue <- psiTrue[, -1]
psiTrue <- round(psiTrue, 1)
rowSums(psiTrue)
originNames <- LETTERS[1:nOriginSites]
targetNames <- as.character(1:nTargetSites)
dimnames(psiTrue) <- list(originNames, targetNames)
geoBias <- OVENdata$geo.bias
geoVCov <- OVENdata$geo.vcov
# Changed these from PABU estimates so we'd average at least one GL per
# population in first and third scenario
originRelAbund <- c(1/30, 2/30, 27/30)
originSites <- OVENdata$originSites
targetSites <- OVENdata$targetSites[2:3, ]
originSitesPABU

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

sf::st_crs(originSitesPABU)
genPops <- simGeneticPops(popBoundaries = list(originSitesPABU[1, ],
                                               originSitesPABU[2, ],
                                               originSitesPABU[3, ]),
                          popNames = originNamesPABU, res = c(10000, 10000),
                          verbose = 1)

S <- vector("list", nScenarios)
S[[1]] <- S[[2]] <- S[[3]] <- S[[4]] <- S[[5]] <- S[[6]] <- S[[7]] <- 1
p <- vector("list", nScenarios)
for (i in 1:nScenarios)
  p[[i]] <- list(1, 1)

sampleSize <- vector("list", nScenarios)
sampleSize[[1]] <- list(30, NULL)
sampleSize[[2]] <- list(NULL, 300)
sampleSize[[3]] <- list(30, 300)
sampleSize[[4]] <- list(c(10, 10, 10), NULL)
sampleSize[[5]] <- list(NULL, rep(50, nTargetSites))
sampleSize[[6]] <- list(c(30, 0, 0), c(75, 75, 75, 75, 0, 0))
sampleSize[[7]] <- list(c(5, 0, 0), c(75, 75, 75, 75, 0, 0))
sampleSizeGL <- lapply(sampleSize, function(x) x[[2]] <- NULL)
sampleSizeGeno <- lapply(sampleSize, function(x) x[[2]])

captured <- vector("list", nScenarios)
for (i in 1:nScenarios) {
  if (!is.null(sampleSize[[i]][[1]]))
    captured[[i]] <- rep("origin", sum(sampleSize[[i]][[1]]))
  if (!is.null(sampleSize[[i]][[2]]))
    captured[[i]] <- c(captured[[i]], rep("target", sum(sampleSize[[i]][[2]])))
}

originDist <- OVENdata$originDist
targetDist <- OVENdata$targetDist[-1,-1]
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
sampleSizes <- array(NA, c(nOriginSites, nTargetSites, nScenarios, nSims),
                list(originNames, targetNames, scenarios, NULL))
dataStore <- vector("list", nSims)
for (sim in 1:nSims) {
  cat("Simulation", sim, "of", nSims, "at", date(), " ")
  dataStore[[sim]] <- vector("list", nScenarios)
  for (sc in 1:nScenarios) {
    if (!is.null(sampleSizeGL[[sc]][[1]]))
      data1 <- simGL(psi = psiTrue, originRelAbund = originRelAbund,
                     sampleSize = sampleSizeGL[[sc]],
                     originSites = originSites, targetSites = targetSites,
                     captured = captured[[sc]],
                     geoBias, geoVCov, geoBiasOrigin, geoVCovOrigin,
                     S = S[[sc]], p = p[[sc]],
                     requireEveryOrigin = is.null(sampleSizeGeno[[sc]]))
    else
      data1 <- NULL
    if (!is.null(sampleSizeGeno[[sc]]))
      data2 <- simGeneticData(genPops = genPops, psi = psiTrue,
                              originRelAbund = originRelAbund,
                              sampleSize = sampleSizeGeno[[sc]],
                              originSites = originSites,
                              targetSites = targetSites,
                              captured = captured[[sc]])
    else
      data2 <- NULL

    # test1 <- MigConnectivity:::locSample(rep(T, sum(data1$recaptured)),
    #                                      rep(F, sum(data1$recaptured)),
    #                                      rep(F, sum(data1$recaptured)),
    #                                      rep(F, sum(data1$recaptured)),
    #                                      geoBias = geoBias,
    #                                      geoVCov = geoVCov,
    #                                      points = data1$targetPointsObs,
    #                                      sites = targetSites,
    #                                      resampleProjection = sf::st_crs(targetSites),
    #                                      nSim = 1000, maxTries = 300)
    est1 <- estTransition(originSites, targetSites,
                          data1$originPointsObs, data1$targetPointsObs,
                          originNames = originNames, targetNames = targetNames,
                          nSamples = 200, isGL = TRUE,
                          captured = captured[[sc]][which(data1$recaptured==1)],
                          geoBias = geoBias, geoVCov = geoVCov,
                          geoBiasOrigin = geoBiasOrigin,
                          geoVCovOrigin = geoVCovOrigin,
                          resampleProjection = sf::st_crs(targetSites),
                          nSim = 400, verbose = 0)
    est2 <- estStrength(originDist = originDist, targetDist = targetDist,
                        originRelAbund = originRelAbund,
                        est1)
    est3 <- estMantel(data1$targetPointsObs, data1$originPointsObs, T,
                      geoBias, geoVCov, targetSites, 200, 300, 0,
                      resampleProjection = sf::st_crs(targetSites),
                      geoBiasOrigin = geoBiasOrigin,
                      geoVCovOrigin = geoVCovOrigin, originSites = originSites,
                      captured = captured[[sc]][which(data1$recaptured==1)])
    psiEst[,,sc,sim] <- est1$psi$mean
    psiCI[,,,sc,sim] <- est1$psi$simpleCI
    MCest[sc,sim] <- est2$MC$mean
    MCCI[,sc,sim] <- est2$MC$simpleCI
    MantelEst[sc,sim] <- est3$corr$mean
    MantelCI[,sc,sim] <- est3$corr$simpleCI
    sampleSizes[,,sc,sim] <- table(data1$originAssignment[which(data1$recaptured==1)],
                                   data1$targetAssignment[which(data1$recaptured==1)])
    dataStore[[sim]][[sc]] <- data1
    cat(sc)
    save(psiEst, psiCI, MCest, MCCI, MantelEst, MantelCI, sampleSizes, sim, sc,
         dataStore, file = 'testGL2b.RData')
  }
  cat("\n")
}

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
