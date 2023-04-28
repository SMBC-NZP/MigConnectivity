devtools::install_github("SMBC-NZP/MigConnectivity", ref = "devpsi2")
library(MigConnectivity)

run.number <- 1
originNames <- c("A", "B", "C")
nOriginSites <- length(originNames)
targetNames <- as.character(1:4)
nTargetSites <- length(targetNames)

psiTrue <- matrix(c(0.55, 0.25, 0.15, 0.05,
                    0.08, 0.4, 0.4, 0.12,
                    0.05, 0.1, 0.2, 0.65), nOriginSites, nTargetSites,
                  TRUE, list(originNames, targetNames))
psiTrue
rowSums(psiTrue)
banded <- c(30, 60, 20); names(banded) <- originNames
rTrue <- c(0.1, 0.1, 0.1, 0.1) ; names(rTrue) <- targetNames
nSims <- 20
reencountered <- psiCalc <- psiEstMCMC <- psiEstBoot <- vector("list", nSims)

set.seed(run.number)
for (i in 1:nSims) {
  cat("Simulation", i, "of", nSims, "at", date(), "\n")
  dataCMR <- simCMRData(psiTrue, banded, rTrue)
  #dataTelemetry <- simTelemetryData(psiTrue, banded/10)
  reencountered[[i]] <- dataCMR$reencountered
  print(reencountered[[i]])
  psiCalc[[i]] <- calcTransition(banded = banded, reencountered = reencountered[[i]],
                                 originNames = originNames, targetNames = targetNames,
                                 method = "BFGS")
  print(psiCalc[[i]])
  psiEstMCMC[[i]] <- estTransition(originNames = originNames, targetNames = targetNames,
                              nSamples = 36000, banded = banded,
                              reencountered = reencountered[[i]],
                              method = "MCMC", verbose = 0, nBurnin = 104000)
  print(psiEstMCMC[[i]])
  psiPrior <- (psiTrue * 10)^6
  print(psiPrior)
  psiEstMCMC2 <- estTransition(originNames = originNames, targetNames = targetNames,
                                   nSamples = 36000, banded = banded,
                                   reencountered = reencountered[[i]],
                                   method = "MCMC", verbose = 1, nBurnin = 104000,
                                   psiPrior = psiPrior) #matrix(c(5000,1,1,1,
                                   #                     1,5000,5000,1,
                                   #                     1,1,2,5000), nOriginSites,
                                   #                   nTargetSites, byrow = T))
  psiEstMCMC2
  psiEstMCMC2$psi$mean - psiEstMCMC[[i]]$psi$mean
  psiEstBoot[[i]] <- estTransition(originNames = originNames, targetNames = targetNames,
                              nSamples = 36000, banded = banded,
                              reencountered = reencountered[[i]],
                              method = "bootstrap", verbose = 0)
  save.image(paste0("~/MC work/simConnectivity/simCMRresultsA",
                    run.number, ".RData"))
}
psiErrorBoot <- sapply(psiEstBoot, function(x) x$psi$mean - psiTrue, simplify = "array")
psiErrorMCMC <- sapply(psiEstMCMC, function(x) x$psi$mean - psiTrue, simplify = "array")
psiRhat <- sapply(psiEstMCMC, function(x) max(x$BUGSoutput$summary[,"Rhat"]))
psiConvergence <- psiRhat < 1.1
psiErrorMCMC
(psiBiasMCMC <- apply(psiErrorMCMC, 1:2, mean))
(psiBiasBoot <- apply(psiErrorBoot, 1:2, mean))
(psiMAEMCMC <- apply(psiErrorMCMC, 1:2, function(x) mean(abs(x), na.rm = T)))
(psiMAEBoot <- apply(psiErrorBoot, 1:2, function(x) mean(abs(x), na.rm = T)))
library(coda)
psiListsMCMC <- lapply(psiEstMCMC, function(x) as.mcmc.list((x$BUGSoutput)))
for (i in 1:nSims) {
  if (!psiConvergence[i])
    plot(psiListsMCMC[[i]])
}
