devtools::install_github("SMBC-NZP/MigConnectivity", ref = "devpsi2")
library(MigConnectivity)

run.number <- 8
originNames <- c("A", "B", "C")
nOriginSites <- length(originNames)
targetNames <- as.character(1:4)
nTargetSites <- length(targetNames)

psiTrue <- matrix(c(0.5, 0.25, 0.15, 0.1,
                    0.15, 0.4, 0.25, 0.2,
                    0.1, 0.15, 0.2, 0.55), nOriginSites, nTargetSites,
                  TRUE, list(originNames, targetNames))
psiTrue
rowSums(psiTrue)
banded <- c(3000, 6000, 12000); names(banded) <- originNames
rTrue <- c(0.5, 0.6, 0.4, 0.7) ; names(rTrue) <- targetNames
nSims <- 12
reencountered <- psiCalc <- psiEstMCMC <- psiEstBoot <- vector("list", nSims)

set.seed(run.number)
for (i in 1:nSims) {
  cat("Simulation", i, "of", nSims, "at", date(), "\n")
  dataCMR <- simCMRData(psiTrue, banded, rTrue)
  #dataTelemetry <- simTelemetryData(psiTrue, banded/10)
  reencountered[[i]] <- dataCMR$reencountered
  psiCalc[[i]] <- calcTransition(banded = banded, reencountered = reencountered[[i]],
                                 originNames = originNames, targetNames = targetNames,
                                 method = "BFGS")

  psiEstMCMC[[i]] <- estTransition(originNames = originNames, targetNames = targetNames,
                              nSamples = 36000, banded = banded,
                              reencountered = reencountered[[i]],
                              method = "MCMC", verbose = 0, nBurnin = 104000)

  psiEstBoot[[i]] <- estTransition(originNames = originNames, targetNames = targetNames,
                              nSamples = 36000, banded = banded,
                              reencountered = reencountered[[i]],
                              method = "bootstrap", verbose = 0)
  save.image(paste0("~/MC work/simConnectivity/simCMRresults",
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
