devtools::install_github("SMBC-NZP/MigConnectivity", ref = "devpsi2")
library(MigConnectivity)

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
banded <- c(4000, 8000, 16000); names(banded) <- originNames
rTrue <- c(0.5, 0.6, 0.4, 0.7) ; names(rTrue) <- targetNames
nSims <- 20

set.seed(9001)
for (i in 1:nSims) {
dataCMR <- simCMRData(psiTrue, banded, rTrue)
calcTransition(banded = banded, reencountered = dataCMR$reencountered,
               originNames = originNames, targetNames = targetNames,
               method = "BFGS")
# $psi
#           1         2         3         4
# A 0.5263609 0.2490784 0.1182439 0.1063167
# B 0.1584947 0.4032901 0.2281568 0.2100584
# C 0.1126774 0.1564856 0.1721251 0.5587118
#
# $r
#         1         2         3         4
# 0.4692599 0.5771284 0.4651415 0.6795756


psiEstMCMC <- estTransition(originNames = originNames, targetNames = targetNames,
                            nSamples = 24000, banded = banded,
                            reencountered = dataCMR$reencountered,
                            method = "MCMC", verbose = 1)
psiEstMCMC$psi$mean
#           [,1]      [,2]      [,3]      [,4]
# [1,] 0.5301817 0.2097081 0.1602976 0.0998126
# [2,] 0.1602001 0.3380999 0.3051014 0.1965987
# [3,] 0.1140086 0.1315385 0.2320618 0.5223911
psiEstMCMC$r$mean

psiEstBoot <- estTransition(originNames = originNames, targetNames = targetNames,
                            nSamples = 24000, banded = banded,
                            reencountered = dataCMR$reencountered,
                            method = "bootstrap", verbose = 1)
}
psiEstBoot$psi$mean
#           1         2         3          4
# A 0.5030111 0.2107851 0.1994911 0.08671263
# B 0.1813598 0.3878865 0.2550761 0.17567751
# C 0.1395595 0.1789259 0.2343083 0.44720638
psiEstBoot$r$mean
#         1         2         3         4
# 0.4966539 0.5129865 0.5481648 0.6927423

psiEstMCMC$psi$mean - psiTrue
#            1           2          3             4
# A 0.03018172 -0.04029191 0.01029759 -0.0001874034
# B 0.01020005 -0.06190007 0.05510136 -0.0034013361
# C 0.01400862 -0.01846149 0.03206181 -0.0276089420

psiEstBoot$psi$mean - psiTrue
#             1           2           3           4
# A 0.003011094 -0.03921485 0.049491127 -0.01328737
# B 0.031359841 -0.01211348 0.005076122 -0.02432249
# C 0.039559502  0.02892585 0.034308265 -0.10279362

