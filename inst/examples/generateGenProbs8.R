# TROUBLESHOOTING SIMULATING GENETIC DATA

################################################################################
library(raster)
library(sf)
library(ks)
library(MigConnectivity)
library(VGAM)
# library(mapview)

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

originNames <- LETTERS[1:nOriginSites]
targetNames <- as.character(1:nTargetSites)
dimnames(psiTrue) <- list(originNames, targetNames)
geoBias <- c(22250.92, -32407.14)
geoVCov <- matrix(c(5637340653, -2692682155, -2692682155, 6012336962),
                  nrow = 2, ncol = 2)

originRelAbund <- c(0.927, 0.05055, 0.02245)
#originRelAbund <- rep(1/3, 3)

rev <- MigConnectivity:::reversePsiRelAbund(psiTrue, originRelAbund)

originSites <- sf::st_read("data-raw/PABU_breeding_regions.shp")
targetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")

originSites$originSite <- originNames
targetSites$targetSite <- targetNames

# sf::st_crs(originSitesPABU)

# looks like originSites needs to be projected
originSites <- sf::st_transform(originSites, "ESRI:102010")
targetSites <- sf::st_transform(targetSites, "ESRI:102010")

set.seed(1020)
shapes <- matrix(c(5.5, 0.25, 0.01,
                   0.25, 5.5, 0.01,
                   0.004, 0.006, 5.75), 3, 3, byrow = T,
                 dimnames = list(originNames, originNames))
prop.table(shapes, 1)


calibSampleSize <- table(factor(rep(1:3, c(18, 22, 136)),
                                levels = 1:nTargetSites))

data2 <- simGenData(psi = psiTrue,
                        originRelAbund = originRelAbund,
                        sampleSize = calibSampleSize,
                        shapes = shapes,
                        captured = "target",
                        verbose = 1)
colMeans(data2$genProbs)
summary(data2$genProbs)
summary(originAssignmentPABU[-breeders, ])
mean(apply(data2$genProbs, 1, var))
mean(apply(originAssignmentPABU[-breeders, ], 1, var))
apply(data2$genProbs, 2, var)
apply(originAssignmentPABU[-breeders, ], 2, var)
mean(apply(data2$genProbs, 1, max)>0.99999)
mean(apply(originAssignmentPABU[-breeders, ], 1, max)>0.99999)
mean(apply(data2$genProbs, 1, max)>0.99)
mean(apply(originAssignmentPABU[-breeders, ], 1, max)>0.99)
