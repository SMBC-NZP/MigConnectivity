###############################################################################
# Input values
###############################################################################

#########################
#Input values 1 of 3
# Eight transition probability scenarios
#########################
nScenarios1 <- length(samplePsis)
MC1 <- rep(NA, nScenarios1)
for (i in 1:nScenarios1) {
  MC1[i] <- calcMC(sampleOriginDist[[1]], sampleTargetDist[[1]], samplePsis[[i]], sampleOriginRelN[[1]])
}
names(MC1) <- names(samplePsis)
round(MC1, 6)

#########################
#Input values 2 of 3
# 12 spatial arrangements that result in different distances between regions
# Distance scenarios
## A) Base distances, linear/ linear    1
## B) Distance between breeding sites 2 and 3 doubled
## C) Distance between breeding sites 2 and 3 halved
## D) Distance between breeding sites 3 and 4 doubled
## E) Distance between breeding sites 3 and 4 halved
## F) Breeding sites on square grid/ winter linear   6
## G) Distance between wintering sites 2 and 3 doubled
## H) Distance between wintering sites 2 and 3 halved
## I) Distance between wintering sites 3 and 4 doubled
## J) Distance between wintering sites 3 and 4 halved
## K) Breeding linear, Wintering sites on square grid
## L) Wintering and breeding on square grid  12
#########################

# Get MC strengths
nScenarios2 <- length(sampleOriginPos)
MC2 <- matrix(NA, nScenarios1, nScenarios2)
rownames(MC2) <- names(samplePsis)
colnames(MC2) <- names(sampleOriginPos)
for (i in 1:nScenarios1) {
  for (j in 1:nScenarios2) {
    MC2[i, j] <- calcMC(sampleOriginDist[[j]], sampleTargetDist[[j]], samplePsis[[i]], sampleOriginRelN[[1]])
  }
}

t(round(MC2, 4))

# Different way of comparing results
MC.diff2 <- apply(MC2, 2, "-", MC2[ , 1])
t(round(MC.diff2, 4))

#########################
#Input values 3 of 3
# Changes to relative breeding abundance:
#   1. Base
#   2. Abundance at site B doubled
#   3. Abundance at site B halved
#   4. Abundance at site D doubled
#   5. Abundance at site D halved
# For all eight transition probability matrices and three distance scenarios
#########################

nScenarios3 <- length(sampleOriginRelN)

# Get MC strengths for breeding linear/ winter linear arrangement
MC3 <- matrix(NA, nScenarios1, nScenarios3)
rownames(MC3) <- names(samplePsis)
colnames(MC3) <- names(sampleOriginRelN)
for (i in 1:nScenarios1) {
  for (j in 1) {
    for (k in 1:nScenarios3) {
      MC3[i, k] <- calcMC(sampleOriginDist[[j]], sampleTargetDist[[j]], samplePsis[[i]], sampleOriginRelN[[k]])
    }
  }
}
t(round(MC3, 4)) # linear arrangement

# Get MC strengths for breeding grid/ winter grid arrangement
MC4 <- matrix(NA, nScenarios1, nScenarios3)
rownames(MC4) <- names(samplePsis)
colnames(MC4) <- names(sampleOriginRelN)
for (i in 1:nScenarios1) {
  for (j in nScenarios2) {
    for (k in 1:nScenarios3) {
      MC4[i, k] <- calcMC(sampleOriginDist[[j]], sampleTargetDist[[j]], samplePsis[[i]], sampleOriginRelN[[k]])
    }
  }
}
t(round(MC4, 4)) # grid arrangement

# Get MC strengths for breeding grid, winter linear arrangement
MC5 <- matrix(NA, nScenarios1, nScenarios3)
rownames(MC5) <- names(samplePsis)
colnames(MC5) <- names(sampleOriginRelN)
for (i in 1:nScenarios1) {
  for (j in 6) {
    for (k in 1:nScenarios3) {
      MC5[i, k] <- calcMC(sampleOriginDist[[j]], sampleTargetDist[[j]], samplePsis[[i]], sampleOriginRelN[[k]])
    }
  }
}
t(round(MC5, 4)) # breeding grid, winter linear arrangement

