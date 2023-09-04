# Generate values for input value simulations

#########################
# Eight Psi scenarios
# Breeding as rows, Non-breeding as columns
#########################
samplePsis <- list(matrix(0.25, 4, 4), #"Full Mix"
                matrix(c(rep(0.32, 12), rep(0.04, 4)), 4, 4), #"Avoid One Site"
                diag(nrow=4),  #"Full Connectivity"
                matrix(c(rep(c(0.5, 0), 2, each=2), rep(c(0, 0.5), 2, each=2)),
                       4, 4), #"Half Mix"
                matrix(c(0.55, 0.2, 0.15, 0.1, 0.1, 0.55, 0.2, 0.15, 0.15, 0.1,
                         0.55, 0.2, 0.2, 0.15, 0.1, 0.55), 4, 4, byrow=TRUE), #"Low"
                matrix(c(rep(c(0.75, 0.15, rep(0.05, 3)), 3), 0.75), 4, 4,
                       byrow=TRUE), #"Medium"
                matrix(c(rep(0.25, 12), rep(0, 3), 1), 4, 4, byrow=TRUE), #"Site Pref"
                matrix(c(0.01, 0.49, 0.49, 0.01, 0.49, 0.01, 0.01, 0.49, 0.49,
                         0.01, 0.01, 0.49, 0.01, 0.49, 0.49, 0.01), 4, 4,
                       byrow=TRUE)) #Negative
samplePsis[[6]][4, 1] <- 0.15
psi.scenario.names <- c("Full Mix", "Avoid One Site", "Full Connectivity",
                        "Half Mix", "Low", "Medium", "One Site Preference",
                        "Negative")
names(samplePsis) <- psi.scenario.names
samplePsis <- lapply(samplePsis, provideDimnames, base = list(LETTERS[1:4], as.character(1:4)))
devtools::use_data(samplePsis, overwrite = TRUE)
nOrigin <- sapply(samplePsis, nrow)
nTarget <- sapply(samplePsis, ncol)
nScenariosPsi <- length(samplePsis)

#########################
#Input values 2 of 3
# 12 spatial arrangements that result in different distances between regions
# Distance scenarios
## Distance scenarios
## A) Base distances, linear/ linear    1
## B) Distance between origin sites 2 and 3 doubled
## C) Distance between origin sites 2 and 3 halved
## D) Distance between origin sites 3 and 4 doubled
## E) Distance between origin sites 3 and 4 halved
## F) Origin sites on square grid/ winter linear   6
## G) Distance between target sites 2 and 3 doubled
## H) Distance between target sites 2 and 3 halved
## I) Distance between target sites 3 and 4 doubled
## J) Distance between target sites 3 and 4 halved
## K) Origin linear, Target sites on square grid
## L) Target and origin on square grid  12
#########################

# 2D plane x/y positions of sites
genericOriginPos4 <- matrix(c(1:4, rep(4, 4)), 4, 2)
genericTargetPos4 <- matrix(c(1:4, rep(-4, 4)), 4, 2)
sampleOriginPos <- list(genericOriginPos4,
                          matrix(c(0.5, 1.5, 3.5, 4.5, rep(4, 4)), 4, 2),
                          matrix(c(0.75, 1.75, 2.25, 3.25, rep(4, 4)), 4, 2),
                          matrix(c(0.5, 1.5, 2.5, 4.5, rep(4, 4)), 4, 2),
                          matrix(c(0.75, 1.75, 2.75, 3.25, rep(4, 4)), 4, 2),
                          matrix(c(2, 2, 3, 3, 4, 3, 4, 3), 4, 2),
                          genericOriginPos4, genericOriginPos4,
                          genericOriginPos4, genericOriginPos4,
                          genericOriginPos4,
                          matrix(c(2, 2, 3, 3, 4, 3, 4, 3), 4, 2))
sampleTargetPos <- list(genericTargetPos4, genericTargetPos4,
                        genericTargetPos4, genericTargetPos4,
                        genericTargetPos4, genericTargetPos4,
                        matrix(c(0.5, 1.5, 3.5, 4.5, rep(-4, 4)), 4, 2),
                        matrix(c(0.75, 1.75, 2.25, 3.25, rep(-4, 4)), 4, 2),
                        matrix(c(0.5, 1.5, 2.5, 4.5, rep(-4, 4)), 4, 2),
                        matrix(c(0.75, 1.75, 2.75, 3.25, rep(-4, 4)), 4, 2),
                        matrix(c(2, 2, 3, 3, -4, -3, -4, -3), 4, 2),
                        matrix(c(2, 2, 3, 3, -4, -3, -4, -3), 4, 2))
nScenariosPos <- length(sampleOriginPos)

# 2D plane distances of sites
sampleOriginDist <- sampleTargetDist <- vector("list", nScenariosPos)
for (j in 1:nScenariosPos) {
  sampleOriginDist[[j]] <- sampleTargetDist[[j]] <- matrix(0, 4, 4)
  rownames(sampleOriginDist[[j]]) <- c('A','B','C','D')
  colnames(sampleOriginDist[[j]]) <- c('A','B','C','D')
  rownames(sampleOriginPos[[j]]) <- c('A','B','C','D')
  for (b1 in 2:4) {
    for (b2 in 1:(b1-1)) {
      sampleOriginDist[[j]][b1,b2] <- sqrt((sampleOriginPos[[j]][b1,1] -
                                              sampleOriginPos[[j]][b2,1])^2 +
                                             (sampleOriginPos[[j]][b1,2] -
                                                sampleOriginPos[[j]][b2,2])^2)
      sampleOriginDist[[j]][b2,b1] <- sampleOriginDist[[j]][b1,b2]
    }
  }
  for (w1 in 2:4) {
    for (w2 in 1:(w1-1)) {
      sampleTargetDist[[j]][w1,w2] <- sqrt((sampleTargetPos[[j]][w1,1] -
                                              sampleTargetPos[[j]][w2,1])^2 +
                                             (sampleTargetPos[[j]][w1,2] -
                                                sampleTargetPos[[j]][w2,2])^2)
      sampleTargetDist[[j]][w2,w1] <- sampleTargetDist[[j]][w1,w2]
    }
  }
}
sampleOriginDist
sampleTargetDist
pos.scenario.names <- c("Linear", "B Dist BC*2", "B Dist BC/2",
                        "B Dist CD*2", "B Dist CD/2", "B Grid", "NB Dist 23*2",
                        "NB Dist 23/2", "NB Dist 34*2", "NB Dist 34/2", "NB Grid",
                        "B/NB Grid")
names(sampleTargetDist) <- names(sampleTargetPos) <- pos.scenario.names
names(sampleOriginDist) <- names(sampleOriginPos) <- pos.scenario.names
pos.scenario.names
devtools::use_data(sampleOriginPos, overwrite = TRUE)
devtools::use_data(sampleOriginDist, overwrite = TRUE)
devtools::use_data(sampleTargetPos, overwrite = TRUE)
devtools::use_data(sampleTargetDist, overwrite = TRUE)

#########################
#Input values 3 of 3
# Changes to relative origin abundance:
#   1. Base
#   2. Abundance at site B doubled
#   3. Abundance at site B halved
#   4. Abundance at site D doubled
#   5. Abundance at site D halved
# For all eight transition probability matrices and three distance scenarios
#########################

sampleOriginN <- list(rep(1000, 4), c(1000, 2000, 1000, 1000),
                      c(1000, 500, 1000, 1000), c(1000, 1000, 1000, 2000),
                      c(1000, 1000, 1000, 500))
nScenariosRelN <- length(sampleOriginN)

sampleOriginRelN <- vector("list", nScenariosRelN)
for (k in 1:nScenariosRelN) {
  sampleOriginRelN[[k]] <- sampleOriginN[[k]]/sum(sampleOriginN[[k]])
}

abund.scenario.names <- c("Base", "B Doub", "B Half", "D Doub", "D Half")
names(sampleOriginRelN) <- names(sampleOriginN) <- abund.scenario.names
sampleOriginN
sampleOriginRelN
sampleTotalN <- sapply(sampleOriginN, sum)
sampleTotalN
devtools::use_data(sampleOriginN, overwrite = TRUE)
devtools::use_data(sampleOriginRelN, overwrite = TRUE)

