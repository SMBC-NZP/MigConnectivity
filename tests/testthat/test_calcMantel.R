library(MigConnectivity)
context('Calculate Mantel correlation')

test_that('Ovenbird example produces right Mantel correlation', {
  nAnimals <- nrow(OVENdata$targetPoints)
  originPoints2 <- sf::st_transform(OVENdata$originPoints,4326)
  originDist <- matrix(NA, nAnimals, nAnimals)

  originDist[lower.tri(originDist)] <- 1

  distIndices <- which(!is.na(originDist), arr.ind = TRUE)

  originDist0 <- geosphere::distGeo(sf::st_coordinates(originPoints2[distIndices[,'row'],]),
                                    sf::st_coordinates(originPoints2[distIndices[,'col'],]))

  originDist[lower.tri(originDist)] <- originDist0
  diag(originDist) <- 0
  originDist <- t(originDist)
  originDist[lower.tri(originDist)] <- originDist0

  targetDist <- matrix(NA, nAnimals, nAnimals)

  targetDist[lower.tri(targetDist)] <- 1

  distIndices <- which(!is.na(targetDist), arr.ind = TRUE)

  # project target points to WGS #
  targetPoints2 <- sf::st_transform(OVENdata$targetPoints,4326)

  targetDist0 <- geosphere::distGeo(sf::st_coordinates(targetPoints2[distIndices[,'row'],]),
                                    sf::st_coordinates(targetPoints2[distIndices[,'col'],]))

  targetDist[lower.tri(targetDist)] <- targetDist0
  diag(targetDist) <- 0
  targetDist <- t(targetDist)
  targetDist[lower.tri(targetDist)] <- targetDist0

  rM0 <- calcMantel(originPoints = OVENdata$originPoints, # Capture Locations
                    targetPoints = OVENdata$targetPoints) # Target locations

  expect_equal(rM0$pointCorr,
               0.813679665)
  expect_equal(lower.tri(rM0$originDist),
               lower.tri(originDist))
  expect_equal(lower.tri(rM0$targetDist),
               lower.tri(targetDist))
  expect_equal(calcMantel(originDist = originDist,
                          targetDist = targetDist)$pointCorr,
               0.813679665)
})
