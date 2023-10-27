library(MigConnectivity)
context('Calculate Mantel correlation')

test_that('Ovenbird example produces right Mantel correlation', {
  nAnimals <- nrow(OVENdata$targetPoints)
  originPoints2 <- sf::st_transform(OVENdata$originPoints,4326)
  originDist <- geodist::geodist(sf::st_coordinates(originPoints2),
                                 measure = "geodesic") / 1000

  # project target points to WGS #
  targetPoints2 <- sf::st_transform(OVENdata$targetPoints,4326)

  targetDist <- geodist::geodist(sf::st_coordinates(targetPoints2),
                                 measure = "geodesic") / 1000

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
