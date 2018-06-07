# Function for generating random transition probability matrix (expression of parametric uncertainty)
makePsiRand <- function(model, origin.set, target.set) {
  vcv1 <- model$results$beta.vcv
  beta1 <- model$results$beta$estimate
  new.beta <- MASS::mvrnorm(1, beta1, vcv1)
  full.beta <- rep(NA, length(model$results$beta$estimate))
  new.real <- RMark::get.real(model, "Psi", new.beta,
                              design = model$design.matrix, vcv=F, se=F)
  mat <- RMark::TransitionMatrix(new.real)[origin.set, target.set]
  return(mat)
}

# Calculates probability matrix based on exponential decline with distance
mlogitMat <- function(slope, dist) {
  preMat <- exp(-slope/mean(dist)*dist)
  diag(preMat) <- 0
  nr <- nrow(dist)
  nc <- ncol(dist)
  outMat <- matrix(0, nr, nc)
  for (b in 1:nr) {
    outMat[b,] <- preMat[b,]/(1+sum(preMat[b, ]))
    outMat[b,b] <- 1 - sum(outMat[b, ])
  }
  return(outMat)
}

# Calculates row of probability matrix based on exponential decline with distance
mlogitRow <- function(slope, dist, row.num) {
  preMat <- exp(-slope/mean(dist)*dist)
  diag(preMat) <- 0
  out.row <- preMat[row.num,]/(1+sum(preMat[row.num, ]))
  out.row[row.num] <- 1 - sum(out.row)
  return(out.row)
}

# Crude optimizable function for developing MC pattern based on MC strength
mlogitMC <- function(slope, MC.in, origin.dist, target.dist, origin.rel.abund) {
  nBreeding <- nrow(origin.dist)
  nWintering <- nrow(target.dist)
  psi <- mlogitMat(slope, origin.dist)
  if (any(psi<0))
    return(5*slope^2)
  MC <- calcMC(origin.dist, target.dist, psi, origin.rel.abund)
  return((MC.in - MC)^2)
}

# Samples from GL and/or GPS locations
# Called by estMantel and estMCGlGps
targetSample <- function(isGL, geoBias, geoVCov, targetPoints, animal.sample,
                         nSim = 1000, targetSites = NULL, targetAssignment = NULL,
                         resampleProjection = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs",
                         maxTries = 300) {
  nAnimals <- length(targetPoints)

  # targetDist1 <- matrix(NA, nAnimals, nAnimals)
  #
  # targetDist1[lower.tri(targetDist1)] <- 1
  #
  # distIndices <- which(!is.na(targetDist1), arr.ind = T)
  #
  # # project target points to WGS #
  # targetPoints2 <- sp::spTransform(targetPoints, sp::CRS(WGS84))
  #
  # targetDist0 <- geosphere::distVincentyEllipsoid(targetPoints2[distIndices[,'row'],],
  #                                                 targetPoints2[distIndices[,'col'],])
  #
  # targetDist1[lower.tri(targetDist1)] <- targetDist0

  if (is.null(targetAssignment)) {
    if (!is.null(targetSites))
      targetAssignment <- sp::over(targetPoints, y = targetSites)
    else
      targetAssignment <- rep(0, nAnimals)
  }

  # Storage
  target.sample <- rep(NA, nAnimals)
  target.point.sample <- matrix(NA, nAnimals, 2)
  # Fill in GPS values (no location uncertainty)
  target.sample[which(!isGL[animal.sample])] <- targetAssignment[animal.sample[which(!isGL[animal.sample])]]
  target.point.sample[which(!isGL[animal.sample]), ] <- targetPoints[animal.sample[which(!isGL[animal.sample])],]@coords
  # Which to sample still
  toSample <- which(isGL[animal.sample])
  # In this case, only have to sample each GL animal once, because no restrictions
  if (is.null(targetSites) && any(isGL[animal.sample])) {
    draws <- 1
    geoBias2 <- array(rep(geoBias, length(toSample)), c(2, length(toSample)))
    point.sample <- array(apply(targetPoints@coords[animal.sample[toSample], , drop = FALSE], 1,
                                MASS::mvrnorm, n=1, Sigma=geoVCov),
                          c(2, length(toSample))) - geoBias2
    point.sample <- sp::SpatialPoints(point.sample, proj4string = sp::CRS(resampleProjection))
    target.point.sample[toSample, ]<- t(point.sample@coords)
  }
  else {
    draws <- 0
    while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
      draws <- draws + 1
      # Create giant version of geoBias we can subtract from point.sample
      geoBias2 <- array(rep(geoBias, length(toSample), each = nSim), c(nSim, 2, length(toSample)))
      # Draw random points for GL animals from multivariate normal distribution
      point.sample <- array(apply(targetPoints@coords[animal.sample[toSample], , drop = FALSE], 1,
                                  MASS::mvrnorm, n=nSim, Sigma=geoVCov),
                            c(nSim, 2, length(toSample))) - geoBias2
      # Convert those to SpatialPoints
      point.sample <- apply(point.sample, 3, sp::SpatialPoints,
                            proj4string = sp::CRS(resampleProjection))
      # Find out which sampled points are in a target site
      target.sample0 <- sapply(point.sample, sp::over, y = targetSites)
      # Identify which animals have at least one valid sample point. good.sample
      # will be NA for those that don't.  For those that do, it will location in
      # point.sample first valid point can be found.
      good.sample <- apply(target.sample0, 2, function(x) which(!is.na(x))[1])
      # Fill in target site of valid sampled points
      target.sample[toSample] <- apply(target.sample0, 2, function(x) x[!is.na(x)][1])
      # Fill in target points of valid sampled points
      if (any(!is.na(good.sample)))
        target.point.sample[toSample[!is.na(good.sample)], ]<- t(mapply(function(x, y) y[x]@coords,
                                                                        good.sample[!is.na(good.sample)], point.sample[!is.na(good.sample)]))
      toSample <- which(is.na(target.sample))
    }
    if (!is.null(maxTries) && draws > maxTries)
      stop(paste0('maxTries (',maxTries,') reached during point resampling, exiting. Examine targetSites, geoBias, and geoVcov to determine why so few resampled points fall within targetSites.'))
  }
  return(list(target.sample = target.sample, target.point.sample = target.point.sample,
              draws = draws))
}

# Samples from isotope target points
# Called by estMCisotope
targetSampleIsotope <- function(targetPoints, animal.sample,
                         nSim = 10, targetSites = NULL,
                         resampleProjection = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs",
                         maxTries = 300, pointsAssigned = FALSE,
                         overlay = NULL) {
  # Here is a catch to save Mike from himself
  #if(is.null(nSim)){nSim<-1000}
  if (pointsAssigned) {
    nAnimals <- dim(targetPoints$SingleCell)[3]
    # DETERMINE THE NUMBER OF RANDOM DRAWS FROM ISOSCAPE used in isoAssign
    nrandomDraws <- dim(targetPoints$SingleCell)[1]

    # Stack 3D array into 2D, to make it easier to get different random point samples for each animal
    targetPoints1 <- array(aperm(targetPoints$SingleCell, c(1, 3, 2)), dim = c(nAnimals * nrandomDraws, 2))

    # project target points to WGS #
    targetPoints2 <- sp::SpatialPoints(targetPoints1, proj4string = sp::CRS(projections$WGS84))
  }
  else {
    nAnimals <- dim(targetPoints$probassign)[3]
  }
  target.sample <- rep(NA, nAnimals)
  target.point.sample <- matrix(NA, nAnimals, 2)
  toSample <- 1:nAnimals
  # In this case, only need to draw once, as no limits on where points can be
  if (pointsAssigned) {
    draws <- 1
    samp <- sample.int(nrandomDraws, size = length(toSample), replace = T)
    samp2 <- samp + (toSample - 1) * nrandomDraws
    point.sample <- targetPoints2[samp2]
    # Changed to make sure x,y coords stack correctly
    target.point.sample[toSample,1]<- point.sample@coords[,1]
    target.point.sample[toSample,2]<- point.sample@coords[,2]
  }
  else if (is.null(targetSites)) {
    draws <- 1

  }
  else {
    draws <- 0
    # Make sure targetSites are WGS84
    targetSites <- sp::spTransform(targetSites, sp::CRS(projections$WGS84))

    while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
      draws <- draws + 1
      # Select nSim points for each animal still to be sampled
      if(!exists("nrandomDraws")){nrandomDraws <- 1000}
      samp <- sample.int(nrandomDraws, size = length(toSample) * nSim, replace = T)
      samp2 <- samp + rep(animal.sample[toSample] - 1, each = nSim) * nrandomDraws
      point.sample <- targetPoints2[samp2]
      # Check which points are in target sites
      target.sample0 <- sp::over(point.sample, y = targetSites)
      # Organize into matrix (separate by animal)
      target.sample1 <- matrix(target.sample0, nSim, length(toSample))
      # Identify which animals have a valid point (inside a target site).
      # good.sample, for those that don't, will be NA. For those that do, it
      # will be location in point.sample where first valid point is.
      good.sample <- apply(target.sample1, 2, function(x) which(!is.na(x))[1]) +
        seq(from = 0, by = nSim, length.out = length(toSample))
      # Put in target sites of animals with valid points sampled
      target.sample[toSample] <- apply(target.sample1, 2, function(x) x[!is.na(x)][1])
      # Put in target points of animals with valid points sampled
      if (any(!is.na(good.sample))) {
        target.point.sample[toSample[!is.na(good.sample)],1]<- point.sample[good.sample[!is.na(good.sample)]]@coords[,1]
        target.point.sample[toSample[!is.na(good.sample)],2]<- point.sample[good.sample[!is.na(good.sample)]]@coords[,2]
      }
      toSample <- which(is.na(target.sample))
    }
    if (!is.null(maxTries) && draws > maxTries)
      stop(paste0('maxTries (',maxTries,') reached during point resampling, exiting. Examine targetSites, geoBias, and geoVcov to determine why so few resampled points fall within targetSites.'))
  }
  return(list(target.sample = target.sample, target.point.sample = target.point.sample,
              draws = draws))
}


#' Distance matrix from position matrix
#'
#' @param pos Number of sites by 2 matrix with postions of each site.  If
#'    \code{surface} is 'ellipsoid' or 'sphere', then column 1 should be
#'    longitude and column 2 should be latitude.  If \code{surface} is 'plane',
#'    column 1 can be x-position and column 2 y-position.
#' @param surface Surface to calculate distances on.  Either 'ellipsoid'
#'    (default), 'sphere', or 'plane'.
#'
#' @return Square matrix of distances between sites. If \code{surface} is
#'    'ellipsoid' or 'sphere', then units will be km; if \code{surface} is
#'    'plane', the units will be the same as the \code{pos} units.
#' @export
#'
#' @examples
#' nBreeding <- 100
#' nWintering <- 100
#' breedingPos <- matrix(c(rep(seq(-99, -81, 2), each = sqrt(nBreeding)),
#'                         rep(seq(49, 31, -2), sqrt(nBreeding))),
#'                       nBreeding, 2)
#' winteringPos <- matrix(c(rep(seq(-79, -61, 2), each = sqrt(nWintering)),
#'                          rep(seq(9, -9, -2), sqrt(nWintering))),
#'                        nWintering, 2)
#' head(breedingPos)
#' tail(breedingPos)
#' head(winteringPos)
#' tail(winteringPos)
#'
#' breedDist <- distFromPos(breedingPos, 'ellipsoid')
#' nonbreedDist <- distFromPos(winteringPos, 'ellipsoid')
#' breedDist[1:12, 1:12]
#' breedDist[1:12, c(1,91,100)]
distFromPos <- function(pos, surface = 'ellipsoid') {
  if (!(surface %in% c('plane', 'sphere', 'ellipsoid')))
    stop('surface must be "plane", "sphere", or "ellipsoid".')
  nSites <- nrow(pos)
  dist <- matrix(0, nSites, nSites)
  for (b1 in 2:nSites) {
    for (b2 in 1:(b1-1)) {
      if (surface == 'ellipsoid')
        dist[b1, b2] <- geosphere::distVincentyEllipsoid(pos[b1, ],
                                                         pos[b2, ]) / 1000
      else if (surface == 'sphere')
        dist[b1, b2] <- geosphere::distVincentySphere(pos[b1, ],
                                                      pos[b2, ]) / 1000
      else
        dist[b1, b2] <- sqrt((pos[b1, 1] - pos[b2, 1]) ^ 2 +
                              (pos[b1, 2] - pos[b2, 2]) ^ 2)
      dist[b2, b1] <- dist[b1, b2]
    }
  }
  return(dist)
}
