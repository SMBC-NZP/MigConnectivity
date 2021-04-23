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
targetSample <- function(isGL,
                         geoBias,
                         geoVCov,
                         targetPoints,
                         animal.sample,
                         nSim = 1000,
                         targetSites = NULL,
                         targetAssignment = NULL,
                         resampleProjection = 'ESRI:53027',
                         maxTries = 300) {

  # Target points should come in as sf object #
  # length(targetPoints) works for sp objects #
  # nrow(targetPoints) works for sf objects #
  nAnimals <- nrow(targetPoints)

  if (is.null(targetAssignment)) {
    if (!is.null(targetSites)){
      #Use the provided spatial layers to create targetAssignment
      if(!st_crs(targetSites)==st_crs(targetPoints)){targetPoints <- sf::st_transform(targetPoints, crs = resampleProjection)}
      targetAssignment <- as.numeric(unclass(sf::st_intersects(x = targetPoints, y = targetSites, sparse = TRUE)))
      #if targetAssignment is NA - convert to 0
      targetAssignment[is.na(targetAssignment)] <- 0
    }else
      targetAssignment <- rep(0, nAnimals)
  }

  # Storage
  target.sample <- rep(NA, nAnimals)
  target.point.sample <- matrix(NA, nAnimals, 2)
  # Fill in GPS values (no location uncertainty)
  target.sample[which(!isGL[animal.sample])] <- targetAssignment[animal.sample[which(!isGL[animal.sample])]]
  target.point.sample[which(!isGL[animal.sample]), ] <- sf::st_coordinates(targetPoints[animal.sample[which(!isGL[animal.sample])],])
  # Which to sample still
  toSample <- which(isGL[animal.sample])
  # In this case, only have to sample each GL animal once, because no restrictions
  if (is.null(targetSites) && any(isGL[animal.sample])) {
    draws <- 1
    geoBias2 <- array(rep(geoBias, length(toSample)), c(2, length(toSample)))

    # generate sample and substract the bias
    point.sample <- array(apply(sf::st_coordinates(targetPoints)[animal.sample[toSample], , drop = FALSE],
                                MARGIN = 1,
                                MASS::mvrnorm, n=1, Sigma=geoVCov),
                          # dimensions of the array
                          c(2, length(toSample))) - geoBias2 #subtract the bias

    # name so can
    rownames(point.sample) <- c("x","y")

    point.sample <- sf::st_as_sf(data.frame(t(point.sample)), coords = c("x","y"), crs = resampleProjection)

    target.point.sample[toSample, ]<- t(sf::st_coordinates(point.sample))
  }else {
    draws <- 0
    while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
      draws <- draws + 1
      # Create giant version of geoBias we can subtract from point.sample
      geoBias2 <- array(rep(geoBias, length(toSample), each = nSim), c(nSim, 2, length(toSample)))
      # Draw random points for GL animals from multivariate normal distribution
      point.sample <- array(apply(sf::st_coordinates(targetPoints)[animal.sample[toSample], , drop = FALSE], 1,
                                  MASS::mvrnorm, n=nSim, Sigma=geoVCov),
                            c(nSim, 2, length(toSample))) - geoBias2

      dimnames(point.sample)[[2]] <- c("x","y")

      # Convert those to SpatialPoints

      point.sample <- apply(point.sample, FUN = function(x){sf::st_as_sf(data.frame(x), coords = c("x","y"), crs = resampleProjection)}, MARGIN = 3)
      #point.sample <- lapply(point.sample1, st_as_sf)
      # Find out which sampled points are in a target site
      if(!st_crs(targetSites)==st_crs(point.sample)){targetSites <- sf::st_transform(targetSites, crs = resampleProjection)}

      target.sample0 <- sapply(point.sample, FUN = function(z){as.numeric(unclass(sf::st_intersects(x = z, y = targetSites, sparse = TRUE)))})

      # Identify which animals have at least one valid sample point. good.sample
      # will be NA for those that don't.  For those that do, it will location in
      # point.sample first valid point can be found.
      good.sample <- apply(target.sample0, 2, function(x) which(!is.na(x))[1])
      # Fill in target site of valid sampled points
      target.sample[toSample] <- apply(target.sample0, 2, function(x) x[!is.na(x)][1])
      # Fill in target points of valid sampled points
      if (any(!is.na(good.sample)))
        target.point.sample[toSample[!is.na(good.sample)], ]<- t(mapply(x = good.sample[!is.na(good.sample)],
                                                                        y = point.sample[!is.na(good.sample)],
                                                                        FUN = function(x, y){st_coordinates(y[x,])}))

      toSample <- which(is.na(target.sample))
    }
    if (!is.null(maxTries) && draws > maxTries)
      stop(paste0('maxTries (',maxTries,') reached during point resampling, exiting. Examine targetSites, geoBias, and geoVcov to determine why so few resampled points fall within targetSites.'))
  }
  return(list(target.sample = target.sample, target.point.sample = target.point.sample,
              draws = draws))
}


locSample <- function(isGL,
                      isRaster,
                      isProb,
                      isTelemetry,
                      geoBias = NULL,
                      geoVCov = NULL,
                      points = NULL,
                      raster1 = NULL,
                      singleCell = NULL,
                      pointsInSites = FALSE,
                      overlap1 = NULL,
                      assignment = NULL,
                      sites = NULL,
                      resampleProjection = 'ESRI:53027',
                      nSim = 1000,
                      maxTries = 300) {
  # points should come in as sf object #
  nAnimals <- length(isGL)
  if (!is.null(singleCell)) {
    # DETERMINE THE NUMBER OF RANDOM DRAWS FROM ISOSCAPE used in isoAssign
    nrandomDraws <- dim(targetIntrinsic$SingleCell)[1]

    # Stack 3D array into 2D, to make it easier to get different random point samples for each animal
    intrinsic1 <- array(aperm(singleCell, c(1, 3, 2)),
                        dim = c(nAnimals * nrandomDraws, 2),
                        dimnames = list(NULL, c("x", "y")))

    # project target points to WGS #
    intrinsic2 <- sf::st_as_sf(as.data.frame(intrinsic1), coords = c("x", "y"),
                               crs = MigConnectivity::projections$WGS84)
  }
  if (!pointsInSites & !is.null(raster1))
    # converts raster to matrix of XY then probs
    matvals <- raster::rasterToPoints(raster1)

  if (is.null(assignment)) {
    if (!is.null(sites)){
      #Use the provided spatial layers to create assignment
      if(!sf::st_crs(sites)==sf::st_crs(points)){
        points <- sf::st_transform(points, crs = resampleProjection)
      }
      assignment <- unclass(sf::st_intersects(x = points, y = sites,
                                                    sparse = TRUE))
      #if assignment is not there - convert to 0
      assignment[sapply(assignment, function(x) length(x)==0)] <- 0
      assignment <- array(unlist(assignment), nAnimals)
    }
    else
      assignment <- array(0, nAnimals)
  }

  # Storage
  site.sample <- rep(NA, nAnimals)
  point.sample <- matrix(NA, nAnimals, 2)
  # Fill in GPS values (no location uncertainty)
  if (length(dim(assignment))==2){
    site.sample[which(isTelemetry)] <- apply(assignment[isTelemetry, ], 1,
                                                  which.max)
  }else{
    site.sample[which(isTelemetry)] <- assignment[which(isTelemetry)]}
  if (!is.null(points))
    point.sample[which(isTelemetry), ] <- sf::st_coordinates(points[which(isTelemetry),])
  # Which to sample still
  toSample <- which(!isTelemetry)
  # In this case, only have to sample each GL animal once, because no restrictions
  # if (is.null(sites) && any(isGL)) {
  #   draws <- 1
  #   geoBias2 <- array(rep(geoBias, length(toSample)), c(2, length(toSample)))
  #
  #   # generate sample and subtract the bias
  #   point.sample0 <- array(apply(sf::st_coordinates(points)[toSample, ,
  #                                                           drop = FALSE],
  #                                MARGIN = 1,
  #                                MASS::mvrnorm, n=1, Sigma=geoVCov),
  #                          # dimensions of the array
  #                          c(2, length(toSample))) - geoBias2 #subtract the bias
  #
  #   # name so can
  #   rownames(point.sample0) <- c("x","y")
  #
  #   point.sample0 <- sf::st_as_sf(data.frame(t(point.sample0)), coords = c("x","y"), crs = resampleProjection)
  #
  #   point.sample[toSample, ]<- t(sf::st_coordinates(point.sample0))
  # }else {
  draws <- 0
  while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
    draws <- draws + 1
    toSampleBool <- 1:nAnimals %in% toSample
    if (any(isGL & toSampleBool)) {
      # Create giant version of geoBias we can subtract from point.sample
      geoBias2 <- array(rep(geoBias, sum(isGL & toSampleBool), each = nSim),
                        c(nSim, 2, sum(isGL & toSampleBool)))
      # Draw random points for GL animals from multivariate normal distribution
      point.sample0 <- array(apply(sf::st_coordinates(points)[toSampleBool & isGL,
                                                              , drop = FALSE], 1,
                                  MASS::mvrnorm, n=nSim, Sigma=geoVCov),
                            c(nSim, 2, sum(isGL & toSampleBool))) - geoBias2

      dimnames(point.sample0)[[2]] <- c("x","y")

      # Convert those to SpatialPoints

      point.sample0 <- apply(point.sample0,
                             FUN = function(x){sf::st_as_sf(data.frame(x),
                                                            coords = c("x","y"),
                                                            crs = resampleProjection)},
                             MARGIN = 3)
      #point.sample <- lapply(point.sample1, st_as_sf)
      # Find out which sampled points are in a target site
      if(!sf::st_crs(sites)==sf::st_crs(point.sample0[[1]])){
        sites <- sf::st_transform(sites, crs = resampleProjection)
      }

      site.sample0 <- sapply(point.sample0, FUN = function(z){
        as.numeric(unclass(sf::st_intersects(x = z, y = sites,
                                             sparse = TRUE)))})

      # Identify which animals have at least one valid sample point. good.sample0
      # will be NA for those that don't.  For those that do, it will location in
      # point.sample first valid point can be found.
      good.sample0 <- apply(site.sample0, 2, function(x) which(!is.na(x))[1])
    }
    if (any(isProb & toSampleBool)) {
      #print(assignment)
      site.sample1 <- apply(assignment[isProb & toSampleBool, , drop = FALSE],
                            1, sample.int, n = ncol(assignment), size = nSim,
                            replace = TRUE)
    }
    if (any(isRaster & toSampleBool)) {
      if (pointsInSites && !any(isRaster & toSampleBool & (isGL | isProb))) {
        samp <- sample.int(nrandomDraws, size = sum(isRaster & toSampleBool),
                           replace = T)
        samp2 <- samp + (toSample - 1) * nrandomDraws
        point.sample2 <- intrinsic2[samp2, ]
        # Changed to make sure x,y coords stack correctly
        # target.point.sample[toSample,1]<- point.sample@coords[,1]
        # target.point.sample[toSample,2]<- point.sample@coords[,2]
        if (!is.null(overlap1))
          # Grab the relevant sites
          site.sample2 <- overlap1[samp2]
        good.sample2 <- rep(1, sum(isRaster & toSampleBool))
      }
      else if (is.null(singleCell)) {
        # This draws samples nSamples per animal (faster than looping over nSamples) and fills the xysim with x,y coords
        for(i in 3:ncol(matvals)) {
          animals <- which(toSample == i - 2)
          if (length(animals) > 0) {
            multidraw <- rmultinom(n = length(animals)*nSim, size = 1, prob = matvals[,i])
            point.sample2a <- matvals[which(multidraw == 1, arr.ind = TRUE)[, 1], 1:2]
            point.sample2b <- sf::st_as_sf(as.data.frame(point.sample2a),
                                           crs = MigConnectivity::projections$WGS84)
            # Check which points are in target sites
            if (!is.null(sites)) {
              site.sample2a <- sf::st_intersects(point.sample2b, y = sites)
              good.sample2 <- which(!is.na(target.sample0))
              nToFill <- min(length(good.sample2), length(animals))
              if (nToFill > 0) {
                site.sample2[animals[1:nToFill]] <- site.sample2a[good.sample2[1:nToFill]]
                point.sample2[animals[1:nToFill], 1] <- point.sample2a[good.sample2[1:nToFill], 1]
                point.sample2[animals[1:nToFill], 2] <- point.sample2a[good.sample2[1:nToFill], 2]
              }
            }
          }
        }
      }
      else {
        # Select nSim points for each animal still to be sampled
        samp <- sample.int(nrandomDraws, size = sum(isRaster & toSampleBool) * nSim,
                           replace = T)
        samp2 <- samp + rep(toSample - 1, each = nSim) * nrandomDraws
        point.sample2a <- intrinsic2[samp2, ]
        # Check which points are in target sites
        site.sample2a <- sf::st_intersects(point.sample, y = targetSites)
        # Organize into matrix (separate by animal)
        site.sample2b <- matrix(site.sample2a, nSim, sum(isRaster & toSampleBool))
        # Identify which animals have a valid point (inside a site).
        # good.sample2, for those that don't, will be NA. For those that do, it
        # will be location in point.sample where first valid point is.
        good.sample2 <- apply(site.sample2b, 2, function(x) which(!is.na(x))[1]) +
          seq(from = 0, by = nSim, length.out = sum(isRaster & toSampleBool))
      }
    }
    if (any(isProb & isGL & !isRaster & toSampleBool)){
      if (any(!is.na(good.sample0))){
        compare <- site.sample0[, which(isGL & toSampleBool) %in% which(isProb & !isRaster),
                                     drop = FALSE] ==
          site.sample1[, which(isProb & toSampleBool) %in% which(isGL & !isRaster),
                       drop = FALSE]
        good.sample01 <- apply(compare, 2, function(x) which(!is.na(x) & x)[1])
        if (any(!is.na(good.sample01))){
          for (col in 1:sum(isProb & isGL & !isRaster & toSampleBool))
            site.sample[which(isProb & isGL & !isRaster & toSampleBool)[col]] <-
              site.sample0[good.sample01[col],
                           which(which(isGL & toSampleBool) %in% which(isProb & !isRaster))[col]]
          # Fill in target points of valid sampled points
          point.sample[which(isProb & isGL & !isRaster & toSampleBool)[which(!is.na(good.sample01))], ]<- t(mapply(x = good.sample01[!is.na(good.sample01)],
                                                                   y = point.sample0[which(which(isGL & toSampleBool) %in% which(isProb & !isRaster))[which(!is.na(good.sample01))]],
                                                                   FUN = function(x, y) sf::st_coordinates(y[x,])))
        }
      }
    }
    if (any(isProb & !isGL & !isRaster & toSampleBool)){
      site.sample[isProb & !isGL & toSampleBool] <- site.sample1[1,
                                                                 !(which(isProb) %in% which(isGL | isRaster))]
    }
    if (any(!isProb & isGL & !isRaster & toSampleBool)){
      if (any(!is.na(good.sample0[!(which(isGL) %in% which(isProb | isRaster))]))){
        site.sample[!isProb & isGL & !isRaster & toSampleBool] <- apply(site.sample0[, !(which(isGL & toSampleBool) %in% which(isProb | isRaster)), drop = FALSE],
                                                            2,
                                                            function(x) x[!is.na(x)][1])
        # Fill in target points of valid sampled points
        point.sample[which(!isProb & isGL & !isRaster & toSampleBool)[which(!is.na(good.sample0[!(which(isGL & toSampleBool) %in% which(isProb | isRaster))]))], ]<-
          t(mapply(x = good.sample0[!is.na(good.sample0) & !(which(isGL & toSampleBool) %in% which(isProb | isRaster))],
                   y = point.sample0[which(!(which(isGL & toSampleBool) %in% which(isProb | isRaster)) & (!is.na(good.sample0)))],
                   FUN = function(x, y) sf::st_coordinates(y[x,])))
      }
    }
    if (any(!isProb & !isGL & isRaster & toSampleBool)){
     if (any(!is.na(good.sample2[!(which(isRaster) %in% which(isProb | isGL))]))){
        site.sample[which(!isProb & !isGL & isRaster & toSampleBool)] <-
          apply(site.sample2[, !(which(isRaster & toSampleBool) %in% which(isProb | isGL)),
                             drop = FALSE],
                2,
                function(x) x[!is.na(x)][1])
        # Fill in target points of valid sampled points
        point.sample[which(!isProb & !isGL & isRaster & toSampleBool)[which(!is.na(good.sample2))], ]<-
          t(mapply(x = good.sample12[!is.na(good.sample12)],
                   y = point.sample2[which(which(isRaster & toSampleBool) %in% which(!isProb & !isGL & !is.na(good.sample2)))],
                   FUN = function(x, y) sf::st_coordinates(y[x,])))
      }
    }
    if (any(isProb & !isGL & isRaster & toSampleBool)){
      if (any(!is.na(good.sample2))){
        compare <- site.sample2[, which(isRaster & toSampleBool) %in% which(isProb & !isGL),
                                drop = FALSE] ==
          site.sample1[, which(isProb & toSampleBool) %in% which(!isGL & isRaster),
                       drop = FALSE]
        good.sample12 <- apply(compare, 2, function(x) which(!is.na(x) & x)[1])
        if (any(!is.na(good.sample12))){
          for (col in 1:sum(isProb & !isGL & isRaster & toSampleBool))
            site.sample[which(isProb & !isGL & isRaster & toSampleBool)[col]] <-
              site.sample2[good.sample12[col],
                           which(which(isRaster & toSampleBool) %in% which(isProb & !isGL))[col]]
          # Fill in target points of valid sampled points
          point.sample[which(isProb & !isGL & isRaster & toSampleBool)[which(!is.na(good.sample12))], ]<-
            t(mapply(x = good.sample12[!is.na(good.sample12)],
                     y = point.sample2[which(which(isRaster & toSampleBool) %in% which(isProb & !isGL))[which(!is.na(good.sample12))]],
                     FUN = function(x, y) sf::st_coordinates(y[x,])))
        }
      }
    }
    if (any(!isProb & isGL & isRaster & toSampleBool)){
      if (any(!is.na(good.sample0)) && any(!is.na(good.sample2))){
        compare <- site.sample2[, which(isRaster & toSampleBool) %in% which(!isProb & isGL),
                                drop = FALSE] ==
          site.sample0[, which(isGL & toSampleBool) %in% which(!isProb & isRaster),
                       drop = FALSE]
        good.sample02 <- apply(compare, 2, function(x) which(!is.na(x) & x)[1])
        if (any(!is.na(good.sample02))){
          for (col in 1:sum(!isProb & isGL & isRaster & toSampleBool))
            site.sample[which(!isProb & isGL & isRaster & toSampleBool)[col]] <-
              site.sample2[good.sample02[col],
                           which(which(isRaster & toSampleBool) %in% which(!isProb & isGL))[col]]
          # Need to change this! Some options:
          # 1) Randomly choose between first GL and first raster point
          # 2) Take point halfway between them, unless that's not in site, then do 1
          # 3) Something more clever?
          point.sample[which(!isProb & isGL & isRaster & toSampleBool)[which(!is.na(good.sample02))], ]<-
            t(mapply(x = good.sample02[!is.na(good.sample02)],
                     y = point.sample2[which(which(isRaster & toSampleBool) %in% which(!isProb & isGL))[which(!is.na(good.sample02))]],
                     FUN = function(x, y) sf::st_coordinates(y[x,])))
        }
      }
    }
    if (any(isProb & isGL & isRaster & toSampleBool)){
      if (any(!is.na(good.sample0)) && any(!is.na(good.sample2))){
        compare <- (site.sample2[, which(isRaster & toSampleBool) %in% which(isProb & isGL),
                                drop = FALSE] ==
          site.sample0[, which(isGL & toSampleBool) %in% which(isProb & isRaster),
                       drop = FALSE]) &
          (site.sample1[, which(isProb & toSampleBool) %in% which(isRaster & isGL),
                        drop = FALSE] ==
             site.sample0[, which(isGL & toSampleBool) %in% which(isProb & isRaster),
                          drop = FALSE])
        good.sample012 <- apply(compare, 2, function(x) which(!is.na(x) & x)[1])
        if (any(!is.na(good.sample012))){
          for (col in 1:sum(isProb & isGL & isRaster & toSampleBool))
            site.sample[which(isProb & isGL & isRaster & toSampleBool)[col]] <-
              site.sample2[good.sample012[col],
                           which(which(isRaster & toSampleBool) %in% which(isProb & isGL))[col]]
          # This one might not need changed (not planning to include isProb in
          # estMantel or look at points in estTransition).
          point.sample[which(isProb & isGL & isRaster & toSampleBool)[which(!is.na(good.sample012))], ]<-
            t(mapply(x = good.sample012[!is.na(good.sample012)],
                     y = point.sample2[which(which(isRaster & toSampleBool) %in% which(isProb & isGL))[which(!is.na(good.sample012))]],
                     FUN = function(x, y) sf::st_coordinates(y[x,])))
        }
      }
    }
    toSample <- which(is.na(site.sample))
  }
  if (!is.null(maxTries) && draws > maxTries)
    stop(paste0('maxTries (',maxTries,') reached during point resampling, exiting. Examine targetSites, geoBias, and geoVcov to determine why so few resampled points fall within targetSites.'))
#  }
  return(list(site.sample = site.sample, point.sample = point.sample,
              draws = draws))
}


# Samples from isotope target points
# Called by estMCisotope
targetSampleIsotope <- function(targetIntrinsic, animal.sample,
                         nSim = 10, targetSites = NULL,
                         resampleProjection = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs",
                         maxTries = 300, pointsAssigned = FALSE,
                         targCon = NULL, pointsInSites = FALSE) {
  # Here is a catch to save Mike from himself
  #if(is.null(nSim)){nSim<-1000}
  if (pointsAssigned) {
    nAnimals <- dim(targetIntrinsic$SingleCell)[3]
    # DETERMINE THE NUMBER OF RANDOM DRAWS FROM ISOSCAPE used in isoAssign
    nrandomDraws <- dim(targetIntrinsic$SingleCell)[1]

    # Stack 3D array into 2D, to make it easier to get different random point samples for each animal
    targetIntrinsic1 <- array(aperm(targetIntrinsic$SingleCell, c(1, 3, 2)), dim = c(nAnimals * nrandomDraws, 2))

    # project target points to WGS #
    targetIntrinsic2 <- sp::SpatialPoints(targetIntrinsic1, proj4string = sp::CRS(MigConnectivity::projections$WGS84))
  }
  else {
    nAnimals <- dim(targetIntrinsic$probassign)[3]
  }
  if (!pointsInSites)
    # converts raster to matrix of XY then probs
    matvals <- raster::rasterToPoints(targetIntrinsic$probassign)
  target.sample <- rep(NA, nAnimals)
  target.point.sample <- matrix(NA, nAnimals, 2)
  toSample <- 1:nAnimals
  # In this case, only need to draw once, because already determined that points fall in targetSites
  if (pointsInSites) {
    draws <- 1
    samp <- sample.int(nrandomDraws, size = length(toSample), replace = T)
    samp2 <- samp + (animal.sample[toSample] - 1) * nrandomDraws
    point.sample <- targetIntrinsic2[samp2]
    # Changed to make sure x,y coords stack correctly
    target.point.sample[toSample,1]<- point.sample@coords[,1]
    target.point.sample[toSample,2]<- point.sample@coords[,2]
    if (!is.null(targCon))
      # Grab the relevant target sites
      target.sample[toSample] <- targCon[samp2]
  }
  # In this case, only need to draw once, as no limits on where points can be
  else if (is.null(targetSites)) {
    draws <- 1
    # This draws samples nSamples per animal (faster than looping over nSamples) and fills the xysim with x,y coords
    for(i in 3:ncol(matvals)) {
      animals <- which(animal.sample==i - 2)
      if (length(animals) > 0) {
        multidraw <- rmultinom(n = length(animals), size = 1, prob = matvals[,i])
        target.point.sample[animals,1] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],1]
        target.point.sample[animals,2] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],2]
      }
    }
  }
  else {
    draws <- 0
    # Make sure targetSites are WGS84
    targetSites <- sp::spTransform(targetSites, sp::CRS(MigConnectivity::projections$WGS84))
    while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
      draws <- draws + 1
      if (!pointsAssigned) {
        # This draws samples nSamples per animal (faster than looping over nSamples) and fills the xysim with x,y coords
        for(i in 3:ncol(matvals)) {
          animals <- which(animal.sample[toSample]==i - 2)
          if (length(animals) > 0) {
            multidraw <- rmultinom(n = length(animals)*nSim, size = 1, prob = matvals[,i])
            point.sample <- matvals[which(multidraw == 1, arr.ind = TRUE)[, 1], 1:2]
            point.sample1 <- sp::SpatialPoints(point.sample, proj4string = sp::CRS(MigConnectivity::projections$WGS84))
            # Check which points are in target sites
            target.sample0 <- sp::over(point.sample1, y = targetSites)
            good.sample <- which(!is.na(target.sample0))
            nToFill <- min(length(good.sample), length(animals))
            if (nToFill > 0) {
              target.sample[animals[1:nToFill]] <- target.sample0[good.sample[1:nToFill]]
              target.point.sample[animals[1:nToFill], 1] <- point.sample[good.sample[1:nToFill], 1]
              target.point.sample[animals[1:nToFill], 2] <- point.sample[good.sample[1:nToFill], 2]
            }
          }
        }
      }
      else {
        # Select nSim points for each animal still to be sampled
        samp <- sample.int(nrandomDraws, size = length(toSample) * nSim, replace = T)
        samp2 <- samp + rep(animal.sample[toSample] - 1, each = nSim) * nrandomDraws
        point.sample <- targetIntrinsic2[samp2]
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
#' @param pos Number of sites by 2 matrix with positions of each site.  If
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
        dist[b1, b2] <- geosphere::distGeo(pos[b1, ], pos[b2, ]) / 1000
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
