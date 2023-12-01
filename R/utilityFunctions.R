# Function for generating random transition probability matrix (expression of parametric uncertainty)
makePsiRand <- function(model, origin.set, target.set) {
  vcv1 <- model$results$beta.vcv
  beta1 <- model$results$beta$estimate
  new.beta <- MASS::mvrnorm(1, beta1, vcv1)
  full.beta <- rep(NA, length(model$results$beta$estimate))
  new.real <- RMark::get.real(model, "Psi", new.beta,
                              design = model$design.matrix, vcv=FALSE, se=FALSE)
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
      if(!sf::st_crs(targetSites)==sf::st_crs(targetPoints)){targetPoints <- sf::st_transform(targetPoints, crs = resampleProjection)}
      targetAssignment <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = targetPoints, y = targetSites, sparse = TRUE))))
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
    draws <- 0
    while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
      draws <- draws + 1
      geoBias2 <- array(rep(geoBias, length(toSample), each = nSim), c(nSim, 2, length(toSample)))

      # generate sample and substract the bias
      point.sample <- array(apply(sf::st_coordinates(targetPoints)[animal.sample[toSample], , drop = FALSE],
                                  MARGIN = 1,
                                  MASS::mvrnorm, n=nSim, Sigma=geoVCov),
                            # dimensions of the array
                            c(nSim, 2, length(toSample))) - geoBias2 #subtract the bias

      # name so can
      dimnames(point.sample)[[2]] <- c("x","y")

      point.sample <- apply(point.sample, FUN = function(x){
        sf::st_as_sf(data.frame(x), coords = c("x","y"),
                     crs = resampleProjection)}, MARGIN = 3)
      good.sample <- apply(sapply(point.sample, sf::st_is_valid), 2,
                           function(x) min(which(x)))
      if (any(!is.na(good.sample)))
        target.point.sample[toSample[!is.na(good.sample)], ]<- t(mapply(x = good.sample[!is.na(good.sample)],
                                                                        y = point.sample[!is.na(good.sample)],
                                                                        FUN = function(x, y){sf::st_coordinates(y[x,])}))

      toSample <- which(is.na(target.point.sample[,1]))
    }
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

      # Convert those to sf

      point.sample <- apply(point.sample, FUN = function(x){
        sf::st_as_sf(data.frame(x), coords = c("x","y"),
                     crs = resampleProjection)}, MARGIN = 3)
      #point.sample <- lapply(point.sample1, st_as_sf)
      # Find out which sampled points are in a target site
      if(!sf::st_crs(targetSites)==sf::st_crs(point.sample)){
        targetSites <- sf::st_transform(targetSites, crs = resampleProjection)
      }
      target.sample0 <- sapply(point.sample, FUN = function(z){
        inter <- sf::st_intersects(x = z, y = targetSites,
                                   sparse = TRUE)
        inter[lengths(inter)==0] <- NA
        inter[lengths(inter)>1] <- sapply(inter[lengths(inter)>1],
                                          function(x) x[1])
        as.numeric(unlist(unclass(inter)))})
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
                                                                        FUN = function(x, y){sf::st_coordinates(y[x,])}))

      toSample <- which(is.na(target.sample))
    }
    if (!is.null(maxTries) && draws > maxTries)
      stop(paste0('maxTries (',maxTries,') reached during point resampling, exiting. Examine targetSites, geoBias, and geoVcov to determine why so few resampled points fall within targetSites.'))
  }
  return(list(target.sample = target.sample, target.point.sample = target.point.sample,
              draws = draws))
}

# Function for sampling from geolocator, telemetry, and/or intrinsic location
# uncertainty on either origin or target side
locSample <- function(isGL,
                      isRaster,
                      isProb,
                      isTelemetry,
                      geoBias = NULL,
                      geoVCov = NULL,
                      points = NULL,
                      matvals = NULL,
                      matvals_crs = NULL,
                      singleCell = NULL,
                      pointsInSites = FALSE,
                      overlap1 = NULL,
                      assignment = NULL,
                      sites = NULL,
                      resampleProjection = 'ESRI:53027',
                      nSim = 1000,
                      maxTries = 300) {
  #cat("Starting locSample\n")
  # points should come in as sf object #
  nAnimals <- length(isGL)
  if (!is.null(singleCell)) {
    # DETERMINE THE NUMBER OF RANDOM DRAWS FROM ISOSCAPE used in isoAssign

    # How does targetIntrinsic make it into this function?
    nrandomDraws <- dim(singleCell)[1]

    # Stack 3D array into 2D, to make it easier to get different random point samples for each animal
    intrinsic1 <- array(aperm(singleCell, c(1, 3, 2)),
                        dim = c(nAnimals * nrandomDraws, 2),
                        dimnames = list(NULL, c("x", "y")))

    # project points to WGS #
    intrinsic2 <- sf::st_as_sf(as.data.frame(intrinsic1), coords = c("x", "y"),
                               crs = 4326)
  }

  # If assignment isn't defined, create one from sites and points
  if (is.null(assignment)) {
    if (!is.null(sites) && !is.null(points)){
      #Use the provided spatial layers to create assignment
      if(!sf::st_crs(sites)==sf::st_crs(points)){
        points <- sf::st_transform(points, crs = resampleProjection)
        sites <- sf::st_transform(sites, crs = resampleProjection)
      }
      assignment <- suppressMessages(unclass(sf::st_intersects(x = points,
                                                               y = sites,
                                                    sparse = TRUE)))
      #if assignment is not there - convert to 0
      assignment[sapply(assignment, function(x) length(x)==0)] <- 0
      assignment <- array(unlist(assignment), nAnimals)
    }
    else
      assignment <- array(0, nAnimals)
  }

  # Storage
  site.sample <- rep(NA, nAnimals)
  point.sample <- matrix(NA, nAnimals, 2, dimnames = list(NULL, c('x', 'y')))
  # Fill in telemetry/GPS values (no location uncertainty)
  if (length(dim(assignment))==2){
    site.sample[which(isTelemetry)] <- apply(assignment[isTelemetry, ,
                                                        drop = FALSE], 1,
                                                  which.max)
    if (is.list(site.sample)) {
      site.sample[lengths(site.sample)==0] <- 0
      site.sample <- unlist(site.sample)
    }

  }else{
    site.sample[which(isTelemetry)] <- assignment[which(isTelemetry)]
    site.sample[which(isTelemetry & is.na(site.sample))] <- 0
  }
  if (!is.null(points))
    point.sample[which(isTelemetry), ] <- sf::st_coordinates(points[which(isTelemetry),])
  # Which to sample still
  toSample <- which(!isTelemetry)
  draws <- 0
  # Loop until all animals have valid sample or have exceeded maxTries
  while (length(toSample) > 0 && (is.null(maxTries) || draws < maxTries)) {
    draws <- draws + 1
    # Convert toSample (numbers) into T/F, so can combine with other T/F
    toSampleBool <- 1:nAnimals %in% toSample
    # cat(toSampleBool, '\n')
    # cat(isGL[1:10], '\n')
    # Sample geolocator points
    if (any(isGL & toSampleBool)) {
      #cat("*************", sum(isGL & toSampleBool), "*************\n")
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
      if (!is.null(sites)) {
        # Find out which sampled points are in a target site
        if(!sf::st_crs(sites)==sf::st_crs(point.sample0[[1]])){
          sites <- sf::st_transform(sites, crs = resampleProjection)
        }
        site.sample0 <- sapply(point.sample0, FUN = function(z){
          inter <- sf::st_intersects(x = z, y = sites,
                                     sparse = TRUE)
          inter[lengths(inter)==0] <- NA
          inter[lengths(inter)>1] <- sapply(inter[lengths(inter)>1], function(x) x[1])
          as.numeric(unlist(unclass(inter)))})

        # Identify which animals have at least one valid sample point. good.sample0
        # will be NA for those that don't.  For those that do, it will location in
        # point.sample first valid point can be found.
        good.sample0 <- apply(site.sample0, 2, function(x) which(!is.na(x))[1])
      }
      else {
        site.sample0 <- array(0, c(nSim, length(point.sample0)))
        good.sample0 <- rep(1, length(point.sample0))
      }
      # cat(good.sample0, "\n")
    }
    if (any(isProb & toSampleBool)) {
      site.sample1 <- apply(assignment[isProb & toSampleBool, , drop = FALSE],
                            1, sample.int, n = ncol(assignment), size = nSim,
                            replace = TRUE)
      if (is.null(dim(site.sample1)))
        site.sample1 <- matrix(site.sample1, 1)
    }
    # walk through this section #
    if (any(isRaster & toSampleBool)) {
      # some raster -
      if (pointsInSites && !any(isRaster & toSampleBool & (isGL | isProb))) {
        # IF ALL THE POINTS ARE WITHIN SITES and no animals have two data types
        samp <- sample.int(nrandomDraws,
                           size = sum(isRaster & toSampleBool),
                           replace = TRUE)
      # GRAB LOCATIONS FROM RASTER #
        samp2 <- samp + ((1:nAnimals)[toSampleBool & isRaster] - 1) * nrandomDraws

        point.sample2 <- intrinsic2[samp2, ]

        # Changed to make sure x,y coords stack correctly
        # target.point.sample[toSample,1]<- point.sample@coords[,1]
        # target.point.sample[toSample,2]<- point.sample@coords[,2]

# ANIMALS WITHIN samp2 that don't have raster
# We ULTIMATELY NEED POINTS BECAUSE THEY ONLY HAVE RASTER DATA
        if (!is.null(overlap1))
          # Grab the relevant sites
          site.sample2 <- overlap1[samp, (1:nAnimals)[toSampleBool & isRaster], drop = FALSE]
        good.sample2 <- seq(from = 1, by = nSim,
              length.out = sum(isRaster & toSampleBool))
      }
      else if (is.null(singleCell)) {
        multidraw <- apply(matvals[ , (1:nAnimals)[toSampleBool & isRaster] + 2, drop = FALSE],
                           2, function(x) (rmultinom(prob = x, n = nSim, size = 1)))
        multidraw <- array(multidraw, c(nrow(matvals), nSim,
                                        sum(toSampleBool & isRaster)))
        #point.sample2 <- matvals[which(multidraw == 1, arr.ind = TRUE)[, 1], 1:2]
        point.sample2 <- apply(multidraw, c(2:3), function(x) which(x==1))
        point.sample2 <- array(matvals[point.sample2, 1:2],
                               c(nSim * sum(toSampleBool & isRaster), 2),
                               list(NULL, c("x","y")))
        #dimnames(point.sample2)[[2]] <- c("x","y")

        # Convert those to SpatialPoints
        #cat("------ I MADE IT TO POINTSAMPLE 2 --------\n")
        point.sample2 <- sf::st_as_sf(as.data.frame(point.sample2),
                                      coords = c("x","y"),
                                      #crs = 4326)
                                      crs = sf::st_crs(matvals_crs))

        # point.sample2 <- apply(point.sample2,
        #                        FUN = function(x){
        #                          sf::st_as_sf(data.frame(x),
        #                                       coords = c("x","y"),
        #                             crs = 4326)},
        #                        MARGIN = 2)
        #point.sample <- lapply(point.sample1, st_as_sf)
        if (!is.null(sites)) {
          # Find out which sampled points are in a target site
          if(!sf::st_crs(sites)==sf::st_crs(point.sample2)){
            sites <- sf::st_transform(sites, crs = resampleProjection)
            point.sample2 <- sf::st_transform(point.sample2,
                                              crs = resampleProjection)
          }

          # site.sample2 <- sapply(point.sample2, FUN = function(z){
          #   suppressMessages(as.numeric(unclass(sf::st_intersects(x = z, y = sites,
          #                                        sparse = TRUE))))})
          # Check which points are in target sites
          site.sample2 <- suppressMessages(sf::st_intersects(point.sample2,
                                                             y = sites))

          # Give values without intersection an NA
          len_intersect <- lengths(site.sample2)
          # quick check to ensure that all the points fall exactly in one site #
          if(any(len_intersect)>1){
            warning("Overlapping targetSites or originSites may cause issues\n")
            site.sample2 <- lapply(site.sample2, function (x) x[1])
          }
          site.sample2[lengths(site.sample2)==0] <- NA
          site.sample2 <- unlist(as.numeric(site.sample2))

          # Organize into matrix (separate by animal)
          site.sample2 <- matrix(site.sample2, nSim, sum(isRaster & toSampleBool))

          # Identify which animals have at least one valid sample point. good.sample2
          # will be NA for those that don't.  For those that do, it will location in
          # point.sample first valid point can be found.
          good.sample2 <- apply(site.sample2, 2, function(x) which(!is.na(x))[1])+
            seq(from = 0, by = nSim, length.out = sum(isRaster & toSampleBool))
          #cat(good.sample2, "\n")
        }
        else {
          site.sample2 <- array(0, c(nSim, sum(isRaster & toSampleBool)))
          good.sample2 <- rep(1, sum(isRaster & toSampleBool))
        }
      }
      else {
        # Select nSim points for each animal still to be sampled
        samp <- sample.int(nrandomDraws,
                           size = sum(isRaster & toSampleBool) * nSim,
                           replace = TRUE)

        samp2 <- samp + rep((1:nAnimals)[toSampleBool & isRaster] - 1, each = nSim) * nrandomDraws

        point.sample2 <- intrinsic2[samp2, ]
        if (!is.null(sites)) {
          if(!sf::st_crs(sites)==sf::st_crs(point.sample2)){
            point.sample2 <- sf::st_transform(point.sample2,
                                              crs = resampleProjection)
          }
          # Check which points are in target sites
          site.sample2 <- suppressMessages(sf::st_intersects(point.sample2,
                                                             y = sites))

          # Give values without intersection an NA
          len_intersect <- lengths(site.sample2)
          # quick check to ensure that all the points fall exactly in one site #
          if(any(len_intersect)>1){
            warning("Overlapping targetSites or originSites may cause issues\n")
            site.sample2 <- lapply(site.sample2, function (x) x[1])
          }
          site.sample2[len_intersect==0] <- NA
          site.sample2 <- unlist(as.numeric(site.sample2))

          # Organize into matrix (separate by animal)
          site.sample2 <- matrix(site.sample2, nSim, sum(isRaster & toSampleBool))
          # Identify which animals have a valid point (inside a site).
          # good.sample2, for those that don't, will be NA. For those that do, it
          # will be location in point.sample where first valid point is.


          good.sample2 <- apply(site.sample2, 2, function(x){which(!is.na(x))[1]})+
            seq(from = 0, by = nSim, length.out = sum(isRaster & toSampleBool))
        }
        else {
          site.sample2 <- array(0, c(nSim, sum(isRaster & toSampleBool)))
          good.sample2 <- rep(1, sum(isRaster & toSampleBool))
        }
      }
    }
    if (any(isProb & !isGL & !isRaster & toSampleBool)){
      site.sample[isProb & !isGL & !isRaster & toSampleBool] <-
        site.sample1[1, !(which(isProb & toSampleBool) %in% which(isGL | isRaster)), drop = FALSE]
    }
    if (any(!isProb & isGL & !isRaster & toSampleBool)){
      if (any(!is.na(good.sample0[!(which(isGL) %in% which(isProb | isRaster))]))){
        site.sample[!isProb & isGL & !isRaster & toSampleBool] <-
          apply(site.sample0[, !(which(isGL & toSampleBool) %in% which(isProb | isRaster)), drop = FALSE],
                2,
                function(x) x[!is.na(x)][1])
        # Fill in target points of valid sampled points
        point.sample[which(!isProb & isGL & !isRaster & toSampleBool)[which(!is.na(good.sample0[!(which(isGL & toSampleBool) %in% which(isProb | isRaster))]))], ]<-
          t(mapply(x = good.sample0[!is.na(good.sample0) & !(which(isGL & toSampleBool) %in% which(isProb | isRaster))],
                   y = point.sample0[which(!(which(isGL & toSampleBool) %in% which(isProb | isRaster)) & (!is.na(good.sample0)))],
                   FUN = function(x, y) sf::st_coordinates(y[x,])))
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
          # point.sample[which(isProb & isGL & !isRaster & toSampleBool)[which(!is.na(good.sample01))], ]<-
          #   t(mapply(x = good.sample01[!is.na(good.sample01)],
          #            y = point.sample0[which(which(isGL & toSampleBool) %in% which(isProb & !isRaster))[which(!is.na(good.sample01))]],
          #            FUN = function(x, y) sf::st_coordinates(y[x,])))
        }
      }
    }
    if (any(!isProb & !isGL & isRaster & toSampleBool)){
     if (any(!is.na(good.sample2[!(which(isRaster) %in% which(isProb | isGL))]))){
        site.sample[which(!isProb & !isGL & isRaster & toSampleBool)] <-
          apply(site.sample2[, !(which(isRaster & toSampleBool) %in% which(isProb | isGL)),
                             drop = FALSE],
                2,
                function(x) x[!is.na(x)][1])
        rows.select <- good.sample2[!(which(isRaster & toSampleBool) %in%
                                        which(isProb | isGL) |
                                        is.na(good.sample2))]
        # Fill in target points of valid sampled points
        point.sample[which(!isProb & !isGL & isRaster & toSampleBool)[which(!is.na(good.sample2))], ] <-
          sf::st_coordinates(point.sample2[rows.select, ])
          #[which(!is.na(good.sample2))]
          # t(mapply(x = good.sample2[!(which(isRaster & toSampleBool) %in% which(isProb | isGL))],
          #          y = point.sample2[which(which(isRaster & toSampleBool) %in% which(!isProb & !isGL & !is.na(good.sample2))),,
          #                            drop = FALSE],
          #          FUN = function(x, y) sf::st_coordinates(y[x,])))
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
          # Doesn't get used, not going to include prob with estMantel
          # point.sample[which(isProb & !isGL & isRaster & toSampleBool)[which(!is.na(good.sample12))], ]<-
          #   t(mapply(x = good.sample12[!is.na(good.sample12)],
          #            y = point.sample2[which(which(isRaster & toSampleBool) %in% which(isProb & !isGL))[which(!is.na(good.sample12))], ],
          #            FUN = function(x, y) sf::st_coordinates(y[x,])))
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
        #cat(good.sample02, '\n')
        if (any(!is.na(good.sample02))){
          for (col in 1:sum(!isProb & isGL & isRaster & toSampleBool))
            site.sample[which(!isProb & isGL & isRaster & toSampleBool)[col]] <-
              site.sample2[good.sample02[col],
                           which(which(isRaster & toSampleBool) %in% which(!isProb & isGL))[col]]
          # Just choosing GL point instead of fussing with raster point - generally
          # those are a lot more precise so probably also more accurate.
          point.sample[which(!isProb &
                               isGL &
                               isRaster &
                               toSampleBool)[which(!is.na(good.sample02))],] <-
            t(mapply(
              x = good.sample02[!is.na(good.sample02)],
              y = point.sample0[which(which(isGL & toSampleBool) %in%
                                        which(!isProb & isRaster))[which(!is.na(good.sample02))]],
              FUN = function(x, y)
                sf::st_coordinates(y[x, ])
            ))
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
          # This one might not matter (not planning to include isProb in
          # estMantel or look at points in estTransition).
          # point.sample[which(isProb & isGL & isRaster & toSampleBool)[which(!is.na(good.sample012))]]<-
          #   t(mapply(x = good.sample012[!is.na(good.sample012)],
          #            y = point.sample0[which(which(isGL & toSampleBool) %in% which(isProb & isRaster))[which(!is.na(good.sample012))]],
          #            FUN = function(x, y) sf::st_coordinates(y[x,])))
        }
      }
    }
    toSample <- which(is.na(site.sample))
  }
  if (length(toSample) > 0){
    notfind <- data.frame(Animal = toSample, isGL = isGL[toSample],
                          isRaster = isRaster[toSample], isProb = isProb[toSample],
                          isTelemetry = isTelemetry[toSample])
    return(list(site.sample = site.sample, point.sample = point.sample,
                draws = draws, notfind = notfind))
  }
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
    targetIntrinsic1 <- array(aperm(targetIntrinsic$SingleCell, c(1, 3, 2)),
                              dim = c(nAnimals * nrandomDraws, 2),
                              dimnames = list(NULL, c("x","y")))

    # project target points to WGS #
    targetIntrinsic2 <- sf::st_as_sf(as.data.frame(targetIntrinsic1), coords = c("x", "y"),
                                     crs = 4326)

  }
  else {
    nAnimals <- dim(targetIntrinsic$probassign)[3]
  }
  if (!pointsInSites){
    if (!inherits(targetIntrinsic$probassign, "SpatRaster"))
      targetIntrinsic$probassign <- terra::rast(targetIntrinsic$probassign)
    # converts raster to matrix of XY then probs
    matvals <- terra::as.data.frame(targetIntrinsic$probassign, xy = TRUE)
  }
  target.sample <- rep(NA, nAnimals)
  target.point.sample <- matrix(NA, nAnimals, 2)
  toSample <- 1:nAnimals
  # In this case, only need to draw once, because already determined that points fall in targetSites
  if (pointsInSites) {
    draws <- 1
    samp <- sample.int(nrandomDraws, size = length(toSample), replace = TRUE)
    samp2 <- samp + (animal.sample[toSample] - 1) * nrandomDraws
    point.sample <- targetIntrinsic2[samp2, ]
    # Changed to make sure x,y coords stack correctly
    #target.point.sample[toSample,1]<- point.sample@coords[,1]
    #target.point.sample[toSample,2]<- point.sample@coords[,2]
    target.point.sample[toSample,1] <- sf::st_coordinates(point.sample)[,1]
    target.point.sample[toSample,2] <- sf::st_coordinates(point.sample)[,2]
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
    targetSites <- sf::st_transform(targetSites, 4326)

    while (length(toSample) > 0 && (is.null(maxTries) || draws <= maxTries)) {
      draws <- draws + 1
      if (!pointsAssigned) {
        # This draws samples nSamples per animal (faster than looping over nSamples) and fills the xysim with x,y coords
        for(i in 3:ncol(matvals)) {
          animals <- which(animal.sample[toSample]==i - 2)
          if (length(animals) > 0) {
            multidraw <- rmultinom(n = length(animals)*nSim, size = 1, prob = matvals[,i])

            point.sample <- matvals[which(multidraw == 1, arr.ind = TRUE)[, 1], 1:2]

            # covert to sf
            point.sample1 <- sf::st_as_sf(as.data.frame(point.sample),
                                          coords = c("x","y"),
                                          crs = 4326)

            # Check which points are in target sites

            target.sample0 <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = point.sample1,
                                                                                    y = targetSites,
                                                                                    sparse = TRUE))))

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
        samp <- sample.int(nrandomDraws, size = length(toSample) * nSim,
                           replace = TRUE)
        samp2 <- samp + rep(animal.sample[toSample] - 1, each = nSim) * nrandomDraws
        point.sample <- targetIntrinsic2[samp2,]
        # Check which points are in target sites
        target.sample0 <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = point.sample,
                                             y = targetSites,
                                             sparse = TRUE))))
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
          target.point.sample[toSample[!is.na(good.sample)],1]<- sf::st_coordinates(point.sample[good.sample[!is.na(good.sample)],])[,1]
          target.point.sample[toSample[!is.na(good.sample)],2]<- sf::st_coordinates(point.sample[good.sample[!is.na(good.sample)],])[,2]
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
#'    longitude and column 2 should be latitude. If \code{surface} is 'plane',
#'    column 1 can be x-position and column 2 y-position.
#' @param surface Surface to calculate distances on.  Either 'ellipsoid'
#'    (default), 'sphere', or 'plane'.
#' @param units Units of return distance matrix. If \code{surface} is 'plane',
#'    then this argument is ignored and the return units will be the same as the
#'    \code{pos} units. Options are 'km' (kilometers, default), 'm' (meters),
#'    'miles', and 'nautical miles'.
#'
#' @return Square matrix of distances between sites. If \code{surface} is
#'    'ellipsoid' or 'sphere', then argument \code{units} will determine units;
#'    if \code{surface} is 'plane', the units will be the same as the \code{pos}
#'    units.
#'
#' @note In version 0.4.3 we switched package dependencies from \code{geosphere}
#'    to \code{geodist}. As a result, spherical distances (and possibly
#'    ellipsoid distances) may differ slightly from those calculated with earlier
#'    versions of our package.
#'
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
distFromPos <- function(pos, surface = 'ellipsoid',
                        units = c("km", "m", "miles", "nautical miles")) {
  if (!is.matrix(pos) || dim(pos)[2] != 2)
    stop("pos should be a matrix with two column (e.g., longitude and latitude)")
  units <- match.arg(units)
  if (!(surface %in% c('plane', 'sphere', 'ellipsoid')))
    stop('surface must be "plane", "sphere", or "ellipsoid".')
  if (is.null(colnames(pos)) && surface != 'plane')
    colnames(pos) <- c("lon", "lat")
  div <- ifelse(units == "km", 1000,
                 ifelse(units == "m", 1,
                        ifelse(units == "miles", 1609.34, 1852)))
  if (surface == 'ellipsoid')
    dist <- geodist::geodist(pos, measure = "geodesic") / div
  else if (surface == 'sphere')
    dist <- geodist::geodist(pos, measure = "vincenty") / div
  else {
    nSites <- nrow(pos)
    dist <- matrix(0, nSites, nSites)
    for (b1 in 2:nSites) {
      for (b2 in 1:(b1-1)) {
        dist[b1, b2] <- sqrt((pos[b1, 1] - pos[b2, 1]) ^ 2 +
                              (pos[b1, 2] - pos[b2, 2]) ^ 2)
        dist[b2, b1] <- dist[b1, b2]
      }
    }
  }
  return(dist)
}

reassignInds <- function(dataOverlapSetting = "none",
                         originPoints = NULL, targetPoints = NULL,
                         originAssignment = NULL, targetAssignment = NULL,
                         isGL = FALSE, isTelemetry = FALSE,
                         isRaster = FALSE, isProb = FALSE,
                         captured = "origin",
                         #originRaster = NULL,
                         originRasterXYZ = NULL,
                         originRasterXYZcrs = NULL,
                         originSingleCell = NULL,
                         #targetRaster = NULL,
                         targetRasterXYZ = NULL,
                         targetRasterXYZcrs = NULL,
                         targetSingleCell = NULL,
                         originSites = NULL,
                         targetSites = NULL) {
  if (dataOverlapSetting != "none") {
    stop('dataOverlapSetting "named" not set up yet')
  }
  if (is.null(targetPoints))
    nTargetPoints <- 0
  else
    nTargetPoints <- nrow(targetPoints)
  if (is.null(originPoints))
    nOriginPoints <- 0
  else
    nOriginPoints <- nrow(originPoints)
  if (any(isGL) || any(isTelemetry)) {
    nPoints <- max(nTargetPoints, nOriginPoints)
  }
  else {
    nPoints <- 0
  }
  # capturePointsGL <- T
  # if (any(isGL))
  #   if (sum(captured[isGL]=="origin")>nOriginPoints)
  if (any(isRaster)) {
    if (!is.null(originRasterXYZ)) {
      dimORast <- dim(originRasterXYZ)
      nORast <- dimORast[2] - 2
    }
    else
      nORast <- 0
    if (!is.null(targetRasterXYZ)) {
      dimTRast <- dim(targetRasterXYZ)
      nTRast <- dimTRast[2] - 2
    }
    else
      nTRast <- 0
  }
  else {
    nORast <- nTRast <- 0
  }
  if (any(isProb)) {
    if (is.null(originAssignment)){
      nOAss <- 0
      dimOAss <- NULL
    }
    else {
      dimOAss <- dim(originAssignment)
      if (is.null(dimOAss)) {
        nOAss <- length(originAssignment)
        tt2 <- matrix(0, nOAss, max(length(unique(originAssignment)),
                                    nrow(originSites)))
        for (i in 1:nOAss)
          tt2[i, originAssignment[i]] <- 1
        originAssignment <- tt2
        dimOAss <- dim(originAssignment)
      }
      else
        nOAss <- dimOAss[1]
    }
    if (is.null(targetAssignment)){
      nTAss <- 0
      dimTAss <- NULL
    }
    else {
      dimTAss <- dim(targetAssignment)
      if (is.null(dimTAss)){
        nTAss <- length(targetAssignment)
        tt2 <- matrix(0, nTAss, max(length(unique(targetAssignment)),
                                    nrow(targetSites)))
        for (i in 1:nTAss)
          tt2[i, targetAssignment[i]] <- 1
        targetAssignment <- tt2
        dimTAss <- dim(targetAssignment)
      }
      else
        nTAss <- dimTAss[1]
    }
  }
  nTotal <- max(length(isGL), length(isTelemetry), length(isRaster),
                length(isProb))
  if (nTotal < 2) {
    nTotal <- length(captured)
    if (nTotal < 2 &&
        ((any(isGL) + any(isTelemetry) + any(isRaster) + any(isProb)) > 1)) {
      warning("For multiple datasets, we recommend having an entry in captured, isGL, isTelemetry, etc. for each animal")
      nTotal <- nPoints + nORast + nTRast + max(nOAss, nTAss)
    }
  }
  if (length(isGL)==1)
    isGL <- rep(isGL, nTotal)
  if (length(isTelemetry)==1)
    isTelemetry <- rep(isTelemetry, nTotal)
  if (length(isRaster)==1)
    isRaster <- rep(isRaster, nTotal)
  if (length(isProb)==1)
    isProb <- rep(isProb, nTotal)
  if (length(captured)==1)
    captured <- rep(captured, nTotal)
  capturePoint <- rep(FALSE, nTotal)
  nTPremain <- nTargetPoints
  nOPremain <- nOriginPoints
  if (any(isGL)){
    nTPremain <- nTPremain - sum(captured[isGL]!="target")
    nOPremain <- nOPremain - sum(captured[isGL]!="origin")
    if (nTPremain > 0) {
      capturePoint[captured=="target" & isGL] <- TRUE
      nTPremain <- nTPremain - sum(captured[isGL]=="target")
    }
    if (nOPremain > 0) {
      capturePoint[captured=="origin" & isGL] <- TRUE
      nOPremain <- nOPremain - sum(captured[isGL]=="origin")
    }
  }
  if ((nTPremain > 0 || nOPremain > 0) && any(isTelemetry)) {
    nTPremain <- nTPremain - sum(captured[isTelemetry]!="target")
    nOPremain <- nOPremain - sum(captured[isTelemetry]!="origin")
    if (nTPremain > 0) {
      capturePoint[captured=="target" & isTelemetry] <- TRUE
      nTPremain <- nTPremain - sum(captured[isTelemetry]=="target")
    }
    if (nOPremain > 0) {
      capturePoint[captured=="origin" & isTelemetry] <- TRUE
      nOPremain <- nOPremain - sum(captured[isTelemetry]=="origin")
    }
  }
  if ((nTPremain > 0 || nOPremain > 0) && any(isRaster)) {
    if (nTPremain > 0) {
      capturePoint[captured=="target" & isRaster] <- TRUE
      nTPremain <- nTPremain - sum(captured[isRaster]=="target")
    }
    if (nOPremain > 0) {
      capturePoint[captured=="origin" & isRaster] <- TRUE
      nOPremain <- nOPremain - sum(captured[isRaster]=="origin")
    }
  }
  if ((nTPremain > 0 || nOPremain > 0) && any(isProb)) {
    if (nTPremain > 0) {
      capturePoint[captured=="target" & isProb] <- TRUE
      nTPremain <- nTPremain - sum(captured[isProb]=="target")
    }
    if (nOPremain > 0) {
      capturePoint[captured=="origin" & isProb] <- TRUE
      nOPremain <- nOPremain - sum(captured[isProb]=="origin")
    }
  }
  if (nTargetPoints > 0 || nOriginPoints > 0) {
    if (is.null(targetPoints)) {
      targetPoints2 <- NULL
    }
    else {
      if (nTargetPoints < nTotal) {
        dummyPoint <- targetPoints[1, ]
        targetPoints2 <- NULL
        place <- 0
        for (i in 1:nTotal) {
          if (isGL[i] || isTelemetry[i] || (captured[i] == "target" &&
                                            capturePoint[i])) {
            place <- place + 1
            targetPoints2 <- rbind(targetPoints2,
                                   targetPoints[place, ])
          }
          else
            targetPoints2 <- rbind(targetPoints2,
                                   dummyPoint)
        }
      }
      else
        targetPoints2 <- targetPoints
    }
    if (is.null(originPoints)) {
      originPoints2 <- NULL
    }
    else {
      if (nOriginPoints < nTotal) {
        dummyPoint <- originPoints[1, ]
        originPoints2 <- NULL
        place <- 0
        for (i in 1:nTotal) {
          if (isGL[i] || isTelemetry[i] || captured[i] == "origin" &&
              capturePoint[i]) {
            place <- place + 1
            originPoints2 <- rbind(originPoints2,
                                   originPoints[place, ])
          }
          else
            originPoints2 <- rbind(originPoints2,
                                   dummyPoint)
        }
      }
      else
        originPoints2 <- originPoints
    }
  }
  else {
    originPoints2 <- originPoints; targetPoints2 <- targetPoints
  }
  if (any(isRaster)) {
    if (!is.null(originRasterXYZ)) {
      if (all(isRaster & captured != "origin")) {
        originRasterXYZ2 <- originRasterXYZ
        originRaster2 <- originRaster
      }
      else {
        originRasterXYZ2 <- originRasterXYZ[, 1:2]
        column <- 2
        for (i in 1:nTotal) {
          if (isRaster[i] && captured[i] != "origin") {
            column <- column + 1
            originRasterXYZ2 <- cbind(originRasterXYZ2,
                                      originRasterXYZ[, column])
          }
          else {
            originRasterXYZ2 <- cbind(originRasterXYZ2,
                                      array(1/dimORast[1], c(dimORast[1], 1)))
          }
        }
        if (!is.null(originSingleCell)) {
          dimOCell <- dim(originSingleCell)
          dummyVals <- c(originSingleCell[ , , 1])
          originSingleCell2 <- NULL
          place <- 0
          for (i in 1:nTotal) {
            if (isRaster[i] && captured[i] != "origin") {
              place <- place + 1
              originSingleCell2 <- c(originSingleCell2, originSingleCell[,,place])
            }
            else
              originSingleCell2 <- c(originSingleCell2, dummyVals)

          }
          originSingleCell <- array(originSingleCell2,
                                    c(dimOCell[1], dimOCell[2], nTotal),
                                    list(NULL, c("Longitude", "Latitude"), NULL))
        }
      }
    }
    else{
      originRasterXYZ2 <- NULL
    }
    if (!is.null(targetRasterXYZ)) {
      if (all(isRaster & captured != "target")) {
        targetRasterXYZ2 <- targetRasterXYZ
      }
      else {
        xvalues <- range(targetRasterXYZ[,1])
        yvalues <- range(targetRasterXYZ[,2])
        targetRasterXYZ2 <- targetRasterXYZ[, 1:2]
        column <- 2
        for (i in 1:nTotal) {
          if (isRaster[i] && captured[i] != "target") {
            column <- column + 1
            targetRasterXYZ2 <- cbind(targetRasterXYZ2,
                                      targetRasterXYZ[, column])
          }
          else{
            targetRasterXYZ2 <- cbind(targetRasterXYZ2,
                                      array(1/dimTRast[1], c(dimTRast[1], 1)))
          }
        }
        if (!is.null(targetSingleCell)) {
          dimTCell <- dim(targetSingleCell)
          dummyVals <- c(targetSingleCell[ , , 1])
          targetSingleCell2 <- NULL
          place <- 0
          for (i in 1:nTotal) {
            if (isRaster[i] && captured[i] != "target") {
              place <- place + 1
              targetSingleCell2 <- c(targetSingleCell2, targetSingleCell[,,place])
            }
            else
              targetSingleCell2 <- c(targetSingleCell2, dummyVals)

          }
          targetSingleCell <- array(targetSingleCell2,
                                    c(dimTCell[1], dimTCell[2], nTotal))
        }
      }
    }
    else {
      targetRasterXYZ2 <- NULL
    }
  }
  else {
    originRasterXYZ2 <- originRasterXYZ; targetRasterXYZ2 <- targetRasterXYZ
  }
  if (any(isProb)) {
    if (nTAss == 0) {
      targetAssignment2 <- NULL
    }
    else {
      if (nTAss < nTotal) {
        dummyAss <- array(1/dimTAss[2], c(1, dimTAss[2]))
        targetAssignment2 <- NULL
        place <- 0
        for (i in 1:nTotal) {
          if (isProb[i] && (captured[i] != "target" || !capturePoint[i])) {
            place <- place + 1
            targetAssignment2 <- rbind(targetAssignment2,
                                       targetAssignment[place, ])
          }
          else {
            if (captured[i] == "target" && capturePoint[i]){
              ass <- sf::st_intersects(x = targetPoints[i,], y = targetSites,
                                       sparse = FALSE)
              if (sum(ass)<1){
                warning("Not all target capture locations are within targetSites. Assigning to closest site.\n",
                        "Affects animal: ", i)
                ass <- matrix(0, 1, dimTAss[2])
                ass[sf::st_nearest_feature(x = targetPoints[i,],
                                           y = targetSites)] <- 1
              }
              else if (sum(ass) > 1){
                warning("Overlapping targetSites may cause issues\n")
                ass0 <- sf::st_intersects(x = targetPoints[i,], y = targetSites,
                                          sparse = TRUE)
                ass <- matrix(0, 1, dimTAss[2])
                ass[ass0[1]] <- 1
              }
              targetAssignment2 <- rbind(targetAssignment2, ass)
            }

            else
              targetAssignment2 <- rbind(targetAssignment2, dummyAss)
          }
        }
      }
      else
        targetAssignment2 <- targetAssignment
    }
    if (is.null(dimOAss)) {
      originAssignment2 <- NULL
    }
    else {
      if (dimOAss[1] < nTotal) {
        dummyAss <- array(1/dimOAss[2], c(1, dimOAss[2]))
        originAssignment2 <- NULL
        place <- 0
        for (i in 1:nTotal) {
          if (isProb[i] && (captured[i] != "origin" || !capturePoint[i])) {
            place <- place + 1
            originAssignment2 <- rbind(originAssignment2,
                                       originAssignment[place, ])
          }
          else {
            if (captured[i] == "origin" && capturePoint[i]){
              ass <- sf::st_intersects(x = originPoints[i,], y = originSites,
                                       sparse = FALSE)
              if (sum(ass)<1){
                warning("Not all origin capture locations are within originSites. Assigning to closest site.\n",
                        "Affects animal: ", i)
                ass <- matrix(0, 1, dimOAss[2])
                ass[sf::st_nearest_feature(x = originPoints[i,],
                                              y = originSites)] <- 1

              }
              else if (sum(ass) > 1){
                warning("Overlapping originSites may cause issues\n")
                ass0 <- sf::st_intersects(x = originPoints[i,], y = originSites,
                                          sparse = TRUE)
                ass <- matrix(0, 1, dimOAss[2])
                ass[ass0[1]] <- 1
              }
              originAssignment2 <- rbind(originAssignment2, ass)
            }
            else
              originAssignment2 <- rbind(originAssignment2,
                                         dummyAss)
          }

        }
      }
      else
        originAssignment2 <- originAssignment
    }
  }
  else {
    originAssignment2 <- originAssignment
    targetAssignment2 <- targetAssignment
  }
  return(list(originPoints = originPoints2, targetPoints = targetPoints2,
              originAssignment = originAssignment2,
              targetAssignment = targetAssignment2,
              isGL = isGL, isTelemetry = isTelemetry, isRaster = isRaster,
              isProb = isProb, captured = captured,
              # originRaster = originRaster2,
              originRasterXYZ = originRasterXYZ2,
              originSingleCell = originSingleCell,
              # targetRaster = targetRaster2,
              targetRasterXYZ = targetRasterXYZ2,
              targetSingleCell = targetSingleCell))
}

# Return randomly sampled points from sites, based on assignment
randomPoints <- function(sites, assignment, geomName = "geometry") {
  nSites <- nrow(sites)
  assignmentSum <- c(table(factor(assignment, levels = 1:nSites)))
  points <- sf::st_sample(sites, assignmentSum)
  points <- sf::st_as_sf(points)
  points <- points[rank(assignment, ties.method = "first"), ]
  names(points)[1] <- geomName
  sf::st_geometry(points) <- geomName
  return(points)
}

# Summary functions
summarize3D <- function(array3D, nSites1, nSites2, names1,
                         names2, point = NULL, alpha = 0.05) {
  meanA <- apply(array3D, 2:3, mean)
  medianA <- apply(array3D, 2:3, median)
  seA <- apply(array3D, 2:3, sd)
  simpleCIA <- apply(array3D, 2:3, quantile, probs = c(alpha/2, 1-alpha/2),
                       na.rm=TRUE)
  dimnames(simpleCIA)[[1]] <- c("lower", "upper")
  matrixA <- array(c(array3D), c(dim(array3D)[1], nSites1 * nSites2),
                      list(NULL, paste(rep(names1, nSites2),
                                       rep(names2, each = nSites1),
                                       sep = "#")))
  mcmcA <- coda::as.mcmc(matrixA)
  hpdCI <- coda::HPDinterval(mcmcA, 1-alpha)
  hpdCI <- array(hpdCI, c(nSites1, nSites2, 2),
                 list(names1, names2, c("lower", "upper")))
  hpdCI <- aperm(hpdCI, c(3, 1, 2))
  bcCIA <- array(NA, dim = c(2, nSites1, nSites2),
                   dimnames = list(c("lower", "upper"), names1, names2))
  for (i in 1:nSites1) {
    for (j in 1:nSites2) {
      psi.z0 <- qnorm(sum(array3D[, i, j] < meanA[i, j], na.rm = TRUE) /
                        length(which(!is.na(array3D[, i, j]))))
      bcCIA[ , i, j] <- quantile(array3D[, i, j],
                                   pnorm(2 * psi.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                   na.rm=TRUE, names = FALSE)
    }
  }
  return(list(sample = array3D, mean = meanA, se = seA,
              simpleCI = simpleCIA, bcCI = bcCIA, hpdCI = hpdCI,
              median = medianA, point = point))
}

summarizeAbund <- function(abund, nSites1, names1, pointAbund = NULL,
                           alpha = 0.05) {
  meanAbund <- apply(abund, 2, mean, na.rm=TRUE)
  medianAbund <- apply(abund, 2, median, na.rm=TRUE)
  seAbund <- apply(abund, 2, sd, na.rm=TRUE)
  # Calculate confidence intervals using quantiles of sampled MC
  simpleCIAbund <- apply(abund, 2, quantile,
                         probs = c(alpha/2, 1-alpha/2),
                         na.rm=TRUE)
  dimnames(simpleCIAbund)[[1]] <- c("lower", "upper")
  bcCIAbund <- array(NA, dim = c(2, nSites1),
                     dimnames = list(c("lower", "upper"), names1))
  for (i in 1:nSites1) {
    abund.z0 <- qnorm(sum(abund[, i] < meanAbund[i], na.rm = TRUE) /
                        length(which(!is.na(abund[, i]))))
    bcCIAbund[ , i] <- quantile(abund[, i],
                                pnorm(2 * abund.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                na.rm=TRUE, names = FALSE)
  }
  return(list(sample = abund, mean = meanAbund, se = seAbund,
              simpleCI = simpleCIAbund, bcCI = bcCIAbund, median = medianAbund,
              point = pointAbund))
}

# Currently called through reverseTransition (in calcConnectivity.R)
reverseEstTransition <- function(psi, originRelAbund, pi = NULL,
                                  originSites=NULL, targetSites=NULL,
                                  originNames = NULL, targetNames = NULL,
                                  nSamples = 1000, row0 = 0, alpha = 0.05) {
  if (is.null(psi)) {
    if (is.array(pi)) {
      if (length(dim(pi))!=3)
        stop('Pi should either be 2-(for fixed migratory combination probabilities) or 3-dimensional array')
      piIn <- pi
    }
    else if (inherits(pi, "estPi")) {
      if (is.null(originNames)) {
        originNames <- pi$input$originNames
      }
      if (is.null(targetNames)) {
        targetNames <- pi$input$targetNames
      }
      piIn <- pi
      pi <- pi$pi$sample
    }
    nOriginSites <- dim(pi)[2]
    nTargetSites <- dim(pi)[3]
    if (is.null(originNames)) {
      originNames <- dimnames(pi)[[2]]
    }
    if (is.null(targetNames)) {
      targetNames <- dimnames(pi)[[3]]
    }
    piBase <- apply(pi, 2:3, mean)
    if (dim(pi)[1]>=nSamples)
      piSamples <- round(seq(from = 1, to = dim(pi)[1],
                             length.out = nSamples))
    else
      piSamples <- sample.int(dim(pi)[1], nSamples, replace = TRUE)
    pointRev <- reverseTransition(pi = piBase)
    pointGamma <- pointRev$gamma
    pointTargetAbund <- pointRev$targetRelAbund
    pointPsi <- pointRev$psi
    pointOriginAbund <- pointRev$originRelAbund
    gamma <- array(NA, c(nSamples, nTargetSites, nOriginSites),
                   dimnames = list(NULL, targetNames, originNames))
    targetRelAbund <- array(NA, c(nSamples, nTargetSites),
                            dimnames = list(NULL, targetNames))
    psi <- array(NA, c(nSamples, nOriginSites, nTargetSites),
                dimnames = list(NULL, originNames, targetNames))
    originRelAbund <- array(NA, c(nSamples, nOriginSites),
                            dimnames = list(NULL, originNames))
    for (i in 1:nSamples) {
      piNew <- pi[piSamples[i],,]
      # Reverse for new pi
      sampleRev <- reverseTransition(pi = piNew)
      gamma[i, , ] <- sampleRev$gamma
      targetRelAbund[i, ] <- sampleRev$targetRelAbund
      psi[i, , ] <- sampleRev$psi
      originRelAbund[i, ] <- sampleRev$originRelAbund
    }
    gammaReturn <- summarize3D(gamma, nTargetSites, nOriginSites, targetNames,
                               originNames, pointGamma, alpha)
    targetAbundReturn <- summarizeAbund(targetRelAbund, nTargetSites,
                                        targetNames, pointTargetAbund, alpha)
    psiReturn <- summarize3D(psi, nOriginSites, nTargetSites, originNames,
                             targetNames, pointPsi, alpha)
    originAbundReturn <- summarizeAbund(originRelAbund, nOriginSites,
                                        originNames, pointOriginAbund, alpha)
    rev <- list(gamma = gammaReturn,
                targetRelAbund = targetAbundReturn,
                psi = psiReturn, originRekAbund = originAbundReturn,
                input = list(pi = piIn,
                             originSites = originSites, targetSites = targetSites,
                             originNames = originNames, targetNames = targetNames,
                             nSamples = nSamples, row0 = row0, alpha = alpha))
    class(rev) <- c("estGamma", "estTargetRelAbund", "estPsi",
                    "estOriginRelAbund", "estMigConnectivity")
    return(rev)
  }
  if (is.matrix(psi)) {
    psiFixed <- TRUE
    psiVCV <- NULL
    nOriginSites <- dim(psi)[1]
    nTargetSites <- dim(psi)[2]
    if (is.null(originNames)) {
      originNames <- dimnames(psi)[[1]]
    }
    if (is.null(targetNames)) {
      targetNames <- dimnames(psi)[[2]]
    }
    psiBase <- psi
    psiIn <- psi
  }
  else if (inherits(psi, "mark")) {
    psiFixed <- FALSE
    if (!is.numeric(originSites) || !is.numeric(targetSites))
      stop('Must specify which RMark Psi parameters represent origin and target sites')
    nOriginSites <- length(originSites)
    nTargetSites <- length(targetSites)
    psiVCV <- psi$results$beta.vcv
    psiBase <- RMark::TransitionMatrix(RMark::get.real(psi, "Psi",
                                                       se=TRUE))[originSites,
                                                                 targetSites]
    if (any(diag(psi$results$beta.vcv) < 0))
      stop("Can't sample model, negative beta variances")
    psiIn <- psi
  }
  else if (is.array(psi)) {
    if (length(dim(psi))!=3)
      stop('Psi should either be 2-(for fixed transition probabilities) or 3-dimensional array')
    psiFixed <- FALSE
    nOriginSites <- dim(psi)[2]
    nTargetSites <- dim(psi)[3]
    if (is.null(originNames)) {
      originNames <- dimnames(psi)[[2]]
    }
    if (is.null(targetNames)) {
      targetNames <- dimnames(psi)[[3]]
    }
    psiBase <- apply(psi, 2:3, mean)
    psiVCV <- NULL
    if (dim(psi)[1]>=nSamples)
      psiSamples <- round(seq(from = 1, to = dim(psi)[1],
                              length.out = nSamples))
    else
      psiSamples <- sample.int(dim(psi)[1], nSamples, replace = TRUE)
    psiIn <- psi
  }
  else if (inherits(psi, "estPsi") || inherits(psi, "estMC")) {
    psiFixed <- FALSE
    if (is.null(originNames)) {
      originNames <- psi$input$originNames
    }
    if (is.null(targetNames)) {
      targetNames <- psi$input$targetNames
    }
    psiIn <- psi
    psi <- psi$psi$sample
    nOriginSites <- dim(psi)[2]
    nTargetSites <- dim(psi)[3]
    psiBase <- apply(psi, 2:3, mean)
    psiVCV <- NULL
    if (dim(psi)[1]>=nSamples)
      psiSamples <- round(seq(from = 1, to = dim(psi)[1],
                              length.out = nSamples))
    else
      psiSamples <- sample.int(dim(psi)[1], nSamples, replace = TRUE)
  }
  originRelAbundIn <- originRelAbund
  if (coda::is.mcmc(originRelAbund) || coda::is.mcmc.list(originRelAbund)) {
    originRelAbund <- as.matrix(originRelAbund)
  }
  if (is.matrix(originRelAbund) && dim(originRelAbund)[1]>1) {
    abundFixed <- FALSE
    if (dim(originRelAbund)[2]>nOriginSites)
      abundParams <- paste('relN[', 1:nOriginSites, ']', sep='')
    else if (dim(originRelAbund)[2]==nOriginSites)
      abundParams <- 1:nOriginSites
    else
      stop('Number of origin sites must be constant between psi and abundance')
    if (dim(originRelAbund)[1] >= nSamples)
      abundRows <- round(seq(from = row0 + 1, to = dim(originRelAbund)[1],
                                    length.out = nSamples))
    else
      abundRows <- sample((row0 + 1):dim(originRelAbund)[1], nSamples,
                          replace = TRUE)
    originRelAbund <- as.matrix(originRelAbund[abundRows, abundParams])
    abundBase <- colMeans(originRelAbund)
  }
  else {
    abundFixed <- TRUE
    if (length(originRelAbund)!=nOriginSites)
      stop('Number of origin sites must be constant between psi and abundance')
    abundBase <- originRelAbund
  }
  pointRev <- reversePsiRelAbund(psi = psiBase, originRelAbund = abundBase)
  pointGamma <- pointRev$gamma
  pointAbund <- pointRev$targetRelAbund
  pointPi <- pointRev$pi
  gamma <- array(NA, c(nSamples, nTargetSites, nOriginSites),
                 dimnames = list(NULL, targetNames, originNames))
  targetRelAbund <- array(NA, c(nSamples, nTargetSites),
                          dimnames = list(NULL, targetNames))
  pi <- array(NA, c(nSamples, nOriginSites, nTargetSites),
              dimnames = list(NULL, originNames, targetNames))
  for (i in 1:nSamples) {
    # Generate random transition probability matrices
    if (psiFixed)
      psiNew <- psiBase
    else if (is.null(psiVCV))
      psiNew <- psi[psiSamples[i],,]
    else
      psiNew <- makePsiRand(psi, originSites, targetSites)
    # Point estimates of breeding densities
    if (abundFixed)
      abundNew <- abundBase
    else
      abundNew <- originRelAbund[i, abundParams]
    # Reverse for new psis and relative breeding densities
    sampleRev <- reversePsiRelAbund(psi = psiNew, originRelAbund = abundNew)
    gamma[i, , ] <- sampleRev$gamma
    targetRelAbund[i, ] <- sampleRev$targetRelAbund
    pi[i, , ] <- sampleRev$pi
  }
  gammaReturn <- summarize3D(gamma, nTargetSites, nOriginSites, targetNames,
                             originNames, pointGamma, alpha)
  targetAbundReturn <- summarizeAbund(targetRelAbund, nTargetSites,
                                      targetNames, pointAbund, alpha)
  piReturn <- summarize3D(pi, nOriginSites, nTargetSites, originNames,
                           targetNames, pointPi, alpha)
  rev <- list(gamma = gammaReturn,
              targetRelAbund = targetAbundReturn,
              pi = piReturn,
              input = list(psi = psiIn, originRelAbund = originRelAbundIn,
                           originSites = originSites, targetSites = targetSites,
                           originNames = originNames, targetNames = targetNames,
                           nSamples = nSamples, row0 = row0, alpha = alpha))
  class(rev) <- c("estGamma", "estTargetRelAbund", "estPi", "estMigConnectivity")
  return(rev)

}

assignmentToMatrix <- function(assignment, nSites) {
  nInds <- length(assignment)
  assignMatrix <- array(0, c(nInds, nSites))
  for (i in 1:nInds)
    if (!is.na(assignment[i]))
      assignMatrix[i, assignment[i]] <- 1
  return(assignMatrix)
}

rasterToProb <- function(originSites = NULL, targetSites = NULL,
                         originAssignment = NULL, targetAssignment = NULL,
                         isRaster = FALSE, isProb = FALSE,
                         captured = "origin",
                         originRaster = NULL, targetRaster = NULL) {
  failTarget <- failOrigin <- FALSE
  nInds <- length(captured)
  if (is.null(originSites))
    failOrigin <- TRUE
  else
    nOriginSites <- nrow(originSites)
  if (is.null(targetSites))
    failTarget <- TRUE
  else
    nTargetSites <- nrow(targetSites)
  inds <- isRaster & !isProb
  capTarget <- captured=="target"
  capOrigin <- captured=="origin"
  whichTarget <- which(inds & !capTarget)
  whichOrigin <- which(inds & !capOrigin)
  if (length(whichTarget>0) && !is.null(targetSites)) {
    targetSites <- sf::st_transform(targetSites, terra::crs(targetRaster))
    targetSum <- terra::extract(targetRaster,
                                terra::vect(targetSites), fun = sum,
                                na.rm = TRUE, ID = FALSE)
    targetSum <- t(targetSum[, whichTarget, drop = FALSE])
    targetSum <- sweep(targetSum, 1, rowSums(targetSum), "/")
    if (is.null(targetAssignment)){
      targetAssignment <- array(0, c(nInds, nTargetSites))
    }
    else if (!is.array(targetAssignment) || length(dim(targetAssignment))<2)
      targetAssignment <- assignmentToMatrix(targetAssignment, nTargetSites)
    targetAssignment[whichTarget, ] <- targetSum
    isRaster[whichTarget] <- FALSE
    isProb[whichTarget] <- TRUE
  }
  if (length(whichOrigin>0) && !is.null(originSites)) {
    originSites <- sf::st_transform(originSites, terra::crs(originRaster))
    originSum <- terra::extract(originRaster,
                                terra::vect(originSites), fun = sum,
                                na.rm = TRUE, ID = FALSE)
    originSum <- t(originSum[, whichOrigin, drop = FALSE])
    originSum <- sweep(originSum, 1, rowSums(originSum), "/")
    if (is.null(originAssignment)){
      originAssignment <- array(0, c(nInds, nOriginSites))
    }
    else if (!is.array(originAssignment) || length(dim(originAssignment))<2)
      originAssignment <- assignmentToMatrix(originAssignment, nOriginSites)
    originAssignment[whichOrigin, ] <- originSum
    isRaster[whichOrigin] <- FALSE
    isProb[whichOrigin] <- TRUE
  }
  return(list(originAssignment = originAssignment,
              targetAssignment = targetAssignment,
              failOrigin = failOrigin, failTarget = failTarget,
              isRaster = isRaster, isProb = isProb))
}

assignRasterStats <- function(theRaster) {
  if (is.null(theRaster)){
    PointsAssigned <- FALSE
    SingleCell <- NULL
    RasterXYZ <- NULL
    RasterXYZcrs <- NULL
    outRaster <- NULL
  }
  else {
    if (inherits(theRaster, "isoAssign")) {
      if (!inherits(theRaster$probassign, "SpatRaster"))
        theRaster$probassign <- terra::rast(theRaster$probassign)
      PointsAssigned <- !is.null(theRaster$SingleCell) &&
                                  !is.null(dim(theRaster$SingleCell))
      RasterXYZ <- terra::as.data.frame(theRaster$probassign, xy = TRUE)
      RasterXYZcrs <- terra::crs(theRaster$probassign)
      SingleCell <- theRaster$SingleCell
      outRaster <- theRaster$probassign
    }
    else {
      if (!inherits(theRaster, "SpatRaster")) {
        theRaster <- terra::rast(theRaster)
      }
      if (is.na(terra::crs(theRaster))){
        stop("Please provide a projection (crs) for each raster\n")
      }
      RasterXYZ <- terra::as.data.frame(theRaster, xy = TRUE)
      RasterXYZcrs <- terra::crs(theRaster)
      SingleCell <- NULL
      PointsAssigned <- FALSE
      outRaster <- theRaster
    }
  }
  return(list(Raster = outRaster, PointsAssigned = PointsAssigned,
              SingleCell = SingleCell, RasterXYZ = RasterXYZ,
              RasterXYZcrs = RasterXYZcrs))
}

# function to make an odds ratio (likely vs unlikely) assignment
oddsFun <- function(x, odds = odds){
  predict(smooth.spline(x = cumsum(sort(x)),
                        sort(x),
                        spar = 0.1),(1-odds))$y
}

propSpatRaster <- function(assignments) {
  if (terra::nlyr(assignments)<2){
    assign2prob  <- assignments /
      unlist(c(terra::global(assignments, fun = "sum", na.rm = T)))
    return(assign2prob)
  }
  summy <- terra::global(assignments, "sum", na.rm = TRUE)
  test <- lapply(summy[,1],
                 function(x) terra::rast(nrows = terra::nrow(assignments),
                                         ncols = terra::ncol(assignments),
                                         xmin = terra::xmin(assignments),
                                         xmax = terra::xmax(assignments),
                                         ymin = terra::ymin(assignments),
                                         ymax = terra::ymax(assignments),
                                         vals = x))
  test <- terra::rast(test)
  assign2prob <- assignments / test
  return(assign2prob)
}
