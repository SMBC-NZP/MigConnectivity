#' @import stats

###############################################################################
# Estimate MC from abundance and transition probability estimates.
###############################################################################
estMCCmrAbund <- function(originDist, targetDist, originRelAbund, psi,
                          sampleSize = NULL,
                          originSites=NULL, targetSites=NULL,
                          nSamples = 1000, row0 = 0, verbose=0,
                          alpha = 0.05,
                          approxSigTest = F, sigConst = 0) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  absAbund <- !is.null(sampleSize)
  if (coda::is.mcmc(originRelAbund)) {
    abundFixed <- FALSE
    abundParams <- paste('relN[', 1:nOrigin, ']', sep='')
    abundBase <- colMeans(originRelAbund[row0 + 1:nSamples, abundParams])
  }
  else {
    abundFixed <- TRUE
    if (length(originRelAbund)!=nOrigin)
      stop('Number of origin sites must be constant between distance matrix and abundance')
    abundBase <- originRelAbund
  }
  if (is.matrix(psi)) {
    psiFixed <- TRUE
    if (nrow(psi)!=nOrigin || ncol(psi)!=nTarget)
      stop('Size of psi matrix must be consistant with distance matrices')
    psiBase <- psi
  }
  else {
    psiFixed <- FALSE
    if (!is.numeric(originSites) || !is.numeric(targetSites))
      stop('Must specify which RMark Psi parameters represent origin and target sites')
    psiBase <- RMark::TransitionMatrix(RMark::get.real(psi, "Psi",
                                                       se=TRUE))[originSites,
                                                              targetSites]
    if (any(diag(psi$results$beta.vcv) < 0))
      stop("Can't sample model, negative beta variances")
  }
  pointMC <- ifelse(absAbund,
                    calcMC(originDist, targetDist, originRelAbund = abundBase,
                           psi = psiBase, sampleSize = sampleSize),
                    calcMC(originDist, targetDist, originRelAbund = abundBase,
                           psi = psiBase))
  sampleMC <- vector('numeric', nSamples)
  psi.array <- array(0, c(nSamples, nOrigin, nTarget))
  for (i in 1:nSamples) {
    if (verbose > 1 || verbose == 1 && i %% 100 == 0)
      cat("\tSample", i, "of", nSamples, "at", date(), "\n")
    # Generate random transition probability matrices
    # resampling from each files, 1000 samples from each of the 100 files
    # using the estimates of psi and uncertainty from the mark file to generate a new estimate of psi
    # choose a random psi from that interval of the distribution
    if (psiFixed)
      psiNew <- psiBase
    else
      psiNew <- makePsiRand(psi, originSites, targetSites)
    psi.array[i, , ] <- psiNew
    # Point estimates of breeding densities
    if (abundFixed)
      abundNew <- abundBase
    else
      abundNew <- originRelAbund[row0 + i, abundParams]
    # Calculate MC for new psis and relative breeding densities
    sampleMC[i] <- ifelse(absAbund,
                          calcMC(originDist, targetDist, originRelAbund = abundNew,
                                 psi = psiNew, sampleSize = sampleSize),
                          calcMC(originDist, targetDist, originRelAbund = abundNew,
                                 psi = psiNew))
    if (verbose > 1 || verbose == 1 && i %% 10 == 0)
      cat(" MC mean:", mean(sampleMC, na.rm=TRUE),
          "SD:", sd(sampleMC, na.rm=TRUE),
          "low quantile:", quantile(sampleMC, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(sampleMC, 1-alpha/2, na.rm=TRUE), "\n")
  }
  meanMC <- mean(sampleMC, na.rm=TRUE)
  medianMC <- median(sampleMC, na.rm=TRUE)
  seMC <- sd(sampleMC, na.rm=TRUE)
  # Calculate confidence intervals using quantiles of sampled MC
  simpleCI <- quantile(sampleMC, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                       names = F)
  z0 <- qnorm(sum((sampleMC)<meanMC)/nSamples)
  bcCI <- quantile(sampleMC, pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, type = 8, names = F)
  MC.mcmc <- coda::as.mcmc(sampleMC) # Ha!
  hpdCI <- as.vector(coda::HPDinterval(MC.mcmc, 1-alpha))
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (pointMC > sigConst)
      simpleP <- sum(sampleMC < sigConst) / nSamples
    else
      simpleP <- sum(sampleMC > sigConst) / nSamples
    if (simpleP == 0)
      simpleP <- 0.5 / nSamples
    bcP <- pnorm(qnorm(simpleP) - 2 * z0)
    if (pointMC < sigConst)
      bcP <- 1 - bcP
  }
  return(list(sampleMC=sampleMC, samplePsi = psi.array, pointPsi = psiBase,
              pointMC=pointMC, meanMC=meanMC,
              medianMC=medianMC, seMC=seMC, simpleCI=simpleCI,
              bcCI=bcCI, hpdCI=hpdCI, simpleP = simpleP, bcP = bcP,
              sampleCorr = NULL, pointCorr = NULL,
              meanCorr = NULL, medianCorr = NULL, seCorr=NULL,
              simpleCICorr=NULL, bcCICorr=NULL, inputSampleSize = sampleSize,
              alpha = alpha, sigConst = sigConst))
}

###############################################################################
# Resampling of uncertainty from geolocators and/or GPS data
###############################################################################
estMCGlGps <- function(originDist, targetDist, originRelAbund, isGL,
                       geoBias, geoVCov, targetPoints, targetSites,
                       sampleSize = NULL,
                       targetAssignment=NULL,
                       originPoints=NULL, originSites=NULL,
                       originAssignment=NULL, originNames=NULL,
                       targetNames=NULL, nBoot = 1000, verbose=0,
                       nSim = 1000, calcCorr=TRUE, alpha = 0.05,
                       approxSigTest = F, sigConst = 0,
            resampleProjection = MigConnectivity::projections$EquidistConic,
                       maxTries = 300) {

  # Input checking and assignment
  if (!(verbose %in% 0:3))
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = number of draws")
  if (length(geoBias)!=2 && (any(isGL) || calcCorr))
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in resampleProjection units, default meters)")
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = F)) && (any(isGL) || calcCorr))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in resampleProjection units, default meters)")
  if ((is.null(originPoints) || is.null(originSites)) &&
      is.null(originAssignment))
    stop("Need to define either originAssignment or originSites and originPoints")
  if (calcCorr && is.null(originPoints))
    stop('If calcCorr is TRUE, need to define originPoints')
  if ((is.null(targetPoints) || is.null(targetSites)) &&
      is.null(targetAssignment))
    stop("Need to define either targetAssignment or targetSites and targetPoints")
  nAnimals <- max(length(targetPoints), length(targetAssignment))
  if (length(isGL)==1)
    isGL <- rep(isGL, nAnimals)
  if(class(originSites)=="SpatialPolygonsDataFrame"){
    originSites <- sp::SpatialPolygons(originSites@polygons,proj4string=originSites@proj4string)}
  if(class(targetSites)=="SpatialPolygonsDataFrame"){
    targetSites <- sp::SpatialPolygons(targetSites@polygons,proj4string=targetSites@proj4string)}
  if (is.null(originAssignment))
    originAssignment <- sp::over(originPoints, originSites)
  if (is.null(targetAssignment))
    targetAssignment <- sp::over(targetPoints, targetSites)
  nOriginSites <- length(unique(originAssignment))
  nTargetSites <- ifelse(is.null(targetSites), nrow(targetDist), length(targetSites))
  if (length(targetPoints)!=nAnimals && length(targetAssignment)!=nAnimals ||
      length(originAssignment)!=nAnimals)
    stop("isGL should be the same length as originAssignment/originPoints and targetPoints/targetAssignment (number of animals)")
  if (any(is.na(originAssignment)))
    stop("NAs in origin sites (make sure all points fall within polygons)")
  if (length(originRelAbund)!=nOriginSites || sum(originRelAbund)!=1)
    stop('originRelAbund should be vector with length number of origin sites that sums to 1')
  if(!is.null(originPoints))
    if(is.na(raster::projection(originPoints)) || is.na(raster::projection(originSites))) {
      stop('Coordinate system definition needed for originSites & originPoints')
    }
  if(is.na(raster::projection(targetSites)) || is.na(raster::projection(targetPoints))){
    stop('Coordinate system definition needed for targetSites & targetPoints')
  }
  targetPoints <- sp::spTransform(targetPoints, sp::CRS(resampleProjection))
  targetSites <- sp::spTransform(targetSites, sp::CRS(resampleProjection))
  if (is.null(targetNames))
    targetNames <- names(targetSites)
  if (is.null(originNames))
    originNames <- names(originSites)
  if (is.null(sampleSize))
    sampleSize <- nAnimals
  if (dim(originDist)!=rep(nOriginSites,2) ||
      dim(targetDist)!=rep(nTargetSites,2))
    stop('Distance matrices should be square with same number of sites of each type as assignments/points (with distances in meters)')
  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))


  MC <- corr <- rep(NA, nBoot)

  # determine the number of animals from the input data
  nAnimals <- length(originAssignment)

  # Point estimate of MC
  pointSites <- array(0, c(nOriginSites, nTargetSites),
                       dimnames = list(originNames, targetNames))

  for(i in 1:nAnimals)
      pointSites[originAssignment[i], targetAssignment[i]] <-
      pointSites[originAssignment[i], targetAssignment[i]] + 1

  pointPsi <- prop.table(pointSites, 1)

  pointMC <- calcMC(originDist, targetDist, originRelAbund, pointPsi,
                    sampleSize = sampleSize)

  if (calcCorr) {
    targetDist1 <- matrix(NA, nAnimals, nAnimals)

    targetDist1[lower.tri(targetDist1)] <- 1

    distIndices <- which(!is.na(targetDist1), arr.ind = TRUE)

    # project target points to WGS #
    targetPoints2 <- sp::spTransform(targetPoints, sp::CRS(MigConnectivity::projections$WGS84))

    targetDist0 <- geosphere::distVincentyEllipsoid(targetPoints2[distIndices[,'row'],],targetPoints2[distIndices[,'col'],])

    targetDist1[lower.tri(targetDist1)] <- targetDist0

    originPoints2 <- sp::spTransform(originPoints, sp::CRS(MigConnectivity::projections$WGS84))

    originDistStart <- matrix(geosphere::distVincentyEllipsoid(originPoints2[rep(1:nAnimals, nAnimals)], originPoints2[rep(1:nAnimals,each=nAnimals)]),
                                    nAnimals, nAnimals)

    originDist1 <- originDistStart[1:nAnimals, 1:nAnimals]

    pointCorr <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = TRUE)$correlation
  }

  boot <- 1
  while (boot <= nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have birds from every origin site
    origin.sample <- 'Filler' # Start with one origin site
    while (length(unique(origin.sample)) < nOriginSites) { #2
      # Sample individual animals with replacement
      animal.sample <- sample.int(nAnimals, replace=TRUE)
      # Get origin points for those animals
      if (calcCorr)
        origin.point.sample <- originPoints[animal.sample]
      # Get origin population for each animal sampled
      origin.sample <- originAssignment[animal.sample]
    }
    tSamp <- targetSample(isGL = isGL, geoBias = geoBias, geoVCov = geoVCov,
                          targetPoints = targetPoints, animal.sample = animal.sample,
                          targetSites = targetSites, targetAssignment = targetAssignment,
                          resampleProjection = resampleProjection, nSim = nSim,
                          maxTries = maxTries)
    target.sample <- tSamp$target.sample
    target.point.sample <- tSamp$target.point.sample
    if (verbose > 2)
      cat(' ', tSamp$draws, 'draws (of length', nSim, 'and of', maxTries, 'possible).\n')
    # Now that we have breeding and non-breeding site for point, add to transition count matrix
    sites <- table(origin.sample, target.sample)
    sites.array[boot, as.integer(rownames(sites)), as.integer(colnames(sites))] <- sites
    # Create psi matrix as proportion of those from each breeding site that went to each NB site
    psi.array[boot, , ] <- prop.table(sites.array[boot, , ], 1)
    # Calculate MC from that psi matrix
    MC[boot] <- calcMC(originDist, targetDist, originRelAbund,
                       psi.array[boot, , ], sampleSize)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=TRUE), "SD:", sd(MC, na.rm=TRUE),
          "low quantile:", quantile(MC, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=TRUE), "\n")
    if (calcCorr) {
      originDist1 <- originDistStart[animal.sample, animal.sample]
      target.point.sample <- sp::SpatialPoints(target.point.sample,sp::CRS(resampleProjection))
      target.point.sample2 <- sp::spTransform(target.point.sample,sp::CRS(MigConnectivity::projections$WGS84))
      targetDist0 <- geosphere::distVincentyEllipsoid(target.point.sample2[distIndices[,'row']],
                                           target.point.sample2[distIndices[,'col']])
      targetDist1[lower.tri(targetDist1)] <- targetDist0
      corr[boot] <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = TRUE)$correlation
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=TRUE), "SD:", sd(corr, na.rm=TRUE),
            "low quantile:", quantile(corr, alpha/2, na.rm=TRUE),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=TRUE), "\n")
    }
    if (!is.na(MC[boot]))
      boot <- boot + 1
  }
  MC.z0 <- qnorm(sum(MC<mean(MC, na.rm = T), na.rm = T)/length(which(!is.na(MC))))
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, type = 8, names = F)
  MC.mcmc <- coda::as.mcmc(MC) # Ha!
  hpdCI <- as.vector(coda::HPDinterval(MC.mcmc, 1-alpha))
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (pointMC > sigConst)
      simpleP <- sum(MC < sigConst) / nBoot
    else
      simpleP <- sum(MC > sigConst) / nBoot
    if (simpleP == 0)
      simpleP <- 0.5 / nBoot
    bcP <- pnorm(qnorm(simpleP) - 2 * MC.z0)
    if (pointMC < sigConst)
      bcP <- 1 - bcP
  }
  if (calcCorr) {
    meanCorr <- mean(corr, na.rm=TRUE)
    medianCorr <- median(corr, na.rm=TRUE)
    seCorr <- sd(corr, na.rm=TRUE)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                             names = F)
    corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                           na.rm=TRUE, type = 8, names = F)
  } else
    pointCorr <- meanCorr <- medianCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  return(list(sampleMC = MC, samplePsi = psi.array,
              pointPsi = pointPsi, pointMC = pointMC, meanMC = mean(MC, na.rm=TRUE),
              medianMC = median(MC, na.rm=TRUE), seMC = sd(MC, na.rm=TRUE),
              simpleCI = quantile(MC, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                                  type = 8, names = F),
              bcCI = bcCI, hpdCI = hpdCI, simpleP = simpleP, bcP = bcP,
              sampleCorr = corr, pointCorr = pointCorr,
              meanCorr = meanCorr, medianCorr = medianCorr, seCorr=seCorr,
              simpleCICorr=simpleCICorr, bcCICorr=bcCICorr,
              inputSampleSize = sampleSize,
              alpha = alpha, sigConst = sigConst))
}

###############################################################################
#
#
# Resampling of uncertainty from intrinsic markers (Stable-Isotopes)
#
#
###############################################################################
estMCisotope <- function(targetDist=NULL,
                         originRelAbund,
                         targetIntrinsic,
                         targetSites = NULL,
                         sampleSize = NULL,
                         # targetAssignment=NULL,
                         originPoints=NULL,
                         originSites=NULL,
                         originDist = NULL,
                         originAssignment=NULL,
                         originNames=NULL,
                         targetNames=NULL,
                         nBoot = 1000,
                         verbose=0,
                         nSim = NULL,
                         calcCorr=TRUE,
                         alpha = 0.05,
                         approxSigTest = F,
                         sigConst = 0,
                         resampleProjection = MigConnectivity::projections$WGS84,
                         maxTries = 300) {

  # Input checking and assignment
  if (!(verbose %in% 0:3)){
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = number of draws")}

  if ((is.null(originPoints) || is.null(originSites)) && is.null(originAssignment)){
    stop("Need to define either originAssignment or originSites and originPoints")}

  if (calcCorr && is.null(originPoints)){
    stop('If calcCorr is TRUE, need to define originPoints')}

  if (is.null(targetIntrinsic)){
    stop("Need to define targetIntrinsic")}

  if (is.null(targetSites))
    targetSites <- targetIntrinsic$targetSites

  if (is.null(targetDist))
    targetDist <- distFromPos(rgeos::gCentroid(targetSites, byid = TRUE)@coords)


  if (!inherits(targetIntrinsic, 'isoAssign'))
    stop("targetIntrinsic should be output of isoAssign when isIntrinsic == TRUE")

  pointsAssigned <- !(is.null(targetIntrinsic$SingleCell) || is.na(targetIntrinsic$SingleCell))

  if(class(originSites)=="SpatialPolygonsDataFrame"){
    originSites <- sp::SpatialPolygons(originSites@polygons,proj4string=originSites@proj4string)}
  if(class(targetSites)=="SpatialPolygonsDataFrame"){
    targetSites <- sp::SpatialPolygons(targetSites@polygons,proj4string=targetSites@proj4string)}
  if (is.null(originAssignment))
    originAssignment <- sp::over(originPoints, originSites)

  nAnimals <- ifelse(pointsAssigned, dim(targetIntrinsic$SingleCell)[3], dim(targetIntrinsic$probassign)[3])
  targetSites <- sp::spTransform(targetSites, sp::CRS(resampleProjection))

  if (length(originAssignment)!=nAnimals)
    stop("originAssignment/originPoints should be the same length as targetIntrinsic (number of animals)")

  pointsInSites <- FALSE
  if (pointsAssigned && !is.null(targetSites)) {
    nSamples <- dim(targetIntrinsic$SingleCell)[1]
    targCon <- array(NA, c(nSamples, nAnimals))
    for(i in 1:nSamples) {
      targCon[i, ] <- sp::over(sp::SpatialPoints(t(targetIntrinsic$SingleCell[i, , ]),
                                         proj4string = sp::CRS(targetSites@proj4string@projargs)),
                           targetSites)
    }
    if (!any(is.na(targCon)))
      pointsInSites <- TRUE
    else if (verbose > 0)
      cat('Single cell assignment points supplied, but some points (proportion',
        sum(is.na(targCon))/length(targCon), ') not in targetSites\n')
  }
  else
    targCon <- NULL

  nOriginSites <- length(unique(originAssignment))
  nTargetSites <- ifelse(is.null(targetSites), nrow(targetDist), length(targetSites))


  if (any(is.na(originAssignment)))
    stop("NAs in origin sites (make sure all points fall within polygons)")

  if (length(originRelAbund)!=nOriginSites || sum(originRelAbund)!=1)
    stop('originRelAbund should be vector with length number of origin sites that sums to 1')

  if(!is.null(originPoints))
    if(is.na(raster::projection(originPoints)) || is.na(raster::projection(originSites))) {
      stop('Coordinate system definition needed for originSites & originPoints')
    }



  if (is.null(targetNames))
    targetNames <- names(targetSites)

  if (is.null(originNames))
    originNames <- names(originSites)

  if (is.null(sampleSize))
    sampleSize <- nAnimals

  if (dim(originDist)!=rep(nOriginSites,2) ||
      dim(targetDist)!=rep(nTargetSites,2))
    stop('Distance matrices should be square with same number of sites of each type as assignments/points (with distances in meters)')


  if (calcCorr) {
    targetDist1 <- matrix(NA, nAnimals, nAnimals)

    targetDist1[lower.tri(targetDist1)] <- 1

    distIndices <- which(!is.na(targetDist1), arr.ind = TRUE)

    originPoints2 <- sp::spTransform(originPoints, sp::CRS(MigConnectivity::projections$WGS84))

    originDistStart <- matrix(geosphere::distVincentyEllipsoid(originPoints2[rep(1:nAnimals, nAnimals)], originPoints2[rep(1:nAnimals,each=nAnimals)]),
                                    nAnimals, nAnimals)

  }



  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))


  MC <- corr <- rep(NA, nBoot)


  boot <- 1
  while (boot <= nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have birds from every origin site
    origin.sample <- 'Filler' # Start with one origin site
    while (length(unique(origin.sample)) < nOriginSites) { #2
      # Sample individual animals with replacement
      animal.sample <- sample.int(nAnimals, replace=TRUE)
      # Get origin points for those animals
      if (calcCorr)
        origin.point.sample <- originPoints[animal.sample]
      # Get origin population for each animal sampled
      origin.sample <- originAssignment[animal.sample]
    }
    if (verbose > 2)
      cat("Origin sample complete, making target sample with",
          ifelse(pointsAssigned, "points assigned,", "no points assigned,"),
          ifelse(pointsInSites, "all points in sites\n", "not all points in sites\n"))
    # Resample from points for each animal
    tSamp <- targetSampleIsotope(targetIntrinsic = targetIntrinsic,
                                 animal.sample = animal.sample,
                                 targetSites = targetSites,
                                 resampleProjection = resampleProjection, nSim = nSim,
                                 maxTries = maxTries, pointsAssigned = pointsAssigned,
                                 targCon = targCon, pointsInSites = pointsInSites)

    target.sample <- tSamp$target.sample
    target.point.sample <- tSamp$target.point.sample
    if (verbose > 2 & !pointsInSites)
      cat(' ', tSamp$draws, 'draws (of length', nSim, 'and of', maxTries, 'possible).\n')
    # Now that we have breeding and non-breeding site for point, add to transition count matrix
    sites <- table(origin.sample, target.sample)
    sites.array[boot, as.integer(rownames(sites)), as.integer(colnames(sites))] <- sites
    # Create psi matrix as proportion of those from each breeding site that went to each NB site
    psi.array[boot, , ] <- prop.table(sites.array[boot, , ], 1)
    # Calculate MC from that psi matrix
    MC[boot] <- calcMC(originDist, targetDist, originRelAbund,
                       psi.array[boot, , ], sampleSize)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=TRUE), "SD:", sd(MC, na.rm=TRUE),
          "low quantile:", quantile(MC, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=TRUE), "\n")
    if (calcCorr) {
      originDist1 <- originDistStart[animal.sample, animal.sample]
      target.point.sample <- sp::SpatialPoints(target.point.sample,sp::CRS(resampleProjection))
      target.point.sample2 <- sp::spTransform(target.point.sample,sp::CRS(MigConnectivity::projections$WGS84))
      targetDist0 <- geosphere::distVincentyEllipsoid(target.point.sample2[distIndices[,'row']],
                                                      target.point.sample2[distIndices[,'col']])
      targetDist1[lower.tri(targetDist1)] <- targetDist0
      corr[boot] <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = TRUE)$correlation
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=TRUE), "SD:", sd(corr, na.rm=TRUE),
            "low quantile:", quantile(corr, alpha/2, na.rm=TRUE),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=TRUE), "\n")
    }
    if (!is.na(MC[boot]))
      boot <- boot + 1
  }
  MC.z0 <- qnorm(sum(MC<mean(MC, na.rm = T), na.rm = T)/length(which(!is.na(MC))))
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                   na.rm=TRUE, type = 8, names = F)
  MC.mcmc <- coda::as.mcmc(MC) # Ha!
  hpdCI <- as.vector(coda::HPDinterval(MC.mcmc, 1-alpha))
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (mean(MC, na.rm=TRUE) > sigConst)
      simpleP <- sum(MC < sigConst) / nBoot
    else
      simpleP <- sum(MC > sigConst) / nBoot
    if (simpleP == 0)
      simpleP <- 0.5 / nBoot
    bcP <- pnorm(qnorm(simpleP) - 2 * MC.z0)
    if (mean(MC, na.rm=TRUE) < sigConst)
      bcP <- 1 - bcP
  }
  if (calcCorr) {
    meanCorr <- mean(corr, na.rm=TRUE)
    medianCorr <- median(corr, na.rm=TRUE)
    seCorr <- sd(corr, na.rm=TRUE)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                             names = F)
    corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                         na.rm=TRUE, type = 8, names = F)
  } else
    pointCorr <- meanCorr <- medianCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  return(list(sampleMC = MC, samplePsi = psi.array,
              pointPsi = NA, pointMC = NA, meanMC = mean(MC, na.rm=TRUE),
              medianMC = median(MC, na.rm=TRUE), seMC = sd(MC, na.rm=TRUE),
              simpleCI = quantile(MC, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                                  type = 8, names = F),
              bcCI = bcCI, hpdCI = hpdCI, simpleP = simpleP, bcP = bcP,
              sampleCorr = corr, pointCorr = NA,
              meanCorr = meanCorr, medianCorr = medianCorr, seCorr=seCorr,
              simpleCICorr=simpleCICorr, bcCICorr=bcCICorr,
              inputSampleSize = sampleSize,
              alpha = alpha, sigConst = sigConst))
}

###############################################################################
#' Estimate MC from abundance and/or transition probability estimates OR
#' geolocator and/or GPS data OR intrinsic markers.
#'
#' Resampling of uncertainty for MC from RMark psi matrix estimates and/or JAGS
#' relative abundance MCMC samples OR SpatialPoints geolocators and/or GPS
#' data OR intrinsic markers such as isotopes.
#'
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'  matrix.
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'  matrix.  Optional for intrinsic data.
#' @param originRelAbund Relative abundance estimates at B origin sites. Either
#'  a numeric vector of length B that sums to 1 or an mcmc object with
#'  \code{nSamples} rows  and columns including 'relN[1]' through 'relN[B]'.
#'  Currently, an mcmc object doesn't work with geolocator, GPS, or intrinsic
#'  data.
#' @param psi Transition probabilities between B origin and W target sites.
#'  Either a matrix with B rows and W columns where rows sum to 1 or a MARK
#'  object with estimates of transition probabilities.  If you are estimating
#'  MC from GPS, geolocator, or intrinsic data, leave this as NULL.
#' @param sampleSize Total sample size of animals that psi will be estimated from.
#'    Should be the number of animals released in one of the origin sites and
#'    observed in one of the target sites.  Optional, but recommended, unless
#'    you are estimating MC from GPS, geolocator, or intrinsic data (in which
#'    case the function can calculate it for you).
#' @param originSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are origin.  If using GPS, geolocator, or
#'  intrinsic data, this can be the geographic definition of sites in the
#'  release season.
#' @param targetSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are target.  If using GPS, geolocator, or
#'  intrinsic data, this must be the geographic definition of sites in the
#'  non-release season.  Optional for intrinsic data; if left out, the function
#'  will use the \code{targetSites} defined in \code{targetIntrinsic}.
#' @param originPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the release location of an animal.
#' @param targetPoints For GL or GPS data, a \code{SpatialPoints} object, with
#'    length number ofanimals tracked.  Each point indicates the point estimate
#'    location in the non-release season.
#' @param originAssignment Assignment of \code{originPoints} to release season
#'    sites. Integer vector with length number of animals tracked. Optional,
#'    but if using GL or GPS data, either \code{originAssignment} or
#'    \code{originSites} and \code{originPoints} should be defined.
#' @param targetAssignment Optional. Point estimate assignment of
#'    \code{targetPoints} to non-release season sites. Integer vector with
#'    length number of animals tracked.
#' @param originNames Optional. Vector of names for the release season sites.
#' @param targetNames Optional. Vector of names for the non-release season
#'    sites.
#' @param nSamples Number of times to resample \code{psi} and/or
#'    \code{originRelAbund} OR number of times to resample \code{targetPoints}
#'    for intrinsic data OR number of bootstrap runs for GL or GPS data. In
#'    the last two cases, animals are sampled with replacement for each. For all,
#'    the purpose is to estimate sampling uncertainty.
#' @param nSim Tuning parameter for GL or intrinsic data. Affects only the
#'    speed; 1000 seems to work well with our GL data and 10 for our intrinsic
#'    data, but your results may vary.  Should be integer > 0.
#' @param isGL Indicates whether or which animals were tracked with geolocators.
#'    Should be either single TRUE or FALSE value, or vector with length of
#'    number of animals tracked, with TRUE for animals in
#'    \code{targetPoints} with geolocators and FALSE for animals with GPS.
#' @param geoBias For GL data, vector of length 2 indicating expected bias
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters).
#' @param geoVCov For GL data, 2x2 matrix with expected variance/covariance
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters).
#' @param row0 If \code{originRelAbund} is an mcmc object, this can be set
#'  to 0 (default) or any greater integer to specify where to stop ignoring
#'  samples ("burn-in").
#' @param verbose 0 (default) to 3. 0 prints no output during run. 1 prints
#'  a line every 100 samples or bootstraps and a summary every 10.  2 prints a
#'  line and summary every sample or bootstrap. 3 also prints the number of
#'  draws (for tuning nSim for GL/intrinsic data only).
#' @param calcCorr In addition to MC, should function also estimate Mantel
#'    correlation between release and non-release locations (GPS or GL data
#'    only)?  Default is FALSE.
#' @param alpha Level for confidence/credible intervals provided.
#' @param approxSigTest Should function compute approximate one-sided
#'    significance tests (p-values) for MC from the bootstrap?  Default is
#'    FALSE.
#' @param sigConst Value to compare MC to in significance test.
#'    Default is 0.
#' @param resampleProjection Projection when sampling from geolocator
#'    bias/error. This projection needs units = m. Default is Equidistant
#'    Conic. The default setting preserves distances around latitude = 0 and
#'    longitude = 0. Other projections may work well, depending on the location
#'    of \code{targetSites}.  Ignored unless data are geolocator or GPS.
#' @param maxTries Maximum number of times to run a single GL/intrinsic
#'    bootstrap before exiting with an error.  Default is 300.  Set to NULL to
#'    never stop.  Thisparameter was added to prevent GL setups where some
#'    sample points never land on target sites from running indefinitely.
#' @param targetIntrinsic For intrinsic tracking data, the results of
#'    \code{isoAssign} or a similar function, of class \code{intrinsicAssign}.
#' @param isIntrinsic Logical indicating whether the animals are tracked via
#'    intrinsic marker (e.g. isotopes) or not.  Currently estMC will only estimate
#'    connectivity for all intrinsically marked animals or all extrinsic (e.g.,
#'    bands, GL, or GPS), so isIntrinsic should be a single TRUE or FALSE.
#'
#' @return \code{estMC} returns a list with elements:
#' \describe{
#'   \item{\code{sampleMC}}{\code{nSamples} sampled values for
#'      MC. Provided to allow the user to compute own summary statistics.}
#'   \item{\code{samplePsi}}{Array of sampled values for psi. \code{nSamples} x
#'      [number of origin sites] x [number of target sites]. Provided to allow
#'      the user to compute own summary statistics.}
#'   \item{\code{pointPsi}}{Simple point estimate of psi matrix, not accounting
#'      for sampling error. NULL when \code{isIntrinsic == TRUE}.}
#'   \item{\code{pointMC}}{Simple point estimate of MC, using the point
#'      estimates of \code{psi} and \code{originRelAbund}, not accounting
#'      for sampling error. NULL when \code{isIntrinsic == TRUE}.}
#'   \item{\code{meanMC, medianMC}}{Mean and median of \code{sampleMC}.
#'      Estimates of MC incorporating parametric uncertainty.}
#'   \item{\code{seMC}}{Standard error of MC, estimated from SD of
#'      \code{sampleMC}.}
#'   \item{\code{simpleCI}}{\code{1 - alpha} confidence interval for MC,
#'      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'      \code{sampleMC}.}
#'   \item{\code{bcCI}}{Bias-corrected \code{1 - alpha} confidence interval
#'      for MC.  Preferable to \code{simpleCI} when \code{meanMC} is the
#'      best estimate of MC. \code{simpleCI} is preferred when
#'      \code{medianMC} is a better estimator. When \code{meanMC==medianMC},
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of \code{sampleMC},
#'      where z0 is the proportion of \code{sampleMC < meanMC}.}
#'   \item{\code{hpdCI}}{\code{1 - alpha} credible interval for MC,
#'      estimated using the highest posterior density (HPD) method.}
#'   \item{\code{simpleP}}{Approximate p-value for MC, estimated as the
#'      proportion of bootstrap iterations where MC < \code{sigConst} (or MC >
#'      \code{sigConst} if \code{pointMC < sigConst}).  Note that if the
#'      proportion is 0, a default value of 0.5 / \code{nSamples} is provided,
#'      but this is best interpreted as p < 1 / \code{nSamples}.  NULL when
#'      \code{approxSigTest==FALSE}.}
#'   \item{\code{bcP}}{Approximate bias-corrected p-value for MC, estimated as
#'      \code{pnorm(qnorm(simpleP) - 2 * z0)}, where z0 is the proportion of
#'      \code{sampleMC < meanMC}.  May be a better approximation of the p-value
#'      than \code{simpleP}, but many of the same limitations apply.  NULL when
#'      \code{approxSigTest==FALSE}.}
#'   \item{\code{sampleCorr}}{\code{nBoot} sampled values for continuous
#'      correlation. Provided to allow the user to compute own summary
#'      statistics.  NULL when \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#'   \item{\code{pointCorr}}{Simple point estimate of continuous correlation,
#'      using \code{originPoints} and \code{targetPoints}, not accounting
#'      for sampling error. NULL when \code{calcCorr==FALSE} or
#'      \code{!is.null(psi)} or \code{isIntrinsic == TRUE}.}
#'   \item{\code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr}}{Summary
#'      statistics for continuous correlation bootstraps.  NULL when
#'      \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#'   \item{\code{inputSampleSize}}{If \code{sampleSize} was provided, this is
#'      that.  If not, it is either the calculated sample size (if that can be
#'      done), or left at NULL.  Useful to determining whether calcMC used the
#'      main MC formula or that for MC(R).}
#' }
#' @example inst/examples/estMCExamples.R
#' @seealso \code{\link{calcMC}}, \code{\link{projections}}, \code{\link{isoAssign}}
#' @export
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513 - 524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
#'
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. In revision. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of Mexico.

estMC <- function(originDist, targetDist = NULL, originRelAbund, psi = NULL,
                  sampleSize = NULL,
                  originSites = NULL, targetSites = NULL,
                  originPoints = NULL, targetPoints = NULL,
                  originAssignment = NULL, targetAssignment = NULL,
                  originNames = NULL, targetNames = NULL,
                  nSamples = 1000, nSim = ifelse(isTRUE(isIntrinsic), 10, 1000), isGL = FALSE,
                  geoBias = NULL, geoVCov = NULL, row0 = 0,
                  verbose = 0,  calcCorr = FALSE, alpha = 0.05,
                  approxSigTest = FALSE, sigConst = 0,
            resampleProjection = MigConnectivity::projections$EquidistConic,
                  maxTries = 300, targetIntrinsic = NULL,
                  isIntrinsic = FALSE) {
  if (is.null(psi)) {
    if(isIntrinsic) {
      mc <- estMCisotope(targetDist = targetDist,
                         originRelAbund = originRelAbund,
                         targetIntrinsic = targetIntrinsic,
                         targetSites = targetSites, sampleSize = sampleSize,
                         # targetAssignment = targetAssignment,
                         originPoints=originPoints, originSites=originSites,
                         originDist = originDist,
                         originAssignment = originAssignment,
                         originNames = originNames, targetNames = targetNames,
                         nBoot = nSamples, verbose = verbose, nSim = nSim,
                         calcCorr=calcCorr, alpha = alpha,
                         approxSigTest = approxSigTest, sigConst = sigConst,
                         maxTries = maxTries)
    }
    else {
      mc <- estMCGlGps(isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                       originRelAbund=originRelAbund, sampleSize = sampleSize,
                       originDist=originDist,
                       targetDist=targetDist,
                       targetPoints=targetPoints, targetSites=targetSites,
                       targetAssignment=targetAssignment,
                       originPoints=originPoints, originSites=originSites,
                       originAssignment=originAssignment,
                       originNames=originNames, targetNames=targetNames,
                       nBoot = nSamples, verbose=verbose,
                       nSim = nSim, calcCorr=calcCorr, alpha = alpha,
                       approxSigTest = approxSigTest, sigConst = sigConst,
                       resampleProjection = resampleProjection,
                       maxTries = maxTries)
    }
  }
  else {
    mc <- estMCCmrAbund(originRelAbund = originRelAbund,
                        sampleSize = sampleSize, psi = psi,
                        originDist = originDist, targetDist = targetDist,
                        originSites=originSites, targetSites=targetSites,
                        nSamples = nSamples, row0 = row0, verbose=verbose,
                        alpha = alpha, approxSigTest = approxSigTest,
                        sigConst = sigConst)
  }
  if (calcCorr & is.null(psi))
    class(mc) <- c("estMC", "estMantel", "estMigConnectivity")
  else
    class(mc) <- c("estMC", "estMigConnectivity")
  return(mc)
}

#' Estimate Mantel correlation (rM) from geolocator and/or GPS data.
#'
#' Resampling of uncertainty for rM from SpatialPoints geolocators and/or GPS
#' data.
#'
#' @param targetPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the point estimate location in
#'    the non-release season.
#' @param originPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the release location of an animal.
#' @param isGL Indicates whether or which animals were tracked with geolocators.
#'    Should be either single TRUE or FALSE value, or vector with length of
#'    number of animals tracked, with TRUE for animals in
#'    \code{targetPoints} with geolocators and FALSE for animals with GPS.
#' @param geoBias For GL data, vector of length 2 indicating expected bias
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters).
#' @param geoVCov For GL data, 2x2 matrix with expected variance/covariance
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters).
#' @param targetSites A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#'    object indicating valid target location(s).  Not needed unless you want
#'    to mask out certain areas (e.g. water).
#' @param nBoot Number of bootstrap runs. Animals are sampled with replacement for each,
#'    to estimate sampling uncertainty.
#' @param nSim Tuning parameter for GL data. Affects only the speed; 1000 seems
#'    to work well with our data.  Should be integer > 0.
#' @param verbose 0 (default) to 3. 0 prints no output during run. 1 prints
#'    a line every 100 bootstraps.  2 prints a line every bootstrap.
#'    3 also prints the number of draws (for tuning nSim for GL data only).
#' @param alpha Level for confidence/credible intervals provided.
#' @param resampleProjection Projection when sampling from geolocator
#'    bias/error. This projection needs units = m. Default is Equidistant
#'    Conic. The default setting preserves distances around latitude = 0 and
#'    longitude = 0. Other projections may work well, depending on the location
#'    of \code{targetSites}.
#' @param maxTries Maximum number of times to run a single GL bootstrap before
#'    exiting with an error.  Default is 300.  Set to NULL to never stop.  This
#'    parameter was added to prevent GL setups where some sample points never
#'    land on target sites from running indefinitely.
#'
#' @return \code{estMantel} returns a list with elements:
#' \describe{
#'   \item{\code{sampleCorr}}{\code{nBoot} sampled values for Mantel
#'      correlation. Provided to allow the user to compute own summary
#'      statistics.}
#'   \item{\code{pointCorr}}{Simple point estimate of Mantel correlation,
#'      using \code{originPoints} and \code{targetPoints}.}
#'   \item{\code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr}}{Summary
#'      statistics for Mantel correlation bootstraps.}
#' }
#' @export
#'
#' @examples
#' rM1 <- estMantel(isGL=OVENdata$isGL,#Logical vector: light-level GL(T)/GPS(F)
#'                  geoBias = OVENdata$geo.bias, # Geolocator location bias
#'                  geoVCov = OVENdata$geo.vcov, # Location covariance matrix
#'                  targetSites = OVENdata$targetSites, # Non-breeding target sites
#'                  originPoints = OVENdata$originPoints, # Capture Locations
#'                  targetPoints = OVENdata$targetPoints, # Target locations
#'                  verbose = 1,   # output options
#'                  nBoot = 100, # This is set low for example
#'                  resampleProjection = raster::projection(OVENdata$targetSites))
#' str(rM1)
#' @seealso \code{\link{estMC}}
#'
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513 - 524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}

estMantel <- function(targetPoints, originPoints, isGL, geoBias = NULL,
                      geoVCov = NULL, targetSites = NULL, nBoot = 1000,
                      nSim = 1000, verbose=0, alpha = 0.05,
            resampleProjection = MigConnectivity::projections$EquidistConic,
            maxTries = 300) {

  # Input checking and assignment
  if (!(verbose %in% 0:3))
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = every animal")
  if (length(geoBias)!=2 && any(isGL))
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in meters)")
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = F)) && any(isGL))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in meters)")
  nAnimals <- length(targetPoints)
  if (length(isGL)==1)
    isGL <- rep(isGL, nAnimals)
  if(is.na(raster::projection(originPoints))) {
    stop('Coordinate system definition needed for originPoints')
  }
  if(is.na(raster::projection(targetPoints))){
    stop('Coordinate system definition needed for targetPoints')
  }
  targetPoints <- sp::spTransform(targetPoints, sp::CRS(resampleProjection))
  if(!is.null(targetSites)){
    if(class(targetSites)=="SpatialPolygonsDataFrame"){
      targetSites <- sp::SpatialPolygons(targetSites@polygons)}
    if(is.na(raster::projection(targetSites))){
      stop('Coordinate system definition needed for targetSites')
    }
    targetSites <- sp::spTransform(targetSites, sp::CRS(resampleProjection))
  }

  corr <- rep(NA, nBoot)

  targetDist1 <- matrix(NA, nAnimals, nAnimals)

  targetDist1[lower.tri(targetDist1)] <- 1

  distIndices <- which(!is.na(targetDist1), arr.ind = TRUE)

  # project target points to WGS #
  targetPoints2 <- sp::spTransform(targetPoints, sp::CRS(MigConnectivity::projections$WGS84))

  targetDist0 <- geosphere::distVincentyEllipsoid(targetPoints2[distIndices[,'row'],],targetPoints2[distIndices[,'col'],])

  targetDist1[lower.tri(targetDist1)] <- targetDist0

  originPoints2 <- sp::spTransform(originPoints, sp::CRS(MigConnectivity::projections$WGS84))

  originDistStart <- matrix(geosphere::distVincentyEllipsoid(originPoints2[rep(1:nAnimals, nAnimals)], originPoints2[rep(1:nAnimals,each=nAnimals)]),
                            nAnimals, nAnimals)

  originDist1 <- originDistStart[1:nAnimals, 1:nAnimals]

  pointCorr <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = TRUE)$correlation

  for (boot in 1:nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Sample individual animals with replacement
    animal.sample <- sample.int(nAnimals, replace=TRUE)
    origin.point.sample <- originPoints[animal.sample]
    tSamp <- targetSample(isGL = isGL, geoBias = geoBias, geoVCov = geoVCov,
                          targetPoints = targetPoints, animal.sample = animal.sample,
                          targetSites = targetSites,
                          resampleProjection = resampleProjection,
                          maxTries = maxTries)
    target.sample <- tSamp$target.sample
    target.point.sample <- tSamp$target.point.sample
    if (verbose > 2)
      cat(' ', tSamp$draws, 'draws (of length', nSim, 'and of', maxTries, 'possible).\n')

    originDist1 <- originDistStart[animal.sample, animal.sample]
    target.point.sample <- sp::SpatialPoints(target.point.sample,sp::CRS(resampleProjection))
    target.point.sample2 <- sp::spTransform(target.point.sample,sp::CRS(MigConnectivity::projections$WGS84))
    targetDist0 <- geosphere::distVincentyEllipsoid(target.point.sample2[distIndices[,'row']],
                                                    target.point.sample2[distIndices[,'col']])
    targetDist1[lower.tri(targetDist1)] <- targetDist0
    corr[boot] <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = TRUE)$correlation
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" Correlation mean:", mean(corr, na.rm=TRUE), "SD:", sd(corr, na.rm=TRUE),
          "low quantile:", quantile(corr, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(corr, 1-alpha/2, na.rm=TRUE), "\n")
  }
  meanCorr <- mean(corr, na.rm=TRUE)
  medianCorr <- median(corr, na.rm=TRUE)
  seCorr <- sd(corr, na.rm=TRUE)
  simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                           names = F)
  corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
  bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, type = 8, names = F)
  return(structure(list(sampleCorr = corr, pointCorr = pointCorr,
                        meanCorr = meanCorr, medianCorr = medianCorr,
                        seCorr=seCorr, simpleCICorr=simpleCICorr,
                        bcCICorr=bcCICorr, alpha = alpha),
                   class = c("estMantel", "estMigConnectivity")))
}

###############################################################################
#' Grab (from https://github.com/SMBC-NZP/MigConnectivity) example RMark
#' transition probability estimates obtained from simulated data
#'
#' Get a dataset containing RMark transition probability estimates from
#' simulated mark-recapture-recovery data from Cohen et al. (2014).  These all
#' represent the intermediate scenario for all settings (moderate connectivity,
#' low re-encounter, 100,000 banded in each breeding area).  Each estimate can
#' be used in \code{estMC} function to estimate MC with uncertainty.  Requires
#' internet connection.
#'
#' @param number Integer 1 - 100, which simulation and RMark estimate you want
#'
#' @return RMark object
#' @export
#' @seealso \code{\link{estMC}}
getCMRexample <- function(number = 1) {
  obj.name <- paste0('psiB.enc2.band100.', number)
  file.name <- paste0('out_', obj.name, '.rds')
  url1 <- paste0('https://github.com/SMBC-NZP/MigConnectivity/blob/master/data-raw/', file.name, '?raw=true')
  temp <- paste(tempdir(), file.name, sep = '/')
  utils::download.file(url1, temp, mode = 'wb')
  fm <- readRDS(temp)
  unlink(temp)
  return(fm)
}


#' Pairwise differences between two or more independent MC estimates
#'
#' Estimates mean (and median) differences in MC, and includes measures of
#' uncertainty (SE and CI).  For those measures of uncertainty to be accurate,
#' only apply this function to MC estimates where all data sources are
#' independent (e.g., different species).
#'
#' @param estimates List of at leat two MC estimates, provided by the estMC
#'    function. If this is a named list (recommended), the function will use
#'    these names in labeling the differences.
#' @param nSamples A positive integer, number of samples (with replacement)
#'    to draw from each pair of MC estimates (default 100000).  If set to NULL,
#'    compares all MC samples from each pair.
#' @param alpha Level for confidence/credible intervals provided.
#' @param returnSamples Should the function return all the sampled differences?
#'    Defaults to FALSE to reduce storage requirements. Change to TRUE to
#'    compute your own summary statistics.
#'
#' @return  \code{diffMC} returns a list with elements:
#' \describe{
#'   \item{\code{meanDiff, medianDiff}}{Vectors with mean and medians of sampled
#'      differences for each pairwise comparison. Estimates of difference
#'      between MC values incorporating parametric uncertainty.}
#'   \item{\code{seDiff}}{Vector with standard errors of MC differences for each
#'      pairwise comparison, estimated from SD of sampled differences.}
#'   \item{\code{simpleCI}}{Matrix of \code{1 - alpha} confidence intervals for
#'      MC differences, estimated as \code{alpha/2} and \code{1 - alpha/2}
#'      quantiles of \code{sampleMC}.}
#'   \item{\code{bcCI}}{Matrix of bias-corrected \code{1 - alpha} confidence
#'      intervals for MC differences for each pairwise comparison. Preferable
#'      to \code{simpleCI} when \code{meanDiff} is the best estimate of the MC
#'      difference. \code{simpleCI} is preferred when
#'      \code{medianDiff} is a better estimator. When \code{meanDiff==medianDiff},
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of sampled
#'      differences, where z0 is the proportion of \code{sampleDiff < meanDiff}.}
#   \item{\code{hpdCI}}{Matrix of \code{1 - alpha} credible intervals for MC
#      differences for each pairwise comparison, etimated using the highest
#      posterior density (HPD) method.}
#'   \item{\code{sampleDiff}}{Only provided if \code{returnSamples} is TRUE.
#'      List of sampled values for each pairwise MC difference.}
#' }
#' @export
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. In revision. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of Mexico.
#'
# @examples
diffMC <- function(estimates, nSamples = 100000, alpha = 0.05, returnSamples = F) {
  nEst <- length(estimates)
  nComparisons <- choose(nEst, 2)
  nSamplesEst <- sapply(estimates, function(x) length(x$sampleMC))
  diffSamples <- vector('list', nComparisons)
  if (is.null(names(estimates)))
    names(estimates) <- 1:nEst
  comparisons <- matrix(c(sequence(1:(nEst - 1)), rep(2:nEst, 1:(nEst - 1))), nComparisons, 2)
  for (i in 1:nComparisons) {
    if (is.null(nSamples))
      diffSamples[[i]] <- rep(estimates[[comparisons[i, 1]]]$sampleMC, times = nSamplesEst[comparisons[i, 2]]) -
        rep(estimates[[comparisons[i, 2]]]$sampleMC, each = nSamplesEst[comparisons[i, 1]])
    else
      diffSamples[[i]] <- sample(estimates[[comparisons[i, 1]]]$sampleMC, nSamples, replace = T) -
        sample(estimates[[comparisons[i, 2]]]$sampleMC, nSamples, replace = T)
  }
  meanDiff <- sapply(diffSamples, mean, na.rm=TRUE)
  medianDiff <- sapply(diffSamples, median, na.rm=TRUE)
  seDiff <- sapply(diffSamples, sd, na.rm=TRUE)
  simpleCI <- sapply(diffSamples, quantile, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                     type = 8, names = F)
  diff.z0 <- sapply(diffSamples, function(MC) qnorm(sum(MC<mean(MC, na.rm = T), na.rm = T)/length(which(!is.na(MC)))))
  bcCI <- mapply(function(MC, z0) quantile(MC, pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, type = 8, names = F), diffSamples, diff.z0)
  diff.mcmc <- lapply(diffSamples, coda::as.mcmc)
  #hpdCI <- sapply(diff.mcmc, function(MC) as.vector(coda::HPDinterval(MC, 1-alpha)))
  names(diffSamples) <- names(meanDiff) <- paste(names(estimates[comparisons[,1]]),
                                                 '-', names(estimates[comparisons[,2]]))
  names(medianDiff) <- names(seDiff) <- names(diffSamples)
  colnames(simpleCI) <- colnames(bcCI) <- names(diffSamples) #colnames(hpdCI) <-
  sampleDiff <- ifelse(returnSamples, diffSamples, NA)
  return(structure(list(meanDiff = meanDiff, medianDiff = medianDiff,
                        seDiff = seDiff, simpleCI = simpleCI, bcCI = bcCI,
                        sampleDiff = sampleDiff, alpha = alpha),#hpdCI = hpdCI,
                   class = c('diffMC', 'diffMigConnectivity')))
}

#' Pairwise differences between two or more independent Mantel correlation
#' estimates
#'
#' Estimates mean (and median) differences in Mantel correlations (rM), and
#' includes measures of uncertainty (SE and CI).  For those measures of
#' uncertainty to be accurate, only apply this function to rM estimates where
#' all data sources are independent (e.g., different species).
#'
#' @param estimates List of at leat two Mantel correlation estimates, provided
#'    by either the estMC or the estMantel functions. If this is a named list
#'    (recommended), the function will use these names in labeling the
#'    differences.
#' @param nSamples A positive integer, number of samples (with replacement)
#'    to draw from each pair of MC estimates (default 100000).  If set to NULL,
#'    compares all Mantel correlation samples from each pair.
#' @param alpha Level for confidence/credible intervals provided.
#' @param returnSamples Should the function return all the sampled differences?
#'    Defaults to FALSE to reduce storage requirements. Change to TRUE to
#'    compute your own summary statistics.
#'
#' @return  \code{diffMantel} returns a list with elements:
#' \describe{
#'   \item{\code{meanDiff, medianDiff}}{Vectors with mean and medians of sampled
#'      differences for each pairwise comparison. Estimates of difference
#'      between rM values incorporating parametric uncertainty.}
#'   \item{\code{seDiff}}{Vector with standard errors of rM differences for each
#'      pairwise comparison, estimated from SD of sampled differences.}
#'   \item{\code{simpleCI}}{Matrix of \code{1 - alpha} confidence intervals for
#'      rM differences, estimated as \code{alpha/2} and \code{1 - alpha/2}
#'      quantiles of \code{sampleCorr}.}
#'   \item{\code{bcCI}}{Matrix of bias-corrected \code{1 - alpha} confidence
#'      intervals for rM differences for each pairwise comparison. Preferable
#'      to \code{simpleCI} when \code{meanDiff} is the best estimate of the rM
#'      difference. \code{simpleCI} is preferred when
#'      \code{medianDiff} is a better estimator. When \code{meanDiff==medianDiff},
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of sampled
#'      differences, where z0 is the proportion of \code{sampleDiff < meanDiff}.}
#   \item{\code{hpdCI}}{Matrix of \code{1 - alpha} credible intervals for rM
#      differences for each pairwise comparison, etimated using the highest
#      posterior density (HPD) method.}
#'   \item{\code{sampleDiff}}{Only provided if \code{returnSamples} is TRUE.
#'      List of sampled values for each pairwise rM difference.}
#' }
#' @export
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. In revision. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of Mexico.
#'
# @examples
diffMantel <- function(estimates, nSamples = 100000, alpha = 0.05, returnSamples = F) {
  nEst <- length(estimates)
  nComparisons <- choose(nEst, 2)
  nSamplesEst <- sapply(estimates, function(x) length(x$sampleCorr))
  diffSamples <- vector('list', nComparisons)
  if (is.null(names(estimates)))
    names(estimates) <- 1:nEst
  comparisons <- matrix(c(sequence(1:(nEst - 1)), rep(2:nEst, 1:(nEst - 1))), nComparisons, 2)
  for (i in 1:nComparisons) {
    if (is.null(nSamples))
      diffSamples[[i]] <- rep(estimates[[comparisons[i, 1]]]$sampleCorr, times = nSamplesEst[comparisons[i, 2]]) -
        rep(estimates[[comparisons[i, 2]]]$sampleCorr, each = nSamplesEst[comparisons[i, 1]])
    else
      diffSamples[[i]] <- sample(estimates[[comparisons[i, 1]]]$sampleCorr, nSamples, replace = T) -
        sample(estimates[[comparisons[i, 2]]]$sampleCorr, nSamples, replace = T)
  }
  meanDiff <- sapply(diffSamples, mean, na.rm=TRUE)
  medianDiff <- sapply(diffSamples, median, na.rm=TRUE)
  seDiff <- sapply(diffSamples, sd, na.rm=TRUE)
  simpleCI <- sapply(diffSamples, quantile, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                     type = 8, names = F)
  diff.z0 <- sapply(diffSamples, function(MC) qnorm(sum(MC<mean(MC, na.rm = T), na.rm = T)/length(which(!is.na(MC)))))
  bcCI <- mapply(function(MC, z0) quantile(MC, pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, type = 8, names = F), diffSamples, diff.z0)
  diff.mcmc <- lapply(diffSamples, coda::as.mcmc)
  #hpdCI <- sapply(diff.mcmc, function(MC) as.vector(coda::HPDinterval(MC, 1-alpha)))
  names(diffSamples) <- names(meanDiff) <- paste(names(estimates[comparisons[,1]]),
                                                 '-', names(estimates[comparisons[,2]]))
  names(medianDiff) <- names(seDiff) <- names(diffSamples)
  colnames(simpleCI) <- colnames(bcCI) <- names(diffSamples) #colnames(hpdCI) <-
  sampleDiff <- ifelse(returnSamples, diffSamples, NA)
  return(structure(list(meanDiff = meanDiff, medianDiff = medianDiff,
                        seDiff = seDiff, simpleCI = simpleCI, bcCI = bcCI,
                        sampleDiff = sampleDiff, alpha = alpha),
                   class = c('diffMantel', 'diffMigConnectivity')))
}
