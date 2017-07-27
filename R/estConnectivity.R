###############################################################################
# Estimate MC from abundance and transition probability estimates.
#
# Resampling of uncertainty for MC from RMark psi matrix estimates and/or JAGS
# relative abundance MCMC samples.
#
# @param originRelAbund Relative abundances at B origin sites. Either a
#  numeric vector of length B that sums to 1 or an mcmc object with
#  \code{nSamples} rows  and columns including 'relN[1]' through 'relN[B]'.
# @param psi Transition probabilities between B origin and W target sites.
#  Either a matrix with B rows and W columns where rows sum to 1 or a MARK
#  object with estimates of transition probabilities.
# @param originDist Distances between the B origin sites.  Symmetric B by B
#  matrix.
# @param targetDist Distances between the W target sites.  Symmetric W by W
#  matrix.
# @param originSites If \code{psi} is a MARK object, this must be a numeric
#  vector indicating which sites are origin.
# @param targetSites If \code{psi} is a MARK object, this must be a numeric
#  vector indicating which sites are target.
# @param nSamples Number of times to resample \code{psi} and/or
#  \code{originRelAbund}.
# @param row0 If \code{originRelAbund} is an mcmc object, this can be set
#  to 0 (default) or any greater integer to specify where to stop ignoring
#  samples ("burn-in").
# @param verbose 0 (default), 1, or 2. 0 prints no output during run. 1 prints
#  a line every 100 samples.  2 prints a line every sample.
# @param alpha Level for confidence intervals provided.
# @return \code{estMCCmrAbund} returns a list with elements:
# \describe{
#   \item{\code{sampleMC}}{\code{nSamples} sampled values for MC. Provided
#      to allow the user to compute own summary statistics.}
#   \item{\code{pointMC}}{Simple point estimate of MC, using the point
#      estimates of \code{psi} and \code{originRelAbund}.}
#   \item{\code{meanMC, medianMC}}{Mean and median of \code{sampleMC}.
#      Estimates of MC incorporating parametric uncertainty.}
#   \item{\code{seMC}}{Standard error of MC, estimated from SD of
#      \code{sampleMC}.}
#   \item{\code{simpleCI}}{\code{1 - alpha} confidence interval for MC,
#      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#      \code{sampleMC}.}
#   \item{\code{bcCI}}{Bias-corrected \code{1 - alpha} confidence interval
#      for MC.  Preferable to \code{simpleCI} when \code{pointMC} is the
#      best estimate of MC. \code{simpleCI} is preferred when
#      \code{medianMC} is a better estimator. When \code{pointMC==medianMC},
#      these should be equivalent.}
#   \item{\code{hpdCI}}{\code{1 - alpha} credible interval for MC,
#      estimated using the highest posterior density (HPD) method.}
# }
# @example inst/examples/estMCCmrAbundExamples.R
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
                                                       se=T))[originSites,
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
      cat("\tSample",i,"of",nSamples,"\n")
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
  }
  meanMC <- mean(sampleMC)
  medianMC <- median(sampleMC)
  seMC <- sd(sampleMC)
  # Calculate confidence intervals using quantiles of sampled MC
  simpleCI <- quantile(sampleMC, c(alpha/2, 1-alpha/2), na.rm=T, type = 8)
  z0 <- qnorm(sum((sampleMC)<meanMC)/nSamples)
  bcCI <- quantile(sampleMC, pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=T, type = 8)
  MC.mcmc <- coda::as.mcmc(sampleMC) # Ha!
  hpdCI <- coda::HPDinterval(MC.mcmc, 1-alpha)
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (pointMC > sigConst)
      simpleP <- sum(sampleMC < sigConst) / nBoot
    else
      simpleP <- sum(sampleMC > sigConst) / nBoot
    if (simpleP == 0)
      simpleP <- 0.5 / nBoot
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
              simpleCICorr=NULL, bcCICorr=NULL, inputSampleSize = sampleSize))
}

###############################################################################
# Resampling of uncertainty from geolocators and/or GPS data
#
# @param isGL Vector indicating which animals were tracked with geolocators.
#    Should be length of number of animals tracked, with TRUE for animals in
#    \code{targetPoints} with geolocators and FALSE for animals with GPS.
# @param geoBias Vector of length 2 indicating expected bias in longitude and
#    latitude of \code{targetPoints}, in meters.
# @param geoVCov 2x2 matrix with expected variance/covariance in longitude and
#    latitude of \code{targetPoints}, in meters.
# @param originRelAbund Relative abundances at B origin sites. A numeric
#    vector of length B that sums to 1.
# @param originDist Distances between the B origin sites.  Symmetric B by B
#    matrix.
# @param targetDist Distances between the W target sites.  Symmetric W by W
#    matrix.
# @param targetPoints A \code{SpatialPoints} object, with length number of
#    animals tracked.  Each point indicates the point estimate location in
#    the non-release season.
# @param targetSites Geographic definition of sites in the non-release season.
# @param targetAssignment Optional. Point estimate assignment of
#    \code{targetPoints} to non-release season sites. Integer vector with
#    length number of animals tracked.
# @param originPoints A \code{SpatialPoints} object, with length number of
#    animals tracked.  Each point indicates the release location of an animal.
# @param originSites Geographic definition of sites in the release season.
# @param originAssignment Assignment of \code{originPoints} to release season
#    sites. Integer vector with length number of animals tracked. Optional,
#    but either \code{originAssignment} or \code{originSites} and
#    \code{originPoints} should be defined.
# @param originNames Optional. Vector of names for the release season sites.
# @param targetNames Optional. Vector of names for the non-release season
#    sites.
# @param nBoot Number of bootstrap runs. Animals are sampled with replacement
#    for each of these to estimate sampling uncertainty.
# @param verbose Integer 0-3 for level of output during bootstrap: 0 = none,
#    1 = every 10, 2 = every run, 3 = every animal.
# @param nSim Number of times to sample random points for each animal from
#    parametric distribution of non-release season error. Ignored for GPS
#    points (assumed to have no geographic error).
# @param calcCorr In addition to MC, should function also estimate continuous
#    correlation between release and non-release locations?  Default is TRUE.
# @param alpha Level for confidence intervals provided.
# @param projection.dist.calc Projection when sampling from geolocator bias/error.
#       Default Equidistant Conic = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
#
# @return A list with elements:
# \describe{
#   \item{\code{sampleMC}}{\code{nBoot} sampled values for MC. Provided
#      to allow the user to compute own summary statistics.}
#   \item{\code{samplePsi}}{Array of sampled values for psi. \code{nBoot} x
#      [number of origin sites] x [number of target sites]. Provided
#      to allow the user to compute own summary statistics.}
#   \item{\code{pointSites}}{Matrix of point assignment of number of animals
#      to each origin and target site combination.}
#   \item{\code{pointPsi}}{Simple point estimate of psi matrix, as
#      \code{prop.table(pointSites, 1)}.}
#   \item{\code{pointMC}}{Simple point estimate of MC, using \code{pointPsi}
#      and \code{originRelAbund}.}
#   \item{\code{meanMC, medianMC}}{Mean and median of \code{sampleMC}.
#      Estimates of MC incorporating parametric uncertainty.}
#   \item{\code{seMC}}{Standard error of MC, estimated from SD of
#      \code{sampleMC}.}
#   \item{\code{simpleCI}}{\code{1 - alpha} confidence interval for MC,
#      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#      \code{sampleMC}.}
#   \item{\code{bcCI}}{Bias-corrected \code{1 - alpha} confidence interval
#      for MC.  Preferable to \code{simpleCI} when \code{pointMC} is the
#      best estimate of MC. \code{simpleCI} is preferred when
#      \code{medianMC} is a better estimator. When \code{pointMC==medianMC},
#      these should be equivalent.}
#   \item{\code{sampleCorr}}{\code{nBoot} sampled values for continuous
#      correlation. Provided to allow the user to compute own summary
#      statistics.}
#   \item{\code{pointCorr}}{Simple point estimate of continuous correlation,
#      using \code{originPoints} and \code{targetPoints}.}
#   \item{\code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr}}{Summary
#      statistics for continuous correlation bootstraps.}
# }
#
# @examples
estMCGlGps <- function(originDist, targetDist, originRelAbund, isGL,
                       geoBias, geoVCov, targetPoints, targetSites,
                       sampleSize = NULL,
                       targetAssignment=NULL,
                       originPoints=NULL, originSites=NULL,
                       originAssignment=NULL, originNames=NULL,
                       targetNames=NULL, nBoot = 1000, verbose=0,
                       nSim = 1000, calcCorr=T, alpha = 0.05,
                       approxSigTest = F, sigConst = 0,
                       projection.dist.calc = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs") {

  # Input checking and assignment
  if (!(verbose %in% 0:3))
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = every animal")
  if (length(geoBias)!=2 && (any(isGL) || calcCorr))
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in meters)")
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = F)) && (any(isGL) || calcCorr))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in meters)")
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
    originSites <- sp::SpatialPolygons(originSites@polygons)}
  if(class(targetSites)=="SpatialPolygonsDataFrame"){
    targetSites <- sp::SpatialPolygons(targetSites@polygons)}
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
  if(!is.null(originPoints) && is.na(originPoints@proj4string)){
    stop('Coordinate system definition needed for originSites')
  }
  if(is.na(raster::projection(targetSites)) || is.na(raster::projection(targetPoints))){
    stop('Coordinate system definition needed for targetSites & targetPoints')
  }

  if (is.null(sampleSize))
    sampleSize <- nAnimals
  if (dim(originDist)!=rep(nOriginSites,2) ||
      dim(targetDist)!=rep(nTargetSites,2))
    stop('Distance matrices should be square with same number of sites of each type as assignments/points (with distances in meters)')
  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))


  WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

  Dist.proj <- projection.dist.calc

  MC <- corr <- rep(NA, nBoot)

#   if (any(isGL) || calcCorr)
#
#     geoBias2 <- array(rep(geoBias, nAnimals, each = nSim), c(nSim, 2, nAnimals))

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

    distIndices <- which(!is.na(targetDist1), arr.ind = T)

    # project target points to WGS #
    targetPoints2 <- sp::spTransform(targetPoints, sp::CRS(WGS84))

    targetDist0 <- geosphere::distVincentyEllipsoid(targetPoints2[distIndices[,'row'],],targetPoints2[distIndices[,'col'],])

    targetDist1[lower.tri(targetDist1)] <- targetDist0

    originPoints2 <- sp::spTransform(originPoints, sp::CRS(WGS84))

    originDistStart <- matrix(geosphere::distVincentyEllipsoid(originPoints2[rep(1:nAnimals, nAnimals)], originPoints2[rep(1:nAnimals,each=nAnimals)]),
                                    nAnimals, nAnimals)

    originDist1 <- originDistStart[1:nAnimals, 1:nAnimals]

    pointCorr <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = T)$correlation
    #ade4::mantel.rtest(as.dist(originDist1), as.dist(targetDist1),
                 #              nrepet=0)
  }

  for (boot in 1:nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have birds from every origin site
    origin.sample <- 'Filler' # Start with one origin site
    while (length(unique(origin.sample)) < nOriginSites) { #2
      # Sample individual animals with replacement
      animal.sample <- sample.int(nAnimals, replace=T)
      # Get origin points for those animals
      if (calcCorr)
        origin.point.sample <- originPoints[animal.sample]
      # Get origin population for each animal sampled
      origin.sample <- originAssignment[animal.sample]
    }
    target.sample <- rep(NA, nAnimals)
    target.point.sample <- matrix(NA, nAnimals, 2)
    toSample <- which(isGL[animal.sample])
    draws <- 0
    while (length(toSample) > 0) {
      draws <- draws + 1
      geoBias2 <- array(rep(geoBias, length(toSample), each = nSim), c(nSim, 2, length(toSample)))
      point.sample <- array(apply(targetPoints@coords[animal.sample[toSample], , drop = FALSE], 1,
                            MASS::mvrnorm, n=nSim, Sigma=geoVCov),
                            c(nSim, 2, length(toSample))) + geoBias2
      point.sample <- apply(point.sample, 3, sp::SpatialPoints,
                            proj4string = sp::CRS(projection.dist.calc))
      target.sample0 <- sapply(point.sample, sp::over, y = targetSites)
      good.sample <- apply(target.sample0, 2, function(x) which(!is.na(x))[1])
      target.sample[toSample] <- apply(target.sample0, 2, function(x) x[!is.na(x)][1])
      if (any(!is.na(good.sample)))
        target.point.sample[toSample[!is.na(good.sample)], ]<- t(mapply(function(x, y) y[x]@coords,
                                               good.sample[!is.na(good.sample)], point.sample[!is.na(good.sample)]))
      toSample <- which(isGL[animal.sample] & is.na(target.sample))
    }
    target.sample[which(!isGL[animal.sample])] <- targetAssignment[animal.sample[which(!isGL[animal.sample])]]
    if (calcCorr)
      target.point.sample[which(!isGL[animal.sample]), ] <- targetPoints[animal.sample[which(!isGL[animal.sample])],]@coords
#     for(i in 1:nAnimals){
#       if (verbose > 2)
#         cat('\tAnimal', i, 'of', nAnimals)
#       draws <- 0
#       if (isGL[animal.sample[i]]) {
#         while (is.na(target.sample[i])) {
#           draws <- draws + 1
#           # Sample random point for each bird from parametric distribution of NB error
#           point.sample <- sp::SpatialPoints(MASS::mvrnorm(n=nSim, mu=cbind(
#             targetPoints@coords[animal.sample[i],1],
#             targetPoints@coords[animal.sample[i],2]), Sigma=geoVCov)+
#               geoBias2, sp::CRS(projection.dist.calc))
#           # filtered to stay in NB areas (land)
#           target.sample0 <- sp::over(point.sample, targetSites)
#           target.sample[i]<-target.sample0[!is.na(target.sample0)][1]
#         }
#         target.point.sample[i, ]<-point.sample[!is.na(target.sample0)][1]@coords
#       }
#       else { # Assume no location error for GPS
#         target.sample[i] <- targetAssignment[animal.sample[i]]
#         if (calcCorr)
#           target.point.sample[i, ] <- targetPoints[animal.sample[i],]@coords
#       }
#       if (verbose > 2)
#         cat(' ', draws, 'draws\n')
      # Now that we have breeding and non-breeding site for point, add to transition count matrix
    sites <- table(origin.sample, target.sample)
    sites.array[boot, as.integer(rownames(sites)), as.integer(colnames(sites))] <- sites
    # Create psi matrix as proportion of those from each breeding site that went to each NB site
    psi.array[boot, , ] <- prop.table(sites.array[boot, , ], 1)
    # Calculate MC from that psi matrix
    MC[boot] <- calcMC(originDist, targetDist, originRelAbund,
                       psi.array[boot, , ], sampleSize)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=T), "SD:", sd(MC, na.rm=T),
          "low quantile:", quantile(MC, alpha/2, na.rm=T),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=T), "\n")
    if (calcCorr) {
      originDist1 <- originDistStart[animal.sample, animal.sample]
      target.point.sample <- sp::SpatialPoints(target.point.sample,sp::CRS(projection.dist.calc))
      target.point.sample2 <- sp::spTransform(target.point.sample,sp::CRS(WGS84))
      targetDist0 <- geosphere::distVincentyEllipsoid(target.point.sample2[distIndices[,'row']],
                                           target.point.sample2[distIndices[,'col']])
      targetDist1[lower.tri(targetDist1)] <- targetDist0
      corr[boot] <- ncf::mantel.test(originDist1, targetDist1, resamp=0, quiet = T)$correlation
        #ade4::mantel.rtest(as.dist(originDist1),as.dist(targetDist1),
        #                         nrepet=0)
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=T), "SD:", sd(corr, na.rm=T),
            "low quantile:", quantile(corr, alpha/2, na.rm=T),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=T), "\n")
    }
  }
  MC.z0 <- qnorm(sum((MC)<mean(MC))/nBoot)
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=T, type = 8)
  MC.mcmc <- coda::as.mcmc(MC) # Ha!
  hpdCI <- coda::HPDinterval(MC.mcmc, 1-alpha)
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
    meanCorr <- mean(corr)
    medianCorr <- median(corr)
    seCorr <- sd(corr)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=T, type = 8)
    corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                           na.rm=T, type = 8)
  } else
    pointCorr <- meanCorr <- medianCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  return(list(sampleMC = MC, samplePsi = psi.array,
              pointPsi = pointPsi, pointMC = pointMC, meanMC = mean(MC),
              medianMC = median(MC), seMC = sd(MC),
              simpleCI = quantile(MC, c(alpha/2, 1-alpha/2), na.rm=T, type = 8),
              bcCI = bcCI, hpdCI = hpdCI, simpleP = simpleP, bcP = bcP,
              sampleCorr = corr, pointCorr = pointCorr,
              meanCorr = meanCorr, medianCorr = medianCorr, seCorr=seCorr,
              simpleCICorr=simpleCICorr, bcCICorr=bcCICorr,
              inputSampleSize = sampleSize))
}

###############################################################################
#' Estimate MC from abundance and/or transition probability estimates OR
#' geolocator and/or GPS data.
#'
#' Resampling of uncertainty for MC from RMark psi matrix estimates and/or JAGS
#' relative abundance MCMC samples OR SpatialPoints geolocators and/or GPS
#' data.
#'
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'  matrix.
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'  matrix.
#' @param originRelAbund Relative abundance estimates at B origin sites. Either
#'  a numeric vector of length B that sums to 1 or an mcmc object with
#'  \code{nSamples} rows  and columns including 'relN[1]' through 'relN[B]'.
#'  Currently, an mcmc object doesn't work with geolocator or GPS data.
#' @param psi Transition probabilities between B origin and W target sites.
#'  Either a matrix with B rows and W columns where rows sum to 1 or a MARK
#'  object with estimates of transition probabilities.  If you are estimating
#'  MC from GPS or geolocator data, leave this as NULL.
#' @param sampleSize Total sample size of animals that psi will be estimated from.
#'    Should be the number of animals released in one of the origin sites and
#'    observed in one of the target sites.  Optional, but recommended, unless
#'    you are estimating MC from GPS or geolocator data (in which case the
#'    function can calculate it for you).
#' @param originSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are origin.  If using GL or GPS data,
#'  this can be the geographic definition of sites in the release season.
#' @param targetSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are target.  If using GL or GPS data,
#'  this must be the geographic definition of sites in the non-release season.
#' @param originPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the release location of an animal.
#' @param targetPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the point estimate location in
#'    the non-release season.
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
#'    \code{originRelAbund} OR number of bootstrap runs for GL or GPS data. In
#'    the latter case, animals are sampled with replacement for each. For all,
#'    the purpose is to estimate sampling uncertainty.
#' @param nSim Tuning parameter for GL data. Affects only the speed; 1000 seems
#'    to work well with our data.  Should be integer > 0.
#' @param isGL Indicates whether or which animals were tracked with geolocators.
#'    Should be either single TRUE or FALSE value, or vector with length of
#'    number of animals tracked, with TRUE for animals in
#'    \code{targetPoints} with geolocators and FALSE for animals with GPS.
#' @param geoBias For GL data, vector of length 2 indicating expected bias
#'    in longitude and latitude of \code{targetPoints}, in meters.
#' @param geoVCov For GL data, 2x2 matrix with expected variance/covariance
#'    in longitude and latitude of \code{targetPoints}, in meters.
#' @param row0 If \code{originRelAbund} is an mcmc object, this can be set
#'  to 0 (default) or any greater integer to specify where to stop ignoring
#'  samples ("burn-in").
#' @param verbose 0 (default) to 3. 0 prints no output during run. 1 prints
#'  a line every 100 samples or bootstraps.  2 prints a line every sample or
#'  bootstrap. 3 also prints a line every animal (GL or GPS data only).
#' @param calcCorr In addition to MC, should function also estimate continuous
#'    correlation between release and non-release locations (GPS or GL data
#'    only)?  Default is FALSE.
#' @param alpha Level for confidence/credible intervals provided.
#' @param approxSigTest Should function compute approximate one-sided
#'    significance tests (p-values) for MC from the bootstrap?  Default is
#'    FALSE.
#' @param sigConst Value to compare MC to in significance test.
#'    Default is 0.
#' @param projection.dist.calc Projection when sampling from geolocator bias/error.
#'    Default Equidistant Conic = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs" 
#'
#' @return \code{estMC} returns a list with elements:
#' \describe{
#'   \item{\code{sampleMC}}{\code{nSamples} sampled values for
#'      MC. Provided to allow the user to compute own summary statistics.}
#'   \item{\code{samplePsi}}{Array of sampled values for psi. \code{nSamples} x
#'      [number of origin sites] x [number of target sites]. Provided to allow
#'      the user to compute own summary statistics.}
#'   \item{\code{pointPsi}}{Simple point estimate of psi matrix.}
#'   \item{\code{pointMC}}{Simple point estimate of MC, using the point
#'      estimates of \code{psi} and \code{originRelAbund}.}
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
#'      these should be equivalent.  Estimated as the
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
#'      using \code{originPoints} and \code{targetPoints}.  NULL when
#'      \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#'   \item{\code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr}}{Summary
#'      statistics for continuous correlation bootstraps.  NULL when
#'      \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#'   \item{\code{inputSampleSize}}{If \code{sampleSize} was provided, this is
#'      that.  If not, it is either the calculated sample size (if that can be
#'      done), or left at NULL.}
#' }
#' @example inst/examples/estMCExamples.R
#' @export
estMC <- function(originDist, targetDist, originRelAbund, psi = NULL,
                  sampleSize = NULL,
                  originSites = NULL, targetSites = NULL,
                  originPoints = NULL, targetPoints = NULL,
                  originAssignment = NULL, targetAssignment = NULL,
                  originNames = NULL, targetNames = NULL,
                  nSamples = 1000, nSim = 1000, isGL = FALSE,
                  geoBias = NULL, geoVCov = NULL, row0 = 0,
                  verbose = 0,  calcCorr = FALSE, alpha = 0.05,
                  approxSigTest = FALSE, sigConst = 0,
                  projection.dist.calc = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs") {
  if (is.null(psi)) {
    return(estMCGlGps(isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
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
                      projection.dist.calc = projection.dist.calc))
  }
  else {
    return(estMCCmrAbund(originRelAbund = originRelAbund,
                         sampleSize = sampleSize, psi = psi,
                         originDist = originDist, targetDist = targetDist,
                         originSites=originSites, targetSites=targetSites,
                         nSamples = nSamples, row0 = row0, verbose=verbose,
                         alpha = alpha, approxSigTest = approxSigTest,
                         sigConst = sigConst))
  }
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
  download.file(url1, temp, mode = 'wb')
  fm <- readRDS(temp)
  unlink(temp)
  return(fm)
}
