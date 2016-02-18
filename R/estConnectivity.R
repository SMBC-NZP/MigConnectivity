###############################################################################
#' Estimate MC from abundance and transition probability estimates.
#'
#' Resampling of uncertainty for MC from RMark psi matrix estimates and/or JAGS
#' relative abundance MCMC samples.
#'
#' @param originRelAbund Relative abundances at B origin sites. Either a
#'  numeric vector of length B that sums to 1 or an mcmc object with
#'  \code{nSamples} rows  and columns including 'relN[1]' through 'relN[B]'.
#' @param psi Transition probabilities between B origin and W target sites.
#'  Either a matrix with B rows and W columns where rows sum to 1 or a MARK
#'  object with estimates of transition probabilities.
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'  matrix.
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'  matrix.
#' @param originSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are origin.
#' @param targetSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are target.
#' @param nSamples Number of times to resample \code{psi} and/or
#'  \code{originRelAbund}.
#' @param row0 If \code{originRelAbund} is an mcmc object, this can be set
#'  to 0 (default) or any greater integer to specify where to stop ignoring
#'  samples ("burn-in").
#' @param verbose 0 (default), 1, or 2. 0 prints no output during run. 1 prints
#'  a line every 100 samples.  2 prints a line every sample.
#' @param alpha Level for confidence intervals provided.
#' @return \code{estMCCmrAbund} returns a list with elements:
#' \describe{
#'   \item{\code{sampleMC}}{\code{nSamples} sampled values for MC. Provided
#'      to allow the user to compute own summary statistics.}
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
#'      for MC.  Preferable to \code{simpleCI} when \code{pointMC} is the
#'      best estimate of MC. \code{simpleCI} is preferred when
#'      \code{medianMC} is a better estimator. When \code{pointMC==medianMC},
#'      these should be equivalent.}
#'   \item{\code{hpdCI}}{\code{1 - alpha} credible interval for MC,
#'      estimated using the highest posterior density (HPD) method.}
#' }
#' @example inst/examples/estMCCmrAbundExamples.R
estMCCmrAbund <- function(originRelAbund, psi, originDist, targetDist,
                             originSites=NULL, targetSites=NULL,
                             nSamples = 1000, row0 = 0, verbose=0,
                             alpha = 0.05) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
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
  pointMC <- calcMC(originDist, targetDist, psiBase, abundBase)
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
    sampleMC[i] <- calcMC(originDist, targetDist, psiNew, abundNew)
  }
  meanMC <- mean(sampleMC)
  medianMC <- median(sampleMC)
  seMC <- sd(sampleMC)
  # Calculate confidence intervals using quantiles of sampled MC
  simpleCI <- quantile(sampleMC, c(alpha/2, 1-alpha/2), na.rm=T, type = 8)
  z0 <- qnorm(sum((sampleMC)<pointMC)/nSamples)
  bcCI <- quantile(sampleMC, pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=T, type = 8)
  MC.mcmc <- coda::as.mcmc(sampleMC) # Ha!
  hpdCI <- coda::HPDinterval(MC.mcmc, 1-alpha)
  return(list(sampleMC=sampleMC, samplePsi = psi.array, pointPsi = psiBase,
              pointMC=pointMC, meanMC=meanMC,
              medianMC=medianMC, seMC=seMC, simpleCI=simpleCI,
              bcCI=bcCI, hpdCI=hpdCI, sampleCorr = NULL, pointCorr = NULL,
              meanCorr = NULL, medianCorr = NULL, seCorr=NULL,
              simpleCICorr=NULL, bcCICorr=NULL))
}

###############################################################################
#
###############################################################################
#' Resampling of uncertainty from geolocators and/or GPS data
#'
#' @param isGL Vector indicating which animals were tracked with geolocators.
#'    Should be length of number of animals tracked, with TRUE for animals in
#'    \code{targetPoints} with geolocators and FALSE for animals with GPS.
#' @param geoBias Vector of length 2 indicating expected bias in longitude and
#'    latitude of \code{targetPoints}, in meters.
#' @param geoVCov 2x2 matrix with expected variance/covariance in longitude and
#'    latitude of \code{targetPoints}, in meters.
#' @param originRelAbund Relative abundances at B origin sites. A numeric
#'    vector of length B that sums to 1.
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'    matrix.
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'    matrix.
#' @param targetPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the point estimate location in
#'    the non-release season.
#' @param targetSites Geographic definition of sites in the non-release season.
#' @param targetAssignment Optional. Point estimate assignment of
#'    \code{targetPoints} to non-release season sites. Integer vector with
#'    length number of animals tracked.
#' @param originPoints A \code{SpatialPoints} object, with length number of
#'    animals tracked.  Each point indicates the release location of an animal.
#' @param originSites Geographic definition of sites in the release season.
#' @param originAssignment Assignment of \code{originPoints} to release season
#'    sites. Integer vector with length number of animals tracked. Optional,
#'    but either \code{originAssignment} or \code{originSites} and
#'    \code{originPoints} should be defined.
#' @param originNames Optional. Vector of names for the release season sites.
#' @param targetNames Optional. Vector of names for the non-release season
#'    sites.
#' @param nBoot Number of bootstrap runs. Animals are sampled with replacement
#'    for each of these to estimate sampling uncertainty.
#' @param verbose Integer 0-3 for level of output during bootstrap: 0 = none,
#'    1 = every 10, 2 = every run, 3 = every animal.
#' @param nSim Number of times to sample random points for each animal from
#'    parametric distribution of non-release season error. Ignored for GPS
#'    points (assumed to have no geographic error).
#' @param calcCorr In addition to MC, should function also estimate continuous
#'    correlation between release and non-release locations?  Default is TRUE.
#' @param alpha Level for confidence intervals provided.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{sampleMC}}{\code{nBoot} sampled values for MC. Provided
#'      to allow the user to compute own summary statistics.}
#'   \item{\code{samplePsi}}{Array of sampled values for psi. \code{nBoot} x
#'      [number of origin sites] x [number of target sites]. Provided
#'      to allow the user to compute own summary statistics.}
#'   \item{\code{pointSites}}{Matrix of point assignment of number of animals
#'      to each origin and target site combination.}
#'   \item{\code{pointPsi}}{Simple point estimate of psi matrix, as
#'      \code{prop.table(pointSites, 1)}.}
#'   \item{\code{pointMC}}{Simple point estimate of MC, using \code{pointPsi}
#'      and \code{originRelAbund}.}
#'   \item{\code{meanMC, medianMC}}{Mean and median of \code{sampleMC}.
#'      Estimates of MC incorporating parametric uncertainty.}
#'   \item{\code{seMC}}{Standard error of MC, estimated from SD of
#'      \code{sampleMC}.}
#'   \item{\code{simpleCI}}{\code{1 - alpha} confidence interval for MC,
#'      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'      \code{sampleMC}.}
#'   \item{\code{bcCI}}{Bias-corrected \code{1 - alpha} confidence interval
#'      for MC.  Preferable to \code{simpleCI} when \code{pointMC} is the
#'      best estimate of MC. \code{simpleCI} is preferred when
#'      \code{medianMC} is a better estimator. When \code{pointMC==medianMC},
#'      these should be equivalent.}
#'   \item{\code{sampleCorr}}{\code{nBoot} sampled values for continuous
#'      correlation. Provided to allow the user to compute own summary
#'      statistics.}
#'   \item{\code{pointCorr}}{Simple point estimate of continuous correlation,
#'      using \code{originPoints} and \code{targetPoints}.}
#'   \item{\code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr}}{Summary
#'      statistics for continuous correlation bootstraps.}
#' }
#'
# @examples
estMCGlGps <- function(isGL, geoBias, geoVCov, originRelAbund,
                          originDist, targetDist, targetPoints,
                          targetSites, targetAssignment=NULL,
                          originPoints=NULL, originSites=NULL,
                          originAssignment=NULL, originNames=NULL,
                          targetNames=NULL, nBoot = 1000, verbose=0,
                          nSim = 1000, calcCorr=T, alpha = 0.05) {
  # Input checking and assignment
  if (!(verbose %in% 0:3))
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = every animal")
  if (length(geoBias)!=2)
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in meters)")
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = F)))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in meters)")
  if ((is.null(originPoints) || is.null(originSites)) &&
      is.null(originAssignment))
    stop("Need to define either originAssignment or originSites and originPoints")
  if (calcCorr && is.null(originPoints))
    stop('If calcCorr is TRUE, need to define originPoints')
  nAnimals <- length(targetPoints)
  if (length(isGL)==1)
    isGL <- rep(isGL, nAnimals)
  if (is.null(originAssignment))
    originAssignment <- sp::over(originPoints, originSites)
  if (is.null(targetAssignment))
    targetAssignment <- sp::over(targetPoints, targetSites)
  nOriginSites <- length(unique(originAssignment))
  nTargetSites <- length(targetSites)
  if (length(targetPoints)!=nAnimals || length(targetAssignment)!=nAnimals ||
      length(originAssignment)!=nAnimals)
    stop("isGL should be the same length as originAssignment/originPoints and targetPoints/targetAssignment (number of animals)")
  if (any(is.na(originAssignment)))
    stop("NAs in origin sites (make sure all points fall within polygons)")
  if (length(originRelAbund)!=nOriginSites || sum(originRelAbund)!=1)
    stop('originRelAbund should be vector with length number of origin sites that sums to 1')
  if (dim(originDist)!=rep(nOriginSites,2) ||
      dim(targetDist)!=rep(nTargetSites,2))
    stop('Distance matrices should be square with same number of sites of each type as assignments/points (with distances in meters)')
  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))

 WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
 Lambert<-"+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

  MC <- corr <- rep(NA, nBoot)
  geoBias2 <- matrix(rep(geoBias, nSim), nrow=nSim, byrow=T)
  nAnimals <- length(originAssignment)

  # Point estimate of MC
  pointSites <- array(0, c(nOriginSites, nTargetSites),
                       dimnames = list(originNames, targetNames))
  for(i in 1:nAnimals)
    pointSites[originAssignment[i], targetAssignment[i]] <-
      pointSites[originAssignment[i], targetAssignment[i]] + 1
  pointPsi <- prop.table(pointSites, 1)
  pointMC <- calcMC(originDist,targetDist,pointPsi,originRelAbund)

  if (calcCorr) {
    targetDist1 <- matrix(NA, nAnimals, nAnimals)
    targetDist1[lower.tri(targetDist1)] <- 1
    distIndices <- which(!is.na(targetDist1), arr.ind = T)
    targetPoints2 <- sp::spTransform(targetPoints, sp::CRS(WGS84))
    targetDist0 <- geosphere::distVincentyEllipsoid(targetPoints2[distIndices[,'row'],],
                                         targetPoints2[distIndices[,'col'],])
    targetDist1[lower.tri(targetDist1)] <- targetDist0
    originPoints2 <- sp::spTransform(originPoints, sp::CRS(WGS84))
    originDistStart <- matrix(geosphere::distVincentyEllipsoid(originPoints2[
      rep(1:nAnimals, nAnimals)], originPoints2[rep(1:nAnimals,
                                                     each=nAnimals)]),
      nAnimals, nAnimals)
    originDist1 <- originDistStart[1:nAnimals, 1:nAnimals]
    pointCorr <- ade4::mantel.rtest(as.dist(originDist1), as.dist(targetDist1),
                               nrepet=0)
  }
  for (boot in 1:nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have birds from every origin site
    origin.sample <- 'Filler' # Start with one origin site
    while (length(unique(origin.sample))<nOriginSites) {
      # Sample individual animals with replacement
      animal.sample <- sample(1:nAnimals, nAnimals, replace=T)
      # Get origin points for those animals
      origin.point.sample <- originPoints[animal.sample]
      # Get origin population for each animal sampled
      origin.sample <- originAssignment[animal.sample]
    }
    target.sample <- rep(NA, nAnimals)
    target.point.sample <- matrix(NA, nAnimals, 2)
    for(i in 1:nAnimals){
      if (verbose > 2)
        cat('\tAnimal', i, 'of', nAnimals)
      draws <- 0
      if (isGL[animal.sample[i]]) {
        while (is.na(target.sample[i])) {
          draws <- draws + 1
          # Sample random point for each bird from parametric distribution of NB error
          point.sample <- sp::SpatialPoints(t(rmvnorm(nsim=nSim, mu=cbind(
            targetPoints@coords[animal.sample[i],1],
            targetPoints@coords[animal.sample[i],2]), V=geoVCov))+
              geoBias2, sp::CRS(Lambert))
          # filtered to stay in NB areas (land)
          target.sample0 <- sp::over(point.sample, targetSites)
          target.sample[i]<-target.sample0[!is.na(target.sample0)][1]
        }
        target.point.sample[i, ]<-point.sample[!is.na(target.sample0)][1]@coords
      }
      else { # Assume no location error for GPS
        target.sample[i] <- targetAssignment[animal.sample[i]]
        target.point.sample[i, ] <- targetPoints[animal.sample[i],]@coords
      }
      if (verbose > 2)
        cat(' ', draws, 'draws\n')
      # Now that we have breeding and non-breeding site for point, add to transition count matrix
      sites.array[boot, origin.sample[i], target.sample[i]] <-
          sites.array[boot, origin.sample[i], target.sample[i]] + 1
    }
    # Create psi matrix as proportion of those from each breeding site that went to each NB site
    psi.array[boot, , ] <- prop.table(sites.array[boot, , ], 1)
    # Calculate MC from that psi matrix
    MC[boot] <- calcMC(originDist, targetDist, psi.array[boot, , ],
                       originRelAbund)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=T), "SD:", sd(MC, na.rm=T),
          "low quantile:", quantile(MC, alpha/2, na.rm=T),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=T), "\n")
    if (calcCorr) {
      originDist1 <- originDistStart[animal.sample, animal.sample]
      target.point.sample <- sp::SpatialPoints(target.point.sample,CRS(Lambert))
      target.point.sample2 <- sp::spTransform(target.point.sample,CRS(WGS84))
      targetDist0 <- geosphere::distVincentyEllipsoid(target.point.sample2[distIndices[,'row']],
                                           target.point.sample2[distIndices[,'col']])
      targetDist1[lower.tri(targetDist1)] <- targetDist0
      corr[boot] <- ade4::mantel.rtest(as.dist(originDist1),as.dist(targetDist1),
                                 nrepet=0)
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=T), "SD:", sd(corr, na.rm=T),
            "low quantile:", quantile(corr, alpha/2, na.rm=T),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=T), "\n")
    }
  }
  MC.z0 <- qnorm(sum((MC)<pointMC)/nBoot)
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=T, type = 8)
  MC.mcmc <- coda::as.mcmc(MC) # Ha!
  hpdCI <- coda::HPDinterval(MC.mcmc, 1-alpha)
  if (calcCorr) {
    meanCorr <- mean(corr)
    medianCorr <- median(corr)
    seCorr <- sd(corr)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=T, type = 8)
    corr.z0 <- qnorm(sum((corr)<pointCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                           na.rm=T, type = 8)
  } else
    pointCorr <- meanCorr <- medianCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  return(list(sampleMC = MC, samplePsi = psi.array,
              pointPsi = pointPsi, pointMC = pointMC, meanMC = mean(MC),
              medianMC = median(MC), seMC = sd(MC),
              simpleCI = quantile(MC, c(alpha/2, 1-alpha/2), na.rm=T, type = 8),
              bcCI = bcCI, hpdCI = hpdCI, sampleCorr = corr, pointCorr = pointCorr,
              meanCorr = meanCorr, medianCorr = medianCorr, seCorr=seCorr,
              simpleCICorr=simpleCICorr, bcCICorr=bcCICorr))
}

#' @describeIn estMCGlGps Convenience function if all points are from geolocators.
estMCGl <- function(geoBias, geoVCov, originRelAbund, originDist,
                      targetDist, targetPoints, targetSites,
                      targetAssignment=NULL, originPoints=NULL,
                      originSites=NULL, originAssignment=NULL,
                      originNames=NULL, targetNames=NULL, nBoot = 1000,
                      verbose=0, nSim = 1000, calcCorr=T, alpha = 0.05) {
  if (is.null(originPoints) && is.null(originAssignment))
    stop("Need to define either originAssignment or originSites and originPoints")
  else if (is.null(originAssignment))
    nAnimals <- length(originPoints)
  else
    nAnimals <- length(originAssignment)
  isGL <- rep(T, nAnimals)
  return(estMCGlGps(isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                       originRelAbund=originRelAbund,
                       originDist=originDist,
                       targetDist=targetDist,
                       targetPoints=targetPoints, targetSites=targetSites,
                       targetAssignment=targetAssignment,
                       originPoints=originPoints, originSites=originSites,
                       originAssignment=originAssignment,
                       originNames=originNames, targetNames=targetNames,
                       nBoot = nBoot, verbose=verbose, nSim = nSim,
                       calcCorr=calcCorr, alpha = alpha))
}

#' @describeIn estMCGlGps Convenience function if all points are from GPS.
estMCGps <- function(originPoints, targetPoints, originSites,
                       targetSites, originRelAbund, originDist,
                       targetDist, originNames=NULL, targetNames=NULL,
                       nBoot = 1000, verbose=0, calcCorr=T, alpha = 0.05) {
  if (is.null(originPoints) && is.null(originAssignment))
    stop("Need to define either originAssignment or originSites and originPoints")
  else if (is.null(originAssignment))
    nAnimals <- length(originPoints)
  else
    nAnimals <- length(originAssignment)
  isGL <- rep(F, nAnimals)
  geoBias <- rep(0,2); geoVCov <- matrix(0,2,2)
  return(estMCGlGps(isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                       originRelAbund=originRelAbund,
                       originDist=originDist,
                       targetDist=targetDist,
                       targetPoints=targetPoints, targetSites=targetSites,
                       targetAssignment=targetAssignment,
                       originPoints=originPoints, originSites=originSites,
                       originAssignment=originAssignment,
                       originNames=originNames, targetNames=targetNames,
                       nBoot = nBoot, verbose=verbose,
                       nSim = 0, calcCorr=calcCorr, alpha = alpha))
}

###############################################################################
#' Estimate MC from abundance and/or transition probability estimates OR
#' geolocator and/or GPS data.
#'
#' Resampling of uncertainty for MC from RMark psi matrix estimates and/or JAGS
#' relative abundance MCMC samples.
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
#' @return \code{estMC} returns a list with elements:
#' \describe{
#'   \item{\code{sampleMC}}{\code{nSamples} or \code{nBoot} sampled values for
#'      MC. Provided to allow the user to compute own summary statistics.}
#'   \item{\code{samplePsi}}{Array of sampled values for psi. \code{nBoot}
#'      OR \code{nSamples} x [number of origin sites] x [number of target
#'      sites]. Provided to allow the user to compute own summary statistics.}
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
#'      for MC.  Preferable to \code{simpleCI} when \code{pointMC} is the
#'      best estimate of MC. \code{simpleCI} is preferred when
#'      \code{medianMC} is a better estimator. When \code{pointMC==medianMC},
#'      these should be equivalent.}
#'   \item{\code{hpdCI}}{\code{1 - alpha} credible interval for MC,
#'      estimated using the highest posterior density (HPD) method.}
#'   \item{\code{sampleCorr}}{\code{nBoot} sampled values for continuous
#'      correlation. Provided to allow the user to compute own summary
#'      statistics.  NULL when \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#'   \item{\code{pointCorr}}{Simple point estimate of continuous correlation,
#'      using \code{originPoints} and \code{targetPoints}.  NULL when
#'      \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#'   \item{\code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr}}{Summary
#'      statistics for continuous correlation bootstraps.  NULL when
#'      \code{calcCorr==FALSE} or \code{!is.null(psi)}.}
#' }
#' @example inst/examples/estMCExamples.R
estMC <- function(originDist, targetDist, originRelAbund, psi = NULL,
                  originSites = NULL, targetSites = NULL,
                  originPoints = NULL, targetPoints = NULL,
                  originAssignment = NULL, targetAssignment = NULL,
                  originNames = NULL, targetNames = NULL,
                  nSamples = 1000, nSim = 1000, isGL = FALSE,
                  geoBias = NULL, geoVCov = NULL, row0 = 0,
                  verbose = 0,  calcCorr = FALSE, alpha = 0.05) {
  if (is.null(psi)) {
    return(estMCGlGps(isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                      originRelAbund=originRelAbund,
                      originDist=originDist,
                      targetDist=targetDist,
                      targetPoints=targetPoints, targetSites=targetSites,
                      targetAssignment=targetAssignment,
                      originPoints=originPoints, originSites=originSites,
                      originAssignment=originAssignment,
                      originNames=originNames, targetNames=targetNames,
                      nBoot = nSamples, verbose=verbose,
                      nSim = nSim, calcCorr=calcCorr, alpha = alpha))
  }
  else {
    return(estMCCmrAbund(originRelAbund = originRelAbund, psi = psi,
                         originDist = originDist, targetDist = targetDist,
                         originSites=originSites, targetSites=targetSites,
                         nSamples = nSamples, row0 = row0, verbose=verbose,
                         alpha = alpha))
  }
}
