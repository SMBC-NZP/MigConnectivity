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
estMCCmrAbund <- function(originRelAbund, psi, originDist, targetDist,
                             originSites=NULL, targetSites=NULL,
                             nSamples = 1000, row0 = 0, verbose=0,
                             alpha = 0.05) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  if ('numeric' %in% class(originRelAbund)) {
    abundFixed <- TRUE
    if (length(originRelAbund)!=nOrigin)
      stop('Number of origin sites must be constant between distance matrix and abundance')
    abundBase <- originRelAbund
  }
  else {
    abundFixed <- FALSE
    abundParams <- paste('relN[', 1:nOrigin, ']', sep='')
    abundBase <- colMeans(originRelAbund[row0 + 1:nSamples, abundParams])
  }
  if ('matrix' %in% class(psi)) {
    psiFixed <- TRUE
    if (nrow(psi)!=nOrigin || ncol(psi)!=nTarget)
      stop('Size of psi matrix must be consistant with distance matrices')
    psiBase <- psi
  }
  else {
    psiFixed <- FALSE
    if (!('numeric' %in% class(originSites)) ||
        !('numeric' %in% class(targetSites)))
      stop('Must specify which RMark Psi parameters represent origin and target sites')
    psiBase <- TransitionMatrix(get.real(psi, "Psi", se=T))[originSites, targetSites]
    if (any(diag(psi$results$beta.vcv) < 0))
      stop("Can't sample model, negative beta variances")
  }
  pointMC <- calcMC(originDist, targetDist, psiBase, abundBase)
  sampleMC <- vector('numeric', nSamples)
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
      psiNew <- make.psiMat.rand(psi, originSites, targetSites)
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
  MC.mcmc <- as.mcmc(sampleMC) # Ha!
  hpdCI <- HPDinterval(MC.mcmc, 1-alpha)
  return(list(sampleMC=sampleMC, pointMC=pointMC, meanMC=meanMC,
              medianMC=medianMC, seMC=seMC, simpleCI=simpleCI,
              bcCI=bcCI, hpdCI=hpdCI))
}

###############################################################################
# Resampling of uncertainty from geolocators and/or GPS data
###############################################################################
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
  if (dim(geoVCov)!=c(2,2))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in meters)")
  if ((is.null(originPoints) || is.null(originSites)) &&
      is.null(originAssignment))
    stop("Need to define either originAssignment or originSites and originPoints")
  if (calcCorr && is.null(originPoints))
    stop('If calcCorr is TRUE, need to define originPoints')
  nAnimals <- length(isGL)
  if (is.null(originAssignment))
    originAssignment <- over(originPoints, originSites)
  if (is.null(targetAssignment))
    targetAssignment <- over(targetPoints, targetSites)
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
    targetDist <- matrix(NA, nAnimals, nAnimals)
    targetDist[lower.tri(targetDist)] <- 1
    distIndices <- which(!is.na(targetDist), arr.ind = T)
    targetPoints2<-spTransform(targetPoints,CRS(WGS84))
    targetDist0 <- distVincentyEllipsoid(targetPoints2[distIndices[,'row'],],
                                         targetPoints2[distIndices[,'col'],])
    targetDist[lower.tri(targetDist)] <- targetDist0
    originPoints2<-spTransform(originPoints,CRS(WGS84))
    originDistStart <- matrix(distVincentyEllipsoid(originPoints2[
      rep(1:nAnimals, nAnimals)], originPoints2[rep(1:nAnimals,
                                                     each=nAnimals)]),
      nAnimals, nAnimals)
    originDist <- originDistStart[1:nAnimals, 1:nAnimals]
    pointCorr <- mantel.rtest(as.dist(originDist), as.dist(targetDist),
                               nrepet=0)
  }
  for (boot in 1:nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
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
          point.sample<-SpatialPoints(t(rmvnorm(nsim=nSim, mu=cbind(
                          targetPoints@coords[animal.sample[i],1],
                          targetPoints@coords[animal.sample[i],2]), V=geoVCov))+
                            geoBias2, CRS(Lambert))
          # filtered to stay in NB areas (land)
          target.sample0 <- over(point.sample, targetSites)
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
    MC[boot] <- calcMC(originDist,targetDist,psi.array[boot, , ],
                       originRelAbund)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=T), "SD:", sd(MC, na.rm=T),
          "low quantile:", quantile(MC, alpha/2, na.rm=T),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=T), "\n")
    if (calcCorr) {
      originDist <- originDistStart[animal.sample, animal.sample]
      target.point.sample <- SpatialPoints(target.point.sample,CRS(Lambert))
      target.point.sample2<-spTransform(target.point.sample,CRS(WGS84))
      targetDist0 <- distVincentyEllipsoid(target.point.sample2[distIndices[,'row']],
                                           target.point.sample2[distIndices[,'col']])
      targetDist[lower.tri(targetDist)] <- targetDist0
      corr[boot] <- mantel.rtest(as.dist(originDist),as.dist(targetDist),
                                 nrepet=0)
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=T), "SD:", sd(corr, na.rm=T),
            "low quantile:", quantile(corr, alpha/2, na.rm=T),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=T), "\n")
    }
  }
  dimnames(psi.array)[2:3] <- list(originNames, targetNames)
  MC.z0 <- qnorm(sum((MC)<pointMC)/nBoot)
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=T, type = 8)
  if (calcCorr) {
    meanCorr <- mean(corr)
    seCorr <- sd(corr)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=T, type = 8)
    corr.z0 <- qnorm(sum((corr)<pointCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                           na.rm=T, type = 8)
  } else
    pointCorr <- meanCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  return(list(sampleMC = MC, samplePsi = psi.array, pointSites=pointSites,
              pointPsi = pointPsi, pointMC = pointMC, meanMC = mean(MC),
              seMC = sd(MC), simpleCI = quantile(MC, c(alpha/2, 1-alpha/2),
                                                      na.rm=T, type = 8),
              bcCI = bcCI, sampleCorr = corr, pointCorr = pointCorr,
              meanCorr = meanCorr, seCorr=seCorr,
              simpleCICorr=simpleCICorr, bcCICorr=bcCICorr))
}

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

estMCGps <- function(originPoints, targetPoints, originSites,
                       targetSites, originRelAbund, originDist,
                       targetDist, originNames=NULL, targetNames=NULL,
                       nBoot = 1000, verbose=0, nSim = 1000, calcCorr=T) {
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
                       nSim = nSim, calcCorr=calcCorr, alpha = alpha))
}
