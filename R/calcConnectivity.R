###############################################################################
# MC = Migratory connectivity strength function
###############################################################################
#' Migratory connectivity strength function
#'
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'    matrix.
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'    matrix.
#' @param originRelAbund Relative abundances at B origin sites. Numeric vector
#'    of length B that sums to 1.
#' @param psi Transition probabilities between B origin and W target sites.
#'    Matrix with B rows and W columns where rows sum to 1.
#' @param sampleSize Total sample size of animals that psi was calculated from.
#'    Should be the number of animals released in one of the origin sites and
#'    observed in one of the target sites.  Optional, but recommended.
#'
#' @return scalar real value, usually between 0 and 1 (can be negative),
#' indicating the strength of migratory connectivity.
#'
#' If \code{sampleSize} is provided, this function uses the standard (relative
#' abundance and small-sample size corrected) formula for MC.  If not, it uses
#' the MC(R) formula, which only corrects for relative abundance.
#'
#' @example inst/examples/calcMCExamples.R
#' @export
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513 - 524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
calcMC <- function(originDist, targetDist, originRelAbund, psi, sampleSize=NULL) {
  nOrigin <- nrow(psi)
  nTarget <- ncol(psi)
  # Check inputs
  if (nrow(originDist)!=ncol(originDist) || !isSymmetric(unname(originDist)))
    stop('originDist (distances between origin sites) must be a square, symmetric matrix.')
  if (nrow(targetDist)!=ncol(targetDist) || !isSymmetric(unname(targetDist)))
    stop('targetDist (distances between target sites) must be a square, symmetric matrix.')
  if (ncol(originDist)!=nOrigin || ncol(targetDist)!=nTarget ||
      !isTRUE(all.equal(rowSums(psi), rep(1, nOrigin), tolerance = 1e-6,
                        check.attributes = FALSE))) {
    if (ncol(targetDist)!=nOrigin || ncol(originDist)!=nTarget ||
        !isTRUE(all.equal(colSums(psi), rep(1, nTarget), tolerance = 1e-6,
                          check.attributes = FALSE)))
      stop('psi should have [number of origin sites] rows and [number of target sites] columns and rows that sum to 1.')
    else { #Just turn it around
      psi <- t(psi)
      nOrigin <- nrow(psi)
      nTarget <- ncol(psi)
    }
  }
  if (length(originRelAbund)!=nOrigin || !isTRUE(all.equal(sum(originRelAbund), 1,
                                                       tolerance = 1e-6)))
    stop('originRelAbund must be a vector with [number of origin sites] values that sum to 1.')
  if (!is.null(sampleSize)) {
    if (length(sampleSize)==1)
      originAbund <- originRelAbund * sampleSize
    else
      stop('sampleSize must be a single positive number.')
    return(calcMCSmall(originDist = originDist, targetDist = targetDist,
                       psi = psi, originAbund = originAbund))
  }
  # Relative abundance at each combination of origin and target sites
  sumWinRelN <- apply(psi, 2, "*", originRelAbund)
  # Average origin distance between any two given birds.  Using matrix product (%*%) to get all probabilities of two birds origin in each site.
  relAbundProd <- originRelAbund %*% t(originRelAbund)
  mu.bD <- sum(originDist * relAbundProd)
  # SD in origin distance between any two given birds
  sd.bD <- sqrt(sum((originDist - mu.bD)^2 * relAbundProd))
  # Average target distance between any two given birds.  Using Kronecker product (%x%) to get all probabilities of two birds origin and target in each pair of sites.
  M <- sumWinRelN %x% sumWinRelN
  winRelNProd <- matrix(colSums(M), nTarget, nTarget)
  mu.wD <- sum(targetDist * winRelNProd)
  # SD in target distance between any two given birds
  sd.wD <- sqrt(sum((targetDist - mu.wD)^2 * winRelNProd))
  # Correlation between origin and target distances, calculated from formula:
  # sum(probability of each combination of origin sites and target sites *
  # (distance between origin sites - mean origin site distance) * (distance between target sites - mean target site distance)) /
  # (SD in origin site distances * SD in target site distances)
  MC <- c((c(originDist - mu.bD) / sd.bD) %*% M %*% (c(targetDist - mu.wD) / sd.wD))
  return(MC)
}

###############################################################################
# Calculate MC for known actual abundance per site
###############################################################################
calcMCSmall <- function(originDist, targetDist, originAbund, psi) {
  nOrigin <- nrow(psi)
  nTarget <- ncol(psi)
  # Check inputs
 if (length(originAbund)!=nOrigin)
   stop('originAbund must be a vector with [number of origin sites] values >= 1.')
  originRelAbund <- originAbund / sum(originAbund)
  originRelAbund2 <- (matrix(rep(originAbund, each = nOrigin), nOrigin, nOrigin) -
                        diag(nrow = nOrigin)) *
    matrix(rep(originRelAbund, nOrigin), nOrigin, nOrigin) / (sum(originAbund) - 1)
  # Relative abundance at each combination of origin and target sites
  sumWinRelN <- apply(psi, 2, "*", originRelAbund)
  sumWinN <- apply(psi, 2, "*", originAbund)
  subMat <- matrix(0, nOrigin ^ 2, nTarget ^ 2)
  subMat[rep(1:nOrigin, nOrigin) == rep(1:nOrigin, each = nOrigin),
         rep(1:nTarget, nTarget) == rep(1:nTarget, each = nTarget)] <- 1
  M <- (sumWinRelN %x% matrix(1, nOrigin, nTarget)) *
    ((matrix(1, nOrigin, nTarget) %x% sumWinN) - subMat) / (sum(originAbund) - 1) #rep(sumWinRelN, each = nOrigin), nOrigin^2, nTarget^2)
  sumWinRelN2 <- matrix(colSums(M), nTarget, nTarget)
  # Average origin distance between any two given birds.  Using matrix product (%*%) to get all probabilities of two birds origin in each site.
  mu.bD <- sum(originDist * originRelAbund2)
  # SD in origin distance between any two given birds
  sd.bD <- sqrt(sum((originDist - mu.bD)^2 * originRelAbund2))
  # Average target distance between any two given birds.  Using Kronecker product (%x%) to get all probabilities of two birds origin and target in each pair of sites.
  mu.wD <- sum(targetDist * sumWinRelN2)
  # SD in target distance between any two given birds
  sd.wD <- sqrt(sum((targetDist - mu.wD)^2 * sumWinRelN2))
  # Correlation between origin and target distances, calculated from formula:
  # sum(probability of each combination of origin sites and target sites *
  # (distance between origin sites - mean origin site distance) * (distance between target sites - mean target site distance)) /
  # (SD in origin site distances * SD in target site distances)
  MC <- sum(M *
              (matrix(rep(originDist, nTarget^2), nOrigin^2, nTarget^2) - mu.bD) *
              (matrix(rep(targetDist, nOrigin^2), nOrigin^2, nTarget^2,
                      byrow = TRUE) - mu.wD))/(sd.bD*sd.wD)
  return(MC)
}

#' Calculate Mantel correlation (rM) from points and/or distances.
#'
#' Calculation of rM from POINTS geolocators and/or GPS
#' data, not accounting for uncertainty. If you've already calculated
#' distances between points, you can use those instead.
#'
#' @param targetPoints A sf \code{POINTS} object, with length number of animals
#'     tracked.  Each point indicates the point estimate location in the
#'     non-release season.
#' @param originPoints A sf \code{POINTS} object, with length number of animals
#'     tracked.  Each point indicates the release location of an animal.
#' @param targetDist Distances between the target locations of the tracked
#'    animals.  Symmetric matrix with number of animals rows and columns,
#'    although really you only need the lower triangle filled in.
#' @param originDist Distances between the origin locations of the tracked
#'    animals.  Symmetric matrix with number of animals rows and columns,
#'    although really you only need the lower triangle filled in.
#'
#' @return \code{calcMantel} returns a list with elements:
#' \describe{
#'   \item{\code{pointCorr}}{Simple point estimate of Mantel correlation.}
#'   \item{\code{originDist, targetDist}}{Distances between each pair of
#'   \code{originPoints} and each pair of \code{targetPoints}, respectively,
#'   in meters. If you used distances as inputs instead, then these are just
#'   what you fed in.}
#' }
#' @export
#'
#' @examples
#' data(OVENdata) # Ovenbird
#' rM0 <- calcMantel(originPoints = OVENdata$originPoints, # Capture Locations
#'                   targetPoints = OVENdata$targetPoints) # Target locations
#' str(rM0)
#' @seealso \code{\link{estMantel}}, \code{\link{calcMC}}, \code{\link{estMC}}
#'
#' @references
#' Ambrosini, R., A. P. Moller, and N. Saino. 2009. A quantitative measure of
#' migratory connectivity. Journal of Theoretical Biology 257:203-211.
#' \href{https://doi.org/10.1016/j.jtbi.2008.11.019}{doi:10.1016/j.jtbi.2008.11.019}

calcMantel <- function(targetPoints = NULL, originPoints = NULL,
                       targetDist = NULL, originDist = NULL) {

  if (is.null(targetPoints) && is.null(targetDist)){
    stop('Define either targetPoints or targetDist')
  }
  if (is.null(originPoints) && is.null(originDist)){
    stop('Define either originPoints or originDist')
  }
  if (!is.null(targetDist))
    nAnimals <- dim(targetDist)[1]
  else {

    if(!is.null(targetPoints) & !inherits(targetPoints, "sf")){
     targetPoints <- sf::st_as_sf(targetPoints)}
    nAnimals <- nrow(targetPoints)

    if(is.na(sf::st_crs(targetPoints))){
      stop('Coordinate system definition needed for targetPoints')
    }
    # NEED A CHECK HERE TO ENSURE THAT targetPoints and originPoints
    # are sf objects


    targetDist <- matrix(NA, nAnimals, nAnimals)

    targetDist[lower.tri(targetDist)] <- 1

    distIndices <- which(!is.na(targetDist), arr.ind = TRUE)

    # project target points to WGS #
    targetPoints2 <- sf::st_transform(targetPoints,4326)
    targetDist0 <- geosphere::distGeo(sf::st_coordinates(targetPoints2[distIndices[,'row'],]),
                                      sf::st_coordinates(targetPoints2[distIndices[,'col'],]))

    targetDist[lower.tri(targetDist)] <- targetDist0
    diag(targetDist) <- 0
    targetDist <- t(targetDist)
    targetDist[lower.tri(targetDist)] <- targetDist0
  }
  if (is.null(originDist)) {
    if(!is.null(originPoints) & !(inherits(originPoints, "sf"))){
      originPoints <- sf::st_as_sf(originPoints)
    }
    if(is.na(sf::st_crs(originPoints))) {
      stop('Coordinate system definition needed for originPoints')
    }
    originPoints2 <- sf::st_transform(originPoints, 4326)
    originDist <- matrix(NA, nAnimals, nAnimals)

    originDist[lower.tri(originDist)] <- 1

    distIndices <- which(!is.na(originDist), arr.ind = TRUE)
    originDist0 <- geosphere::distGeo(sf::st_coordinates(originPoints2[distIndices[,'row'],]),
                                      sf::st_coordinates(originPoints2[distIndices[,'col'],]))
    originDist[lower.tri(originDist)] <- originDist0
    diag(originDist) <- 0
    originDist <- t(originDist)
    originDist[lower.tri(originDist)] <- originDist0
  }
  pointCorr <- ncf::mantel.test(originDist, targetDist, resamp=0, quiet = TRUE)$correlation
  return(list(pointCorr = pointCorr, originDist = originDist, targetDist = targetDist))
}


###############################################################################
# Mantel function for individuals
###############################################################################
calcStrengthInd <- function(originDist, targetDist, locations, resamp=1000, verbose = 0) {
  nInd <- dim(locations)[1]
  originDist2 <- targetDist2 <- matrix(0, nInd, nInd)
  for (i in 1:(nInd-1)) {
    for (j in (i+1):nInd) {
      originDist2[i,j] <- originDist2[j,i] <- originDist[locations[i,1,1,1], locations[j,1,1,1]]
      targetDist2[i,j] <- targetDist2[j,i] <- targetDist[locations[i,2,1,1], locations[j,2,1,1]]
    }
  }
  return(ncf::mantel.test(originDist2, targetDist2, resamp=resamp, quiet = !verbose))
}

###############################################################################
# Simple approach to estimate psi matrix and MC from simulated (or real) data
# (doesn't include uncertainty)
###############################################################################
calcPsiMC <- function(originDist, targetDist, originRelAbund, locations,
                      verbose=FALSE) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  psiMat <- matrix(0, nOrigin, nTarget)
  nInd <- dim(locations)[1]
  nYears <- dim(locations)[3]
  nMonths <- dim(locations)[4]
  for (i in 1:nInd) {
    if (i %% 1000 == 0 && verbose) #
      cat("Individual", i, "of", nInd, "\n")
    originMat <- locations[i,1,,]
    targetMat <- locations[i,2,,]
    bIndices <- which(!is.na(originMat))
    wIndices <- which(!is.na(targetMat))
    if (length(bIndices) && length(wIndices))
      for (bi in bIndices)
        for (wi in wIndices)
          psiMat[originMat[bi], targetMat[wi]] <- psiMat[originMat[bi], targetMat[wi]] + 1
  }
  psiMat <- apply(psiMat, 2, "/", rowSums(psiMat))
  MC <- calcMC(originDist, targetDist, psi = psiMat,
               originRelAbund = originRelAbund, sampleSize = nInd)
  return(list(psi=psiMat, MC=MC))
}

#' Reverse transition probabilities and origin relative abundance
#'
#' Reverse transition probabilities (psi; sum to 1 for each origin site) and
#' origin relative abundance (originRelAbund; sum to 1 overall) estimates to
#' calculate or estimate target site to origin site transition probabilities
#' (gamma; sum to 1 for each target site), target site relative abundances
#' (targetRelAbund; sum to 1 overall), and origin/target site combination
#' probabilities (pi; sum to 1 overall). If either psi or originRelAbund is an
#' estimate with sampling uncertainty expressed, this function will propagate
#' that uncertainty to provide true estimates of gamma, targetRelAbund, and pi;
#' otherwise (if both are simple point estimates), it will also provide point
#' estimates.
#'
#' Alternatively, can be used to reverse migratory combination (joint)
#' probabilities (pi; sum to 1 overall) to psi, originRelAbund, gamma, and
#' targetRelAbund.
#'
#' @param psi Transition probabilities between B origin and W target sites.
#'  Either a matrix with B rows and W columns where rows sum to 1, an array with
#'  dimensions x, B, and W (with x samples of the transition probability matrix
#'  from another model), an 'estPsi' object (result of calling estTransition),
#'  or a MARK object with estimates of transition probabilities
#' @param originRelAbund Relative abundance estimates at B origin sites. Either
#'  a numeric vector of length B that sums to 1 or an mcmc object with at least
#'  \code{nSamples} rows and columns including 'relN[1]' through 'relN[B]'
#' @param pi Migratory combination (joint) probabilities. Either a matrix with B
#'  rows and W columns where all entries sum to 1, an array with dimensions x,
#'  B, and W, or an 'estPi' object (currently only the results of calling this
#'  function) Either pi or psi and originRelAbund should be specified.
#' @param originSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are origin
#' @param targetSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are target
#' @param originNames Vector of names for the origin sites. If not provided, the
#'  function will try to get them from psi
#' @param targetNames Vector of names for the target sites. If not provided, the
#'  function will try to get them from psi
#' @param nSamples Number of times to resample \code{psi} and/or
#'  \code{originRelAbund}. The purpose is to estimate sampling uncertainty;
#'  higher values here will do so with more precision
#' @param row0 If \code{originRelAbund} is an mcmc object or array, this can be
#'  set to 0 (default) or any greater integer to specify where to stop ignoring
#'  samples (additional "burn-in")
#' @param alpha Level for confidence/credible intervals provided. Default (0.05)
#'  gives 95 percent CI
#'
#' @return If both psi and originRelAbund are simple point estimates,
#' \code{reversePsiRelAbund} returns a list with point estimates of gamma,
#' targetRelAbund, and pi. Otherwise, it returns a list with the elements:
#' \describe{
#'   \item{\code{gamma}}{List containing estimates of reverse transition
#'   probabilities:
#'   \itemize{
#'    \item{\code{sample}} Array of sampled values for gamma. \code{nSamples} x
#'      [number of target sites] x [number of origin sites]. Provided to allow
#'      the user to compute own summary statistics.
#'    \item{\code{mean}} Main estimate of gamma matrix. [number of target sites]
#'      x [number of origin sites].
#'    \item{\code{se}} Standard error of gamma, estimated from SD of
#'      \code{gamma$sample}.
#'    \item{\code{simpleCI}} \code{1 - alpha} confidence interval for gamma,
#'      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'      \code{gamma$sample}.
#'    \item{\code{bcCI}} Bias-corrected \code{1 - alpha} confidence interval
#'      for gamma. May be preferable to \code{simpleCI} when \code{mean} is the
#'      best estimate of gamma. \code{simpleCI} is preferred when
#'      \code{median} is a better estimator. When the mean and median are equal,
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of \code{sample},
#'      where z0 is the proportion of \code{sample < mean}.
#'    \item{\code{median}} Median estimate of gamma matrix.
#'    \item{\code{point}} Simple point estimate of gamma matrix, not accounting
#'      for sampling error.
#'   }
#'   }
#'   \item{\code{targetRelAbund}}{List containing estimates of relative
#'    abundance at target sites. Items within are the same as within gamma,
#'    except for having one fewer dimension.}
#'   \item{\code{pi}}{List containing estimates of origin/target site
#'    combination probabilities (sum to 1). Items within are the same as within
#'    gamma, except for reversing dimensions (same order as psi).}
#'   \item{\code{input}}{List containing the inputs to \code{reversePsiRelAbund}.}
#' }
#' If the input is pi instead of psi and originRelAbund, then pi is not an
#' output, but psi and originRelAbund are. Otherwise the same.
#'
#' @export
#'
#' @example inst/examples/reversePsiRelAbundExamples.R
#'
reverseTransition <- function(psi = NULL, originRelAbund = NULL, pi = NULL,
                               originSites=NULL, targetSites=NULL,
                               originNames = NULL, targetNames = NULL,
                               nSamples = 1000, row0 = 0, alpha = 0.05) {
  if ((is.null(psi) || is.null(originRelAbund)) && is.null(pi))
    stop("Either psi and originRelAbund or pi must be specified")
  if (is.null(psi)){
    if (is.matrix(pi)) {
      if (is.null(originNames)) {
        originNames <- dimnames(pi)[[1]]
      }
      if (is.null(targetNames)) {
        targetNames <- dimnames(pi)[[2]]
      }
      gamma <- t(prop.table(pi, 2))
      targetRelAbund <- colSums(pi)
      psi <- prop.table(pi, 1)
      originRelAbund <- rowSums(pi)
      return(list(psi = psi, originRelAbund = originRelAbund, gamma = gamma,
                  targetRelAbund = targetRelAbund))
    }
    else {
      return(reverseEstTransition(psi = psi, originRelAbund = originRelAbund,
                                  pi = pi,
                                  originSites=originSites,
                                  targetSites=targetSites,
                                  originNames = originNames,
                                  targetNames = targetNames,
                                  nSamples = nSamples, row0 = row0,
                                  alpha = alpha))
    }
  }
  if (is.matrix(psi) && is.numeric(originRelAbund)) {
    nOriginSites <- nrow(psi)
    if (is.null(originNames)) {
      originNames <- dimnames(psi)[[1]]
    }
    if (is.null(targetNames)) {
      targetNames <- dimnames(psi)[[2]]
    }
    if (length(originRelAbund) != nOriginSites ||
        !isTRUE(all.equal(sum(originRelAbund), 1, tolerance = 1e-6)))
      stop('originRelAbund must be a vector with [number of origin sites/number of rows in psi] values that sum to 1.')
    pi <- sweep(psi, 1, originRelAbund, "*")
    dimnames(pi) <- list(originNames, targetNames)
    gamma <- t(prop.table(pi, 2))
    targetRelAbund <- colSums(pi)
    return(list(gamma = gamma, targetRelAbund = targetRelAbund, pi = pi))
  }
  else
    return(reverseEstTransition(psi = psi, originRelAbund = originRelAbund,
                                 originSites=originSites,
                                 targetSites=targetSites,
                                 originNames = originNames,
                                 targetNames = targetNames,
                                 nSamples = nSamples, row0 = row0,
                                 alpha = alpha))
}

# This version based on multi-logit link (others are for testing only, I hope)

divCoefNLL <- function(psi_r, banded, reencountered, counts) {
  nOriginSites <- nrow(reencountered)
  nTargetSites <- ncol(reencountered)
  #print(psi_r)
  psi <- matrix(psi_r[1:(nOriginSites * (nTargetSites - 1))], nOriginSites,
                nTargetSites - 1)
  for (o in 1:nOriginSites)
    psi[o,] <- exp(psi[o,]) / (1 + sum(exp(psi[o, ])))
  psi <- cbind(psi, 1 - rowSums(psi))
  #psi[psi<0] <- 0
  #print(psi)
  r <- stats::plogis(psi_r[(nOriginSites * (nTargetSites - 1)) + 1:nTargetSites])
  #print(r)
  p <- sweep(psi, 2, r, "*")
  p <- cbind(p, 1 - rowSums(p))
  reencountered <- cbind(reencountered, banded - rowSums(reencountered))
  d1 <- d <- rep(0, length = nOriginSites)
  if (any(is.nan(psi)) || any(is.nan(r)) || any(psi<0) || any(r<0) || any(r>1))
    return(0)
  for (o in 1:nOriginSites) {
    d[o] <- stats::dmultinom(x = reencountered[o, ], size = banded[o], prob = p[o, ],
                      log = TRUE)
    if (!is.null(counts))
      d1[o] <- stats::dmultinom(x = counts[o, ], prob = psi[o, ], log = TRUE)
  }
  #print(d)
  #print(d1)
  if (any(r==0) || any(r==1))
    d <- d/10
  return(-sum(d) - sum(d1))
}

divCoefGrad <- function(psi_r, banded, reencountered, counts) {
  nOriginSites <- nrow(reencountered)
  nTargetSites <- ncol(reencountered)
  if (is.null(counts))
    ss <- reencountered
  else
    ss <- reencountered + counts
  lost <- banded - rowSums(reencountered)
  phi <- matrix(psi_r[1:(nOriginSites * (nTargetSites - 1))], nOriginSites,
                nTargetSites - 1)
  psi <- matrix(0, nOriginSites, nTargetSites - 1)
  for (o in 1:nOriginSites)
    psi[o,] <- exp(phi[o,]) / (1 + sum(exp(phi[o, ])))
  psi <- cbind(psi, 1 - rowSums(psi))
  rho <- psi_r[(nOriginSites * (nTargetSites - 1)) + 1:nTargetSites]
  r <- 1 / (1 + exp(-rho))
  dphi <- matrix(0, nOriginSites, nTargetSites - 1)
  for (o in 1:nOriginSites) {
    for (t in 1:(nTargetSites - 1)) {
      dphi[o, t] <- sum(ss[o, ] * psi[o, t]) - ss[o, t] -
        lost[o] * psi[o, t] * (r[nTargetSites] - r[t] * (1 + sum(exp(phi[o, ]))) +
                                 sum(r[-nTargetSites] * exp(phi[o, ]))) /
        (1 + sum((1 - r[-nTargetSites]) * exp(phi[o, ])) - r[nTargetSites])
    }
  }
  drho <- rep(0, nTargetSites)
  for (t in 1:nTargetSites){
    drho[t] <- sum(-reencountered[, t] / (1 + exp(rho[t])) +
                     lost * psi[, t] / (2 + exp(rho[t]) + exp(-rho[t])) /
                     (1 - apply(sweep(psi, 2, r, "*"), 1, sum)))
  }
  return(c(dphi, drho))
}

#deriv(~n*log(psi) + m*log(psi*(exp(rho) / (1 + exp(rho)))), "rho")

#' Calculate psi (transition probabilities between sites in two phases of the
#' annual cycle)
#'
#' Provides simple maximum-likelihood point estimate of transition probabilities
#' that does not include measures of uncertainty. Incorporates detection
#' heterogeneity where appropriate (band/ring return data), but not location
#' uncertainty. Shared primarily for testing; use of \code{\link{estTransition}}
#' is recommended instead.
#'
#' @param banded For band return data, a vector of the number of released
#'  animals from each origin site (including those never reencountered
#'  in a target region)
#' @param reencountered For band return data, a matrix with B rows and W
#'  columns. Number of animals reencountered
#'  on each target site by origin site they came from
#' @param counts Migration data without target-region detection heterogeneity
#'  (i.e., anything but band return data) can be entered one of two ways: either
#'  here or with \code{originAssignment} and \code{targetAssignment}. If here, a
#'  matrix with B rows and W columns with counts of animals observed in each
#'  combination of origin and target site
#' @param originAssignment Assignment of animals (not including band return
#'  data) to origin season sites. A vector of integers (1-B) with length number
#'  of animals tracked. Note that these data can either be entered using this
#'  argument and \code{targetAssignment} or using \code{counts}, but not both
#' @param targetAssignment Assignment of animals (not including band return
#'  data) to target season sites. A vector of integers (1-W) with length number
#'  of animals tracked. Note that these data can either be entered using this
#'  argument and \code{originAssignment} or using \code{counts}, but not both
#' @param originNames Optional, but recommended to keep track. Vector of names
#'  for the origin sites. If not provided, the function will either try to get
#'  these from another input or provide default names (capital letters)
#' @param targetNames Optional, but recommended to keep track. Vector of names
#'  for the target sites. If not provided, the function will either try to get
#'  these from another input or provide default names (numbers)
#' @param method See \code{optim}. "SANN" is slow but reasonably accurate.
#'  "Nelder-Mead" is fast but not accurate here. "BFGS" is also fast but not
#'  very stable here
#'
#' @return \code{calcTransition} returns a list with the element(s):
#' \describe{
#'   \item{\code{psi}}{Matrix with point estimate of transition probabilities}
#'   \item{\code{r}}{Vector containing point estimate of reencounter
#'    probabilities at each target site. Not included unless data includes
#'    band reencounters}
#' }
#'
#' @export
#'
#' @examples
#' nOriginSites <- 4
#' nTargetSites <- 4
#' originNames <- LETTERS[1:nOriginSites]
#' targetNames <- 1:nTargetSites
#' psiTrue <- array(c(0.1, 0.2, 0.3, 0.4,
#'                    0.2, 0.3, 0.4, 0.1,
#'                    0.3, 0.4, 0.1, 0.2,
#'                    0.4, 0.1, 0.2, 0.3),
#'                  c(nOriginSites, nTargetSites),
#'                  dimnames = list(originNames, targetNames))
#' rowSums(psiTrue)
#' rTrue <- c(0.5, 0.05, 0.3, 0.6)
#' banded1 <- c(1000, 2000, 3000, 5000)
#' reencountered1 <- simCMRData(psiTrue, banded1, rTrue)$reencountered
#' psi_r_calc <- calcTransition(banded = banded1,
#'                              reencountered = reencountered1,
#'                              originNames = originNames,
#'                              targetNames = targetNames,
#'                              method = "SANN")
#' psi_r_calc
#' \dontrun{
#' psi_r_mcmc <- estTransition(banded = banded1, reencountered = reencountered1,
#'                             originNames = originNames,
#'                             targetNames = targetNames,
#'                             method = "MCMC",
#'                             nThin = 1, nSamples = 60000, nBurnin = 20000,
#'                             verbose = 1)
#' psi_r_mcmc
#' psi_r_boot <- estTransition(banded = banded1, reencountered = reencountered1,
#'                             originNames = originNames,
#'                             targetNames = targetNames,
#'                             method = "bootstrap",
#'                             nSamples = 1000, verbose = 2)
#' psi_r_boot
#' }
#' @seealso \code{\link{estTransition}}, \code{\link{optim}}
calcTransition <- function(banded = NULL, reencountered = NULL, counts = NULL,
                           originAssignment = NULL, targetAssignment = NULL,
                           originNames = NULL, targetNames = NULL,
                           method = "SANN") {
  if (is.null(counts) && length(originAssignment)>0) {
    if (is.null(originNames))
      nOriginSites <- length(unique(originAssignment))
    else
      nOriginSites <- length(originNames)
    if (is.null(targetNames))
      nTargetSites <- length(unique(targetAssignment))
    else
      nTargetSites <- length(targetNames)
    counts <- as(table(factor(originAssignment, levels = 1:nOriginSites),
                       factor(targetAssignment, levels = 1:nTargetSites)),
                 "matrix")
    dimnames(counts) <- list(originNames, targetNames)
  }
  if (is.null(banded))
    return(list(psi = prop.table(counts, 1)))
  nOriginSites <- nrow(reencountered)
  nTargetSites <- ncol(reencountered)
  bunded <- banded + 0.00001 # In case of zeros
  reuncountered <- reencountered + 0.00001
  startPsiR <- sweep(reuncountered, 1, bunded, "/")# + 0.00001
  # startPsi <- prop.table(startPsiR, 1)
  # if (!is.null(counts))
  #   startPsi <- (startPsi * sum(reencountered) +
  #                  prop.table(counts + 0.00001, 1) * sum(counts)) /
  #   (sum(reencountered) + sum(counts))
  startR <- colSums(startPsiR)
  startR[startR>0.9999] <- 0.9999
  if (is.null(counts))
    startPsi <- sweep(reuncountered, 2, startR, "/")
  else
    startPsi <- sweep(reuncountered, 2, startR, "/") + counts
  startPsi <- prop.table(startPsi, 1)
  # print(startPsi); print(startR)
  # print(banded); print(reencountered)
  startPar <- c(log(sweep(startPsi[, -nTargetSites, drop = FALSE], 1,
                          startPsi[, nTargetSites], "/")),
                stats::qlogis(startR))
  #print(startPar)
  # print(divCoefNLL(startPar, banded = banded, reencountered = reencountered,
  #                  counts = counts))
  opt1 <- optim(fn = divCoefNLL,
                par = startPar,
                gr = divCoefGrad,
                method = method,
                # control = list(trace = 0, REPORT = 1,
                #                ndeps = c(rep(0.001, nOriginSites * (nTargetSites - 1)),
                #                          rep(0.0001, nTargetSites)),
                #                type = 1,
                #                maxit = 5000),
                # lower = rep(0, (nOriginSites + 1) * nTargetSites),
                # upper = rep(1, (nOriginSites + 1) * nTargetSites),
                banded = banded, reencountered = reencountered, counts = counts)
  psi_r <- opt1$par
  psi <- matrix(psi_r[1:(nOriginSites * (nTargetSites - 1))], nOriginSites,
                nTargetSites - 1)
  for (o in 1:nOriginSites)
    psi[o,] <- exp(psi[o,]) / (1 + sum(exp(psi[o, ])))
  psi <- cbind(psi, 1 - rowSums(psi))
  r <- stats::plogis(psi_r[(nOriginSites * (nTargetSites - 1)) + 1:nTargetSites])
  dimnames(psi) <- list(originNames, targetNames)
  names(r) <- targetNames
  return(list(psi = psi, r = r))
}

calcPi <- function(banded = NULL, reencountered = NULL, counts = NULL,
                   originAssignment = NULL, targetAssignment = NULL,
                   originNames = NULL, targetNames = NULL, method = "BFGS",
                   relAbundOriginData = NULL, relAbundTargetData = NULL) {

}


#' @rdname reverseTransition
#' @export
reversePsiRelAbund <- reverseTransition

#' @rdname reverseTransition
#' @export
reverseTransitionRelAbund <- reverseTransition

#' @rdname reverseTransition
#' @export
reversePi <- reverseTransition

#' @rdname calcTransition
#' @export
calcPsi <- calcTransition

#' @rdname calcMC
#' @export
calcStrength <- calcMC


