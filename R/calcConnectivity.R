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
                        check.attributes = F))) {
    if (ncol(targetDist)!=nOrigin || ncol(originDist)!=nTarget ||
        !isTRUE(all.equal(colSums(psi), rep(1, nTarget), tolerance = 1e-6,
                          check.attributes = F)))
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
  MC <- sum(M *
              (matrix(rep(originDist, nTarget^2), nOrigin^2, nTarget^2) - mu.bD) *
              (matrix(rep(targetDist, nOrigin^2), nOrigin^2, nTarget^2, byrow=T) - mu.wD))/(sd.bD*sd.wD)
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
              (matrix(rep(targetDist, nOrigin^2), nOrigin^2, nTarget^2, byrow=T) - mu.wD))/(sd.bD*sd.wD)
  return(MC)
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
calcPsiMC <- function(originDist, targetDist, originRelAbund, locations, verbose=F) {
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

