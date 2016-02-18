###############################################################################
# MC = Migratory connectivity strength function
###############################################################################
#' Migratory connectivity strength function
#'
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'    matrix.
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'    matrix.
#' @param psi Transition probabilities between B origin and W target sites.
#'    Matrix with B rows and W columns where rows sum to 1.
#' @param originRelAbund Relative abundances at B origin sites. Numeric vector
#'    of length B that sums to 1.
#'
#' @return real value between -1 and 1, inclusive (usually between 0 and 1)
#'    indicating the degree of migratory connectivity.
#'
#' @example inst/examples/calcMCExamples.R
calcMC <- function(originDist, targetDist, psi, originRelAbund) {
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
  # Relative abundance at each combination of origin and target sites
  sumWinRelN <- apply(psi, 2, "*", originRelAbund)
  # Average origin distance between any two given birds.  Using matrix product (%*%) to get all probabilities of two birds origin in each site.
  mu.bD <- sum(originDist * (originRelAbund %*% t(originRelAbund)))
  # SD in origin distance between any two given birds
  sd.bD <- sqrt(sum((originDist - mu.bD)^2 * (originRelAbund %*% t(originRelAbund))))
  # Average target distance between any two given birds.  Using Kronecker product (%x%) to get all probabilities of two birds origin and target in each pair of sites.
  mu.wD <- sum(targetDist * matrix(colSums(sumWinRelN %x% sumWinRelN), nTarget, nTarget))
  # SD in target distance between any two given birds
  sd.wD <- sqrt(sum((targetDist - mu.wD)^2 * matrix(colSums(sumWinRelN %x% sumWinRelN), nTarget, nTarget)))
  # Correlation between origin and target distances, calculated from formula:
  # sum(probability of each combination of origin sites and target sites *
  # (distance between origin sites - mean origin site distance) * (distance between target sites - mean target site distance)) /
  # (SD in origin site distances * SD in target site distances)
  MC <- sum((sumWinRelN %x% sumWinRelN) *
              (matrix(rep(originDist, nTarget^2), nOrigin^2, nTarget^2) - mu.bD) *
              (matrix(rep(targetDist, nOrigin^2), nOrigin^2, nTarget^2, byrow=T) - mu.wD))/(sd.bD*sd.wD)
  return(MC)
}

###############################################################################
# rho function for individuals
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
  MC <- calcMC(originDist, targetDist, psiMat, originRelAbund)
  return(list(psi=psiMat, MC=MC))
}

