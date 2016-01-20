###############################################################################
# MC = Migratory connectivity strength function
###############################################################################
#' Migratory connectivity strength function
#'
#' @param originDist
#' @param targetDist
#' @param psi
#' @param originRelN
#'
#' @return real value between -1 and 1, inclusive (usually between 0 and 1)
#'    indicating the degree of migratory connectivity.
#'
# @examples
calcMC <- function(originDist, targetDist, psi, originRelN) {
  nOrigin <- nrow(psi)
  nTarget <- ncol(psi)
  # Check inputs
  if (nrow(originDist)!=ncol(originDist) || !isSymmetric(unname(originDist)))
    stop('originDist (distances between origin sites) must be a square, symmetric matrix.')
  if (nrow(targetDist)!=ncol(targetDist) || !isSymmetric(unname(targetDist)))
    stop('targetDist (distances between target sites) must be a square, symmetric matrix.')
  if (ncol(originDist)!=nOrigin || ncol(targetDist)!=nTarget || any(rowSums(psi)!=1)) {
    if (ncol(targetDist)!=nOrigin || ncol(originDist)!=nTarget || any(colSums(psi)!=1))
      stop('psi should have [number of origin sites] rows and [number of target sites] columns and rows that sum to 1.')
    else { #Just turn it around
      psi <- t(psi)
      nOrigin <- nrow(psi)
      nTarget <- ncol(psi)
    }
  }
  if (length(originRelN)!=nOrigin || sum(originRelN)!=1)
    stop('originRelN must be a vector with [number of origin sites] values that sum to 1.')
  # Relative abundance at each combination of origin and target sites
  sumWinRelN <- apply(psi, 2, "*", originRelN)
  # Average origin distance between any two given birds.  Using matrix product (%*%) to get all probabilities of two birds origin in each site.
  mu.bD <- sum(originDist * (originRelN %*% t(originRelN)))
  # SD in origin distance between any two given birds
  sd.bD <- sqrt(sum((originDist - mu.bD)^2 * (originRelN %*% t(originRelN))))
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
calcStrengthInd <- function(originDist, targetDist, locations, resamp=1000, latlon=F) {
  nInd <- dim(locations)[1]
  originDist2 <- targetDist2 <- matrix(0, nInd, nInd)
  for (i in 1:(nInd-1)) {
    for (j in (i+1):nInd) {
      originDist2[i,j] <- originDist2[j,i] <- originDist[locations[i,1,1,1], locations[j,1,1,1]]
      targetDist2[i,j] <- targetDist2[j,i] <- targetDist[locations[i,2,1,1], locations[j,2,1,1]]
    }
  }
  return(ncf::mantel.test(originDist2, targetDist2, resamp=resamp, latlon=latlon))
}

###############################################################################
# Simple approach to estimate psi matrix and MC from simulated (or real) data
# (doesn't include uncertainty)
###############################################################################
calcPsiMC <- function(originDist, targetDist, originRelN, locations, verbose=F) {
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
  MC <- calcMC(originDist, targetDist, psiMat, originRelN)
  return(list(psi=psiMat, MC=MC))
}

