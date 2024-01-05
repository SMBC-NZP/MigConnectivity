### Dispersal simulation ----
## Utility functions for use in simulations

# Simple approach to estimate psi matrix and MC from simulated (or real) data
# (doesn't include uncertainty).  Only uses one year for computation
calcPsiMC <- function(originDist, targetDist, originRelAbund, locations,
                      years = 1, months = 1, verbose=FALSE) {
  nOrigin <- nrow(originDist)
  nTarget <- nrow(targetDist)
  psiMat <- matrix(0, nOrigin, nTarget)
  nInd <- dim(locations)[1]
  nYears <- dim(locations)[3]
  nMonths <- dim(locations)[4]
  for (i in 1:nInd) {
    if (i %% 1000 == 0 && verbose) #
      cat("Individual", i, "of", nInd, "\n")
    originMat <- locations[i, 1, years, months]
    targetMat <- locations[i, 2, years, months]
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

## Simulation
originNames <- c("A", "B", "C")
nBreeding <- length(originNames) # Number of sites reduced for example speed
targetNames <- as.character(1:4)
nWintering <- length(targetNames)

psi <- matrix(c(0.5, 0.25, 0.15, 0.1,
                0.15, 0.4, 0.25, 0.2,
                0.1, 0.15, 0.2, 0.55), nBreeding, nWintering,
              TRUE, list(originNames, targetNames))
psi
breedingPos <- matrix(c(seq(-99, -93, 3),
                        rep(40, nBreeding)), nBreeding, 2)
winteringPos <- matrix(c(seq(-88, -82, 2),
                         rep(0, nWintering)), nWintering, 2)
breedingPos
winteringPos

breedDist <- distFromPos(breedingPos, 'ellipsoid')
nonbreedDist <- distFromPos(winteringPos, 'ellipsoid')

# Breeding Abundance
breedingN <- rep(50, nBreeding) # Reduced from 5000 for example speed
breedingRelN <- breedingN/sum(breedingN)


# Baseline strength of migratory connectivity
\donttest{
  MC <- calcMC(breedDist, nonbreedDist, breedingRelN, psi, sum(breedingN))
  round(MC, 4)
}
# Other basic simulation parameters

## Dispersal simulations---
set.seed(1516)
nYears <- 4 # Reduced from 15 for example speed
nMonths <- 2 # Each season, reduced from 4 for example speed
Drates <- c(0.04, 0.16) # Rates of dispersal, fewer for example speed
\donttest{

  birdLocDisp <- vector('list', length(Drates))
  Disp.df  <- data.frame(Year=rep(1:nYears, length(Drates)),
                         Rate=rep(Drates, each = nYears), MC = NA)
  for(i in 1:length(Drates)){
    cat('Dispersal Rate', Drates[i], '\n')
    birdLocDisp[[i]] <- simMove(breedingN, breedDist, nonbreedDist, psi, nYears,
                                nMonths, sumDispRate = Drates[i])
    for(j in 1:nYears){
      cat('\tYear', j, '\n')
      temp.results <- calcPsiMC(breedDist, nonbreedDist, breedingRelN,
                                   birdLocDisp[[i]]$animalLoc, years = j)
      Disp.df$MC[j + (i - 1) * nYears] <- temp.results$MC
    }
  } # end i loop

  Disp.df$Year <- Disp.df$Year - 1 #just run once!
  data.frame(Disp.df, roundMC = round(Disp.df$MC, 2),
             nearZero = Disp.df$MC < 0.01)

  # Convert dispersal rates to probabilities of dispersing at least certain
  # distance
  threshold <- 1000
  probFarDisp <- matrix(NA, nBreeding, length(Drates),
                        dimnames = list(NULL, Drates))
  for (i in 1:length(Drates)) {
    for (k in 1:nBreeding) {
      probFarDisp[k, i] <- sum(
        birdLocDisp[[i]]$natalDispMat[k, which(breedDist[k, ]>= threshold)])
    }
  }
  summary(probFarDisp)

  #plot results
  with(subset(Disp.df, Rate == 0.04), plot(Year, MC, "l", col = "blue",
                                           ylim = c(0, NA)))
  lines(Disp.df$Year, Disp.df$MC, col = "darkblue")
  legend("topright", legend = Drates, col = c("blue", "darkblue"), lty = 1)

}
