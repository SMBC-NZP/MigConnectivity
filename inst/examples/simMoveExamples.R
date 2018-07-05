# Dispersal simulation ----
## Functions

###############################################################################
# Utility functions for use in simulations
###############################################################################
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

# Crude optimizable function for developing MC pattern based on MC strength
mlogitMC <- function(slope, MC.in, origin.dist, target.dist, origin.abund, sample.size) {
  nBreeding <- nrow(origin.dist)
  nWintering <- nrow(target.dist)
  psi <- mlogitMat(slope, origin.dist)
  if (any(psi<0))
    return(5*slope^2)
  MC <- calcMC(origin.dist, target.dist, originRelAbund = origin.abund, psi,
               sampleSize = sample.size)
  return((MC.in - MC)^2)
}

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

## Simulation ----
nBreeding <- 100
nWintering <- 100
breedingPos <- matrix(c(rep(seq(-99, -81, 2), each=sqrt(nBreeding)),
                        rep(seq(49, 31, -2), sqrt(nBreeding))), nBreeding, 2)
winteringPos <- matrix(c(rep(seq(-79, -61, 2), each=sqrt(nWintering)),
                         rep(seq(9, -9, -2), sqrt(nWintering))), nWintering, 2)
head(breedingPos)
tail(breedingPos)
head(winteringPos)
tail(winteringPos)

breedDist <- distFromPos(breedingPos, 'ellipsoid')
nonbreedDist <- distFromPos(winteringPos, 'ellipsoid')

# Breeding Abundance
breedingN <- rep(50, nBreeding) # Lower for example speed
breedingRelN <- breedingN/sum(breedingN)

\dontrun{
  # Breeding Abundance
  breedingN <- rep(5000, nBreeding)
  breedingRelN <- breedingN/sum(breedingN)
}

# Set up psi matrix
slope <- 1.923904
\dontrun{
  o <- optimize(mlogitMC, MC.in = 0.25, origin.dist = breedDist,
                target.dist = nonbreedDist, origin.abund = breedingRelN,
                sample.size = sum(breedingN),
                interval = c(0, 10), tol = .Machine$double.eps^0.5)
  o
  slope <- o$minimum
}
psi <- mlogitMat(slope, breedDist)

# Baseline strength of migratory connectivity
MC <- calcMC(breedDist, nonbreedDist, breedingRelN, psi, sum(breedingN))
round(MC, 4)

# Other basic simulation parameters

## Dispersal simulations---
set.seed(1516)
nYears <- 4 # Lower for example speed
nMonths <- 2 # Each season, lower for example speed
Drates <- c(0.04, 0.16)    #rates of dispersal, fewer for example speed
\dontrun{
  nYears <- 15
  nMonths <- 4 # Each season
  Drates <- c(0.02, 0.04, 0.08, 0.16, 0.32, 0.64)    #rates of dispersal

birdLocDisp <- vector('list', length(Drates))
Disp.df  <- data.frame(Year=rep(1:nYears, length(Drates)),
                       Rate=rep(Drates, each = nYears), MC = NA)
for(i in 1:length(Drates)){
  cat('Dispersal Rate', Drates[i], '\n')
  birdLocDisp[[i]] <- simMove(breedingN, breedDist, nonbreedDist, psi, nYears, nMonths,
                              sumDispRate = Drates[i])
  for(j in 1:nYears){
    cat('\tYear', j, '\n')
    temp.results <- calcPsiMC(breedDist, nonbreedDist, breedingRelN,
                                 birdLocDisp[[i]]$animalLoc, years = j)
    Disp.df$MC[j + (i - 1) * nYears] <- temp.results$MC
  }
} # end i loop

Disp.df$Year <- Disp.df$Year - 1 #just run once!
data.frame(Disp.df, roundMC = round(Disp.df$MC, 2), nearZero = Disp.df$MC < 0.01)

  write.csv(Disp.df, "../Disp.df.csv")

# Convert dispersal rates to probabilities of dispersing at least certain distance
threshold <- 1000
probFarDisp <- matrix(NA, nBreeding, length(Drates), dimnames = list(NULL, Drates))
for (i in 1:length(Drates)) {
  for (k in 1:nBreeding) {
    probFarDisp[k, i] <- sum(birdLocDisp[[i]]$natalDispMat[k, which(breedDist[k, ]>= threshold)])
  }
}
summary(probFarDisp)
}


#plot results
\dontrun{
  require(ggplot2)
  require(ggthemes)

  line_set=c(0.5, 0.75, 1, 1.25, 1.75, 2.25)
  png('../dispersal.plot.png', width = 6, height = 6, res = 600, units = 'in')
  ggplot(Disp.df, aes(x=Year, y=MC, size=as.factor(Rate)))+
    geom_line()+
    scale_size_manual(values=line_set)+
    labs(size="Dispersal Rate")+
    scale_y_continuous("MC", limits=c(-.005, 0.26), breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25))+
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15))+
    theme_bw()+
    theme(axis.title=element_text(size=16, face ="bold"),axis.text=element_text(size=14),
          panel.grid.major = element_line(color="grey90"),panel.grid.major.x = element_blank(),
          legend.title=element_text(size=14), legend.text=element_text(size=14))
  dev.off()
}
