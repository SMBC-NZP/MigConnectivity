###############################################################################
# Simulates position of birds by individual, season, year, and month.
# Incorporates migratory connectivity, movement within season, and dispersal
# between seasons.
# Does not incorporate births or deaths.
###############################################################################
#' Simulates position of birds by individual, season, year, and month.
#'
#' Incorporates migratory connectivity, movement within season, and dispersal
#' between seasons. Does not incorporate births or deaths.
#'
#' @param breedingAbund Vector with number of birds to simulate starting at
#'    each breeding site.
#' @param breedingDist Distances between the breeding sites. Symmetric matrix.
#' @param winteringDist Distances between the wintering sites. Symmetric
#'    matrix.
#' @param psi Transition probabilities between B origin and W target sites.
#'    A matrix with B rows and W columns where rows sum to 1.
#' @param nYears Number of years to simulate movement.
#' @param nMonths Number of months per breeding and wintering season.
#' @param winMoveRate Within winter movement rate.  Defaults to 0 (no
#'    movement).
#' @param sumMoveRate Within summer movement rate.  Defaults to 0 (no
#'    movement).
#' @param winDispRate Between winter dispersal rate.  Defaults to 0 (no
#'    dispersal).
#' @param sumDispRate Between summer dispersal rate.  Defaults to 0 (no
#'    dispersal).  Setting this to a value above 0 is equivalent to setting
#'    both natal and breeding dispersal to that same value.
#' @param natalDispRate Natal dispersal rate. Controls the movement of
#'    animals from their birthplace on their first return to the breeding
#'    grounds.  Defaults to 0 (return to the birthplace for all).
#' @param breedDispRate Breeding dispersal rate. Controls the movement of
#'    animals between breeding sites on spring migrations after the first.
#'    Defaults to 0 (return to the same breeding site each year).
#' @param verbose If set to a value > 0, informs the user on the passage
#'    of years and seasons during the simulation. Defaults to 0 (no output
#'    during simulation).
#'
#' @return \code{simMove} returns a list with elements:
#' \describe{
#'   \item{\code{animalLoc}}{\code{sum(breedingAbund)} (number of animals)
#'      by 2 by \code{nYears} by \code{nMonths} array with the simulated
#'      locations of each animal in each month of each season (summer or
#'      winter) of each year.  Values of cells are 1...B (first column) and
#'      1...W (second column) where B is the number of breeding sites and W is
#'      the number of wintering sites.}
#'   \item{\code{breedDispMat}}{B by B matrix of probabilities of breeding
#'      dispersal between each pair of 1...B breeding sites.  Direction is from
#'      row to column, so each row sums to 1.}
#'   \item{\code{natalDispMat}}{B by B matrix of probabilities of natal
#'      dispersal between each pair of 1...B breeding sites.  Direction is from
#'      row to column, so each row sums to 1.}
#'   \item{\code{sumMoveMat}}{B by B matrix of probabilities of within season
#'      movement between each pair of 1...B breeding sites.  Direction is from
#'      row to column, so each row sums to 1.}
#'   \item{\code{winDispMat}}{W by W matrix of probabilities of dispersal
#'      between each pair of 1...W nonbreeding sites.  Direction is from
#'      row to column, so each row sums to 1.}
#'   \item{\code{winMoveMat}}{W by W matrix of probabilities of within season
#'      movement between each pair of 1...W nonbreeding sites.  Direction is
#'      from row to column, so each row sums to 1.}
#' }
#'
#' @export
# @examples
simMove <- function(breedingAbund, breedingDist, winteringDist, psi,
                    nYears = 10, nMonths = 3, winMoveRate = 0,
                    sumMoveRate = 0, winDispRate = 0, sumDispRate = 0,
                    natalDispRate = 0, breedDispRate = 0, verbose = 0)
{
  nSeasons <- 2
  nBreeding <- length(breedingAbund)
  nWintering <- nrow(winteringDist)

  if (sumDispRate>0 && (natalDispRate>0 || breedDispRate>0) &&
      (sumDispRate!=natalDispRate || sumDispRate!=breedDispRate))
    stop("Can't specify summer dispersal separately from breeding or natal dispersal")
  if (sumDispRate<0 || natalDispRate<0 || breedDispRate<0 ||sumMoveRate<0 ||
      winMoveRate<0)
    stop("Can't specify negative movement or dispersal rates")

  # Turn rate terms into probabilities
  if (winMoveRate>0) {
    winMoveMat <- mlogit.mat(1/sqrt(winMoveRate), winteringDist)
  }
  else
    winMoveMat <- NULL
  if (sumMoveRate>0) {
    sumMoveMat <- mlogit.mat(1/sqrt(sumMoveRate), breedingDist)
  }
  else
    sumMoveMat <- NULL
  if (winDispRate>0) {
    winDispMat <- mlogit.mat(1/sqrt(winDispRate), winteringDist)
  }
  else
    winDispMat <- NULL
  if (sumDispRate>0) {
    natalDispRate <- sumDispRate
    breedDispRate <- sumDispRate
  }
  if (natalDispRate>0) {
    natalDispMat <- mlogit.mat(1/sqrt(natalDispRate), breedingDist)
  }
  else
    natalDispMat <- NULL
  if (breedDispRate>0) {
    breedDispMat <- mlogit.mat(1/sqrt(breedDispRate), breedingDist)
  }
  else
    breedDispMat <- NULL

  # Storage of locations
  animalLoc <- array(NA, c(sum(breedingAbund), nSeasons, nYears, nMonths))
  animalLoc[,1,1,1] <- rep(1:nBreeding, breedingAbund)

  # Run simulation
  for (y in 1:nYears) {
    if (verbose>0)
      cat("Year", y, "Summer, ")
    if (nMonths > 1)
      for (sm in 2:nMonths) {
        if (sumMoveRate==0)
          animalLoc[,1,y,sm] <- animalLoc[,1,y,sm-1]
        else
          for (i in 1:sum(breedingAbund))
            animalLoc[i,1,y,sm] <- which(rmultinom(1,1, sumMoveMat[animalLoc[i,1,y,sm-1], ])>0)
      }
    if (verbose>0)
      cat("Fall, ")
    if (y == 1) {
      for (i in 1:sum(breedingAbund))
        animalLoc[i,2,y,1] <- which(rmultinom(1,1, psi[animalLoc[i,1,y,1], ])>0)
    }
    else if (winDispRate==0)
      animalLoc[,2,y,1] <- animalLoc[,2,y-1,1]
    else
      for (i in 1:sum(breedingAbund))
        animalLoc[i,2,y,1] <- which(rmultinom(1,1, winDispMat[animalLoc[i,2,y-1,1], ])>0)
    if (verbose>0)
      cat("Winter, ")
    if (nMonths > 1)
      for (wm in 2:nMonths) {
        if (winMoveRate==0)
          animalLoc[,2,y,wm] <- animalLoc[,2,y,wm-1]
        else
          for (i in 1:sum(breedingAbund))
            animalLoc[i,2,y,wm] <- which(rmultinom(1,1, winMoveMat[animalLoc[i,2,y,wm-1], ])>0)
      }
    if (verbose>0)
      cat("Spring\n")
    if (y == 1 & nYears>1) {
      if (natalDispRate==0)
        animalLoc[,1,y+1,1] <- animalLoc[,1,y,1]
      else
        for (i in 1:sum(breedingAbund))
          animalLoc[i,1,y+1,1] <- which(rmultinom(1,1, natalDispMat[animalLoc[i,1,y,1], ])>0)
    }
    else if (y < nYears) {
      if (breedDispRate==0)
        animalLoc[,1,y+1,1] <- animalLoc[,1,y,1]
      else
        for (i in 1:sum(breedingAbund))
          animalLoc[i,1,y+1,1] <- which(rmultinom(1,1, breedDispMat[animalLoc[i,1,y,1], ])>0)
    }
  }
  return(list(animalLoc = animalLoc, natalDispMat = natalDispMat,
              breedDispMat = breedDispMat, sumMoveMat = sumMoveMat,
              winDispMat = winDispMat, winMoveMat = winMoveMat))
}


### Function for generating simulated count data ----

simCountData <- function (nPops, routePerPop, nYears, alphaPop, beta,
                          sdPop, sdRoute, sdYear){

  nRoutes <- nPops*routePerPop # Total number of routes
  # Generate data structure to hold counts and log (lambda)
  C <- log.expectedCount <- array(NA, dim = c(nYears, nRoutes))

  # Generate covariate values
  year <- 1:nYears
  yr <- (year - (nYears/2))/(nYears/2) # Standardize
  pops <- rep(1:nPops, each=routePerPop)

  # Draw two sets of random effects from their respective distributions
  epsRoute <- rnorm(n = nRoutes, mean = 0, sd = sdRoute)
  epsYear <- rnorm(n = nYears, mean = 0, sd = sdYear)

  # Loop over routes
  for (i in 1:nRoutes){

    # Build up systematic part of the GLM including random effects
    log.expectedCount[,i] <- alphaPop[pops[i]] + beta*yr + epsRoute[i] +
      epsYear
    expectedCount <- exp(log.expectedCount[,i])

    C[,i] <- rpois(n = nYears, lambda = expectedCount)
  }

  return(list(nPops = nPops, nRoutes = nRoutes, nYears = nYears,
              routePerPop = routePerPop,
              alphaPop = alphaPop, pop = pops, epsRoute = epsRoute,
              epsYear = epsYear, beta = beta, year = yr, sdPop = sdPop,
              sdRoute = sdRoute, sdYear = sdYear,
              expectedCount = expectedCount, C = C))
}


modelCountDataJAGS <- function (sim_data, ni, nt, nb, nc) {
  nPops <- sim_data$nPops
  nRoutes <- sim_data$nRoutes
  # Initial values
  jags.inits <- function()list(mu = runif(1,0,2), alpha = runif(nPops, -1,1), beta1 = runif(1,-1,1),
                               tau.alpha = runif(1,0,0.1), tau.noise = runif(1,0,0.1),
                               tau.rte = runif(1,0,0.1), route = runif(nRoutes,-1,1))
  # Parameters to monitor
  params <- c("mu", "alpha", "beta1", "sd.alpha", "sdRoute", "sd.noise", "totalN", "popN", "relN")
#   jags.data <- list(C = sim_data$C, nPops=sim_data$nPops, nRoutes = sim_data$nRoutes,
#                     routePerPop=sim_data$nRoutes/sim_data$nPops,
#                     year = sim_data$year, nYears = length(sim_data$year), pop = sim_data$pop)

  out <- jags(sim_data, inits=jags.inits, params, "sim_Poisson2.txt",
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          working.directory = file.path(R_PACKAGE_DIR, 'inst', paste0('JAGS', R_ARCH)))
  return(as.mcmc(out))
}




