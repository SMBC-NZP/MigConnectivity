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
#' @example inst/examples/simMoveExamples.R
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513-524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
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
    winMoveMat <- mlogitMat(1/sqrt(winMoveRate), winteringDist)
  }
  else
    winMoveMat <- NULL
  if (sumMoveRate>0) {
    sumMoveMat <- mlogitMat(1/sqrt(sumMoveRate), breedingDist)
  }
  else
    sumMoveMat <- NULL
  if (winDispRate>0) {
    winDispMat <- mlogitMat(1/sqrt(winDispRate), winteringDist)
  }
  else
    winDispMat <- NULL
  if (sumDispRate>0) {
    natalDispRate <- sumDispRate
    breedDispRate <- sumDispRate
  }
  if (natalDispRate>0) {
    natalDispMat <- mlogitMat(1/sqrt(natalDispRate), breedingDist)
  }
  else
    natalDispMat <- NULL
  if (breedDispRate>0) {
    breedDispMat <- mlogitMat(1/sqrt(breedDispRate), breedingDist)
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


###############################################################################
# Function for generating simulated count data
###############################################################################
#' Simulates Breeding Bird Survey-style count data
#'
#'
#' @param nPops Number of populations/regions
#' @param routePerPop Vector of length 1 or nPops containing the number of routes (i.e. counts) per population. If length(routePerPop) == 1, number of routes is identical for each population
#' @param nYears Number of years surveys were conducted
#' @param alphaPop Vector of length 1 or nPops containing the log expected number of individuals counted at each route for each population. If length(alphaPop) == 1, expected counts are identical for each population
#' @param beta Coefficient of linear year effect (default = 0)
#' @param sdRoute Standard deviation of random route-level variation
#' @param sdYear Standard deviation of random year-level variation
#'
#' @return \code{simCountData} returns a list containing:
#'  \describe{
#'    \item{\code{nPops}}{Number of populations/regions.}
#'    \item{\code{nRoutes}}{Total number of routes.}
#'    \item{\code{nYears}}{Number of years.}
#'    \item{\code{routePerPop}}{Number of routes per population.}
#'    \item{\code{year}}{Vector of length nYears with standardized year values.}
#'    \item{\code{pop}}{Vector of length nRoutes indicating the population/region in which each route is located.}
#'    \item{\code{alphaPop}}{log expected count for each populations.}
#'    \item{\code{epsRoute}}{realized deviation from alphaPop for each route.}
#'    \item{\code{epsYear}}{realized deviation from alphaPop for each year.}
#'    \item{\code{beta}}{linear year effect.}
#'    \item{\code{sdRoute}}{standard deviation of random route-level variation.}
#'    \item{\code{sdYear}}{standard deviation of random year-level variation.}
#'    \item{\code{expectedCount}}{nRoutes by nYears matrix containing deterministic expected counts.}
#'    \item{\code{C}}{nRoutes by nYears matrix containing observed counts.}
#'   }
#'
#'
#' @export
#' @example inst/examples/simCountExamples.R
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513-524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
#'
#' Link, W. A. and J. R. Sauer. 2002. A hierarchical analysis of population
#' change with application to Cerulean Warblers. Ecology 83: 2832-2840.

simCountData <- function (nPops, routePerPop, nYears, alphaPop, beta = 0, sdRoute, sdYear){

  if(length(routePerPop) == 1){
    nRoutes <- nPops*routePerPop # Total number of routes
    pop <- gl(nPops, routePerPop, nRoutes) # Population index for each route
  }else{
    nRoutes <- sum(routePerPop)
    pop <- as.factor(rep(seq(1:nPops), routePerPop)) # Population index for each route
  }

  if(length(alphaPop) == 1) {
    alphaPop <- rep(alphaPop, nPops)
  }


  # Generate data structure to hold counts and log (lambda)
  C <- log.expectedCount <- array(NA, dim = c(nYears, nRoutes))

  # Generate covariate values
  year <- 1:nYears
  yr <- (year - (nYears/2))/(nYears/2) # Standardize

  # Draw two sets of random effects from their respective distributions
  epsRoute <- rnorm(n = nRoutes, mean = 0, sd = sdRoute)
  epsYear <- rnorm(n = nYears, mean = 0, sd = sdYear)

  # Loop over routes
  for (i in 1:nRoutes){

    # Build up systematic part of the GLM including random effects
    log.expectedCount[,i] <- alphaPop[pop[i]] + beta*yr + epsRoute[i] +
      epsYear
    expectedCount <- exp(log.expectedCount[,i])

    C[,i] <- rpois(n = nYears, lambda = expectedCount)
  }

  return(list(nPops = nPops, nRoutes = nRoutes, nYears = nYears,
              routePerPop = routePerPop, year = yr, pop = pop,
              alphaPop = alphaPop, epsRoute = epsRoute,
              epsYear = epsYear, beta = beta,
              sdRoute = sdRoute, sdYear = sdYear,
              expectedCount = expectedCount, C = C))
}



###############################################################################
# Estimates population-level relative abundance from count data
###############################################################################
#' Estimates population-level relative abundance from count data
#'
#' Uses a Bayesian heirarchical model to estimate relative abundance of regional
#' populations from count-based data (e.g., Breeding Bird Survey)
#'
#' @param count_data List containing the following elements:
#' ' \describe{
#'    \item{\code{C}}{nYears by nRoutes matrix containing the observed number of individuals counted at each route in each year.}
#'    \item{\code{pop}}{Vector of length nRoutes indicating the population/region in which each route is located.}
#'    \item{\code{routePerPop}}{Vector of length 1 or nPops containing the number of routes (i.e. counts) per population. If length(routePerPop) == 1, number of routes is identical for each population.}
#' }
#' @param ni Number of MCMC iterations. Default = 20000.
#' @param nt Thinning rate. Default = 5.
#' @param nb Number of MCMC iterations to discard as burn-in. Default = 5000.
#' @param nc Number of chains. Default = 3.
#'
#' @return \code{modelCountDataJAGS} returns an mcmc object containing posterior samples for each monitored parameter.
#
#'
#' @export
#' @example inst/examples/simCountExamples.R
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513-524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
#'
#' Link, W. A. and J. R. Sauer. 2002. A hierarchical analysis of population
#' change with application to Cerulean Warblers. Ecology 83: 2832-2840.

modelCountDataJAGS <- function (count_data, ni = 20000, nt = 5, nb = 5000, nc = 3) {
  nPops <- length(unique(count_data$pop))
  nRoutes <- dim(count_data$C)[2]
  nYears = dim(count_data$C)[1]
  if(length(count_data$routePerPop) == 1){
    routePerPop = rep(count_data$routePerPop, nPops)
  } else {
    routePerPop = count_data$routePerPop
  }

  # Initial values
  jags.inits <- function()list(mu = runif(1,0,2), alpha = runif(nPops, -1,1), beta1 = runif(1,-1,1),
                               tau.alpha = runif(1,0,0.1), tau.noise = runif(1,0,0.1),
                               tau.rte = runif(1,0,0.1), route = runif(nRoutes,-1,1))
  # Parameters to monitor
  params <- c("mu", "alpha", "beta1", "sd.alpha", "sd.rte", "sd.noise", "totalN", "popN", "relN")

  # Data
  jags.data <- list(C = count_data$C, nPops = length(unique(count_data$pop)), nRoutes = nRoutes,
                    routePerPop = routePerPop,
                    year = seq(from = 0, to = 1, length.out = nYears), nYears = nYears, pop = count_data$pop)


  out <- R2jags::jags(data = jags.data, inits = jags.inits, params,
                      paste0(find.package('MigConnectivity'), "/JAGS/sim_Poisson2.txt"),
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      progress.bar = 'none')

  return(coda::as.mcmc(out))
}
