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

#' Simulate geolocator data
#'
#' @param psi Transition probabilities between B origin and W target sites.
#'    A matrix with B rows and W columns where rows sum to 1.
#' @param originRelAbund Relative abundances at B origin sites. Numeric vector
#'    of length B that sums to 1.
#' @param sampleSize List of length two. The first element is either a vector of
#'    length B with the number of simulated animals to release with geolocators
#'    at each of the B origin sites, a single integer with the total number of
#'    simulated animals to release with geolocators at origin sites (in which
#'    case, the origin sites will be sampled according to the relative
#'    abundance), or NULL if all animals are released at target sites. The
#'    second element is either a vector of length W with the number of simulated
#'    animals to release with geolocators at each of the W target sites, a
#'    single integer with the total number of simulated animals to release with
#'    geolocators at target sites (in which case, the target sites will be
#'    sampled according to their relative abundance), or NULL if all animals are
#'    released at origin sites.
#' @param originSites A polygon spatial layer (sf - MULTIPOLYGON) defining the
#'    geographic representation of sites in the origin season.
#' @param targetSites A polygon spatial layer (sf - MULTIPOLYGON) defining the
#'    geographic representation of sites in the target season.
# @param captured
#' @param geoBias Vector of length 2 indicating expected bias in longitude and
#'    latitude of animals captured and released at origin sites, in
#'    \code{targetSites} units.
#' @param geoVCov 2x2 matrix with expected variance/covariance in longitude and
#'    latitude of animals captured and released at origin sites, in
#'    \code{targetSites} units.
#' @param geoBiasOrigin Vector of length 2 indicating expected bias in longitude
#'    and latitude of animals captured and released at target sites, in
#'    \code{originSites} units.
#' @param geoVCovOrigin 2x2 matrix with expected variance/covariance in
#'    longitude and latitude of animals captured and released at target sites,
#'    in \code{originSites} units.
#' @param S Survival probabilities of released geolocator animals. Either a
#'    matrix with B rows and W columns (if survival depends on both origin site
#'    and target site), a vector of length W (if survival depends only on target
#'    site), or a single number (if survival is the same for all animals).
#'    Default 1 (all animals with geolocators survive a year).
#' @param p Recapture probabilities of released geolocator animals; list of
#'    length two. The first element is either a vector of length B (if recapture
#'    depends on origin site), or a single number (if recapture is the same for
#'    all animals released on origin sites). The second element is either a
#'    vector of length W (if recapture depends on target site), or a single
#'    number (if recapture is the same for all animals released on target
#'    sites). Default list(1, 1) (all animals that survive are recaptured).
#' @param requireEveryOrigin if TRUE, the function will throw an error if it
#'    looks like at least one origin site has no animals released in or
#'    migrating to it, or if it can, keep simulating until representation is
#'    met. This helps estTransition or estMC not throw an error. Default FALSE.
#' @return
#' @export
#'
# @examples
simGLData <- function(psi, originRelAbund, sampleSize,
                  originSites = NULL, targetSites = NULL,
                  #captured = "origin",
                  geoBias = NULL, geoVCov = NULL,
                  geoBiasOrigin=NULL, geoVCovOrigin=NULL,
                  S = 1, p = list(1, 1),
                  requireEveryOrigin = FALSE,
                  verbose = 0) {
  nOriginSites <- nrow(psi)
  nTargetSites <- ncol(psi)
  rev <- MigConnectivity:::reversePsiRelAbund(psi, originRelAbund)
  gamma <- rev$gamma
  targetRelAbund <- rev$targetRelAbund
  if (length(dim(S))<2) {
    S <- matrix(S, nOriginSites, nTargetSites, byrow = T)
  }
  if (length(p[[1]]==1))
    p[[1]] <- rep(p[[1]], nOriginSites)
  if (length(p[[2]]==1))
    p[[2]] <- rep(p[[2]], nTargetSites)
  nAnimals <- sum(sampleSize[[1]]) + sum(sampleSize[[2]])
  captured <- rep(c("origin", "target"),
                  c(sum(sampleSize[[1]]), sum(sampleSize[[2]])))
  # if (length(captured)==1)
  #   captured <- rep(captured, nAnimals)
  targetAssignment <- rep(NA, nAnimals)
  originAssignment <- rep(NA, nAnimals)
  lived <- rep(1, nAnimals)
  recaptured <- rep(1, nAnimals)
  if (length(sampleSize[[1]])>1){
    originAssignment[captured=="origin"] <- rep(1:nOriginSites, sampleSize[[1]])
    if (requireEveryOrigin && is.null(sampleSize[[2]])){
      if (any(sampleSize[[1]] < 1))
        stop("Some origin site or sites have no samples ", sampleSize[[1]])
    }
  }
  if (length(sampleSize[[2]])>1){
    targetAssignment[captured=="target"] <- rep(1:nTargetSites, sampleSize[[2]])
  }

  runWell <- FALSE
  while (!runWell){
    for (a in 1:nAnimals) {
      if (captured[a]=="origin") {
        if (is.na(originAssignment[a]))
          originAssignment[a] <- sample.int(nOriginSites, 1, prob = originRelAbund)
        targetAssignment[a] <- sample.int(nTargetSites, 1, prob = psi[originAssignment[a], ])
        lived[a] <- rbinom(1, 1, S[originAssignment[a], targetAssignment[a]])
        recaptured[a] <- rbinom(1, 1, lived[a] * p[[1]][originAssignment[a]])
      }
      else {
        if (is.na(targetAssignment[a]))
          targetAssignment[a] <- sample.int(nTargetSites, 1, prob = targetRelAbund)
        originAssignment[a] <- sample.int(nOriginSites, 1, prob = gamma[targetAssignment[a], ])
        lived[a] <- rbinom(1, 1, S[originAssignment[a], targetAssignment[a]])
        recaptured[a] <- rbinom(1, 1, lived[a] * p[[2]][targetAssignment[a]])
      }
    }
    oa1 <- sort(unique(originAssignment))
    runWell <- !requireEveryOrigin || any(captured=="target") ||
      isTRUE(all.equal(oa1, 1:nOriginSites))
    if (!runWell && verbose > 0)
      cat(oa1, "\n")
  }


  oaTest <- list(1:2); taTest <- list(2:3)
  while (any(lengths(oaTest)>1) || any(lengths(taTest)>1)) {
    originPointsTrue <- randomPoints(originSites, originAssignment)
    targetPointsTrue <- randomPoints(targetSites, targetAssignment)
    originPointsObs <- array(NA, c(sum(recaptured), 2),
                             dimnames = list(NULL, c("x","y")))
    targetPointsObs <- array(NA, c(sum(recaptured), 2),
                             dimnames = list(NULL, c("x","y")))
    linkerOrigin <- which(recaptured == 1) %in%
      which(captured=="origin" & recaptured == 1)
    linkerTarget <- which(recaptured == 1) %in%
      which(captured=="target" & recaptured == 1)
    if (any(captured=="origin" & recaptured == 1)) {
      geoBias2 <- array(rep(geoBias, sum(captured=="origin" & recaptured == 1)),
                        c(2, sum(captured=="origin" & recaptured == 1)))
      targetPointsObs[linkerOrigin, ] <-
        t(array(apply(sf::st_coordinates(targetPointsTrue)[captured=="origin" & recaptured == 1, , drop = FALSE], 1,
                      MASS::mvrnorm, n=1, Sigma=geoVCov),
                c(2, sum(captured=="origin" & recaptured == 1))) + geoBias2)
      originPointsObs[linkerOrigin, ] <-
        sf::st_coordinates(originPointsTrue)[captured=="origin" & recaptured == 1, ,
                                             drop = FALSE]
    }
    if (any(captured == "target" & recaptured == 1)) {
      geoBias2 <- array(rep(geoBiasOrigin, sum(captured == "target" & recaptured == 1)),
                        c(2, sum(captured == "target" & recaptured == 1)))
      originPointsObs[linkerTarget, ] <-
        t(array(apply(sf::st_coordinates(originPointsTrue)[captured == "target" & recaptured == 1, , drop = FALSE], 1,
                      MASS::mvrnorm, n=1, Sigma=geoVCovOrigin),
                c(2, sum(captured == "target" & recaptured == 1))) + geoBias2)
      targetPointsObs[linkerTarget, ] <-
        sf::st_coordinates(targetPointsTrue)[captured=="target" & recaptured == 1, ,
                                             drop = FALSE]
    }
    targetPointsObs <- sf::st_as_sf(data.frame(targetPointsObs),
                                    coords = c("x","y"),
                                    crs = sf::st_crs(targetSites))
    originPointsObs <- sf::st_as_sf(data.frame(originPointsObs),
                                    coords = c("x","y"),
                                    crs = sf::st_crs(originSites))
    oaTest <- suppressMessages(unclass(sf::st_intersects(x = originPointsObs,
                                                         y = originSites,
                                                         sparse = TRUE)))
    taTest <- suppressMessages(unclass(sf::st_intersects(x = targetPointsObs,
                                                         y = targetSites,
                                                         sparse = TRUE)))
  }
  return(list(originAssignment = originAssignment,
              targetAssignment = targetAssignment,
              originPointsTrue = originPointsTrue,
              targetPointsTrue = targetPointsTrue,
              originPointsObs = originPointsObs,
              targetPointsObs = targetPointsObs,
              lived = lived,
              recaptured = recaptured,
              input = list(psi = psi, originRelAbund = originRelAbund,
                           originSites = originSites,
                           targetSites = targetSites,
                           captured = captured,
                           geoBias = geoBias, geoVCov = geoVCov,
                           geoBiasOrigin = geoBiasOrigin,
                           geoVCovOrigin = geoVCovOrigin,
                           S = S, p = p)))
}



###############################################################################
# Simulate genetic populations
###############################################################################
#' Simulates genetic population rasters from input polygons
#'
#' This code simulates rasters analogous to the genetic surfaces created from
#' the bird genoscape project
#'
#' @param popBoundaries - a list containing polygons that represent the
#'   boundaries of each region
#' @param npts - number of random points used to generate a kernel density
#'   surface within each region
#' @param res - the desired resolution of the output raster
#' @param bufferRegions boolean if \code{TRUE} a buffer is applied to each
#'   region to avoid or increase overlap
#' @param bufferDist - Desired buffer distance in meters. Positive values
#'   enlarge regions, negative buffers make regions smaller. Only used when
#'   \code{bufferRegions=TRUE}.
#' @param popNames - a vector the same length as popBoundaries
#'
#' @return returns a rasterStack probability surface.
#' @export

simGeneticPops <- function(popBoundaries,
                           npts = 1000,
						               res = NULL, #change to match isotope resolution
						               bufferRegions = FALSE,
						               bufferDist = -50000,
						               popNames = NULL,
						               verbose = 0){

  if(is.null(popNames)){
    popNames <- paste0("pop.",1:length(popBoundaries))
  }

  # merge the polygons to make an empty raster for the state-space #

  if(bufferRegions){
    if (verbose > 0)
      cat('Creating buffers ... \n')
    popBoundaries <- lapply(popBoundaries,
                            FUN = function(x){
                              origCRS <- sf::st_crs(x)
                              z <- sf::st_transform(x, "ESRI:102010")
                              z1 <- sf::st_buffer(z, bufferDist)
                              z2 <- sf::st_transform(z1, origCRS)
                              return(z2)})
  }

  if (verbose > 0)
    cat('Preparing data ... \n')

  crdsParams <- do.call(rbind,lapply(popBoundaries,sf::st_bbox))

  if(is.null(res)){
    res <- rep((crdsParams[,3]-crdsParams[,1])/20,2)
  }

  emptyRast <- raster::raster(xmn = min(crdsParams[,1]),
                              xmx = max(crdsParams[,3]),
                              ymn = min(crdsParams[,2]),
                              ymx = max(crdsParams[,4]),
                              res = res)
 # give raster same projection as it was made with
  popBcrs <- sf::st_crs(popBoundaries[[1]], parameters = TRUE)

  crs(emptyRast) <- sp::CRS(popBcrs$proj4string)

  if (verbose > 0)
    cat('Generating KDE ... \n')
  # No buffer applied to edges i.e., some probability overlap
  popKDE <- lapply(popBoundaries,
                   FUN = function(x){
                     # generate random points
                     z <- sf::st_sample(x,
                                        size = npts,
                                        type = "random")
                     z1 <- raster::raster(ks::kde(sf::st_coordinates(z),
                                          h = ks::Hlscv(sf::st_coordinates(z))))
                     # convert to probability
                     z1 <- z1/raster::cellStats(z1, stat = 'sum', na.rm = TRUE)
                     # resample to larger raster
                     y <- raster::resample(z1, emptyRast)
                     y[is.na(y)] <- 0
                     y[y<0] <- 0
                     return(y)})

  # convert the probability of assignment within each population to 1
  # i.e., they can be assigned anywhere within that population if
  # assigned there - following gaiah package

  popBinary <- lapply(popBoundaries,
                      FUN = function(x){
                      popRast <- raster::rasterize(as(x,'Spatial'),
                                           emptyRast,
                                           field = 1,
                                           background=0)
                      return(popRast)
                      })

  if (verbose > 0)
    cat('Stacking rasters ... \n')

  # Stack the rasters #
  genStack <- raster::stack(popKDE)
  popBinary <- raster::stack(popBinary)

  # assign crs to rasters
  crs(genStack) <- crs(popBinary) <- crs(emptyRast)

  # rename the stack to identify the groups
  names(genStack) <- popNames
  names(popBinary) <- popNames

  return(list(genStack = genStack,
              popRast = popBinary))
}

#' Simulate animal genoscape data for estimating migratory connectivity
#'
#' @param genPops
#' @param psi
#' @param originRelAbund
#' @param sampleSize
#' @param originSites
#' @param targetSites
#' @param captured
#' @param requireEveryOrigin
#'
#' @return
#' @export
#'
# @examples
simGeneticData <- function(genPops,
                           psi,
                           originRelAbund,
                           sampleSize,
                           originSites = NULL,
                           targetSites = NULL,
                           captured = "target",
                           requireEveryOrigin = FALSE,
                           verbose = 0) {
  if (verbose > 0)
    cat("Setting up...\n")
  nOriginSites <- nrow(psi)
  nTargetSites <- ncol(psi)
  rev <- MigConnectivity:::reversePsiRelAbund(psi, originRelAbund)
  gamma <- rev$gamma
  targetRelAbund <- rev$targetRelAbund
  nAnimals <- sum(sampleSize)
  if (length(captured)>1)
    captured <- captured[1]
  targetAssignment <- rep(NA, nAnimals)
  originAssignment <- rep(NA, nAnimals)
  if (captured=="origin" && length(sampleSize)>1){
    originAssignment <- rep(1:nOriginSites, sampleSize)
    if (requireEveryOrigin){
      if (any(sampleSize < 1))
        stop("Some origin site or sites have no samples ", sampleSize)
    }
  }
  if (captured=="target" && length(sampleSize)>1){
    targetAssignment <- rep(1:nTargetSites, sampleSize)
  }

  if (verbose > 0)
    cat("Assigning sites to animals...\n")
  runWell <- FALSE
  while (!runWell){
    for (a in 1:nAnimals) {
      if (captured=="origin") {
        if (is.na(originAssignment[a]))
          originAssignment[a] <- sample.int(nOriginSites, 1, prob = originRelAbund)
        targetAssignment[a] <- sample.int(nTargetSites, 1, prob = psi[originAssignment[a], ])
      }
      else {
        if (is.na(targetAssignment[a]))
          targetAssignment[a] <- sample.int(nTargetSites, 1, prob = targetRelAbund)
        originAssignment[a] <- sample.int(nOriginSites, 1, prob = gamma[targetAssignment[a], ])
      }
    }
    runWell <- !requireEveryOrigin || captured=="target" ||
      all.equal(unique(originAssignment), 1:nOriginSites)
  }
  # originAssignment <- sort(originAssignment)
  # targetAssignment <- sort(targetAssignment)

  if (verbose > 0)
    cat("Assigning points to animals...\n")
  originPointsTrue <- randomPoints(originSites, originAssignment)
  targetPointsTrue <- randomPoints(targetSites,targetAssignment)
  if (verbose > 0)
    cat("Assigning genProbs...\n")

  if (captured == "target") {

    genProbs <- raster::extract(genPops$genStack,
                                sf::st_coordinates(originPointsTrue))
    # convert from probability of being from location to probability of being
    # from population
    genProbs <- t(apply(genProbs, 1, FUN=function(x) x * originRelAbund /
        sum(x * originRelAbund)))
  }
  else {
    genProbs <- raster::extract(genPops$genStack,
                                sf::st_coordinates(targetPointsTrue))
    # convert from probability of being from location to probability of being
    # from population
    genProbs <- t(apply(genProbs, 1, FUN=function(x) x * targetRelAbund /
                          sum(x * targetRelAbund)))
  }

  # split out inds to list #
  if (verbose > 0)
    cat("Splitting genProbs to list...\n")
  indGP <- split(genProbs, f=1:nrow(genProbs))

  # use genProbs to create a raster probability surface
  if (verbose > 0)
    cat("Creating raster probability surface...\n")
  Ncell <- raster::cellStats(genPops$popRast, sum)

  indRasts <- stack(
              lapply(indGP,
                     function(x){
                      raster::calc(x*genPops$popRast / Ncell, fun = sum)
                     }))
  if (captured=="target"){
    tsCRS <- sf::st_crs(targetSites, parameters = TRUE)
    crs(indRasts) <- sp::CRS(tsCRS$proj4string) # sf::st_crs(targetSites)
  }else{
    osCRS <- sf::st_crs(originSites, parameters = TRUE)
    crs(indRasts) <- sp::CRS(osCRS$proj4string) #sf::st_crs(originSites)
    }
  return(list(originAssignment = originAssignment,
              targetAssignment = targetAssignment,
              originPointsTrue = originPointsTrue,
              targetPointsTrue = targetPointsTrue,
              genProbs = genProbs,
              genRaster = indRasts,
              input = list(genPops = genPops,
                           psi = psi,
                           originRelAbund = originRelAbund,
                           originSites = originSites,
                           targetSites = targetSites,
                           captured = captured)))
}
