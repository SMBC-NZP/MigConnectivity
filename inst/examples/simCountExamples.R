set.seed(150)

### Set parameters for simulation ----

# Number of populations
nStrata. <- 4
# Number of routes w/i each population (assumed to be balanced)
routePerStrat. <- 30 # reduced from 90 for example speed
# Number of years
nYears. <- 5 # reduced from 10 for example speed
# log(Expected number of birds counted at each route)
alphaStrat. <- 1.95
# standard deviation of normal distribution assumed for route/observer random
# effects
sdRoute. <- 0.6
# standard deviation of normal distribution assumed for year random effects
sdYear. <- 0.18


# Number of simulated datasets to create and model
nsims <- 50 # reduced from 100 for example speed
# Number of MCMC iterations
ni. <- 1000 # reduced from 15000 for example speed
# Number of iterations to thin from posterior (reduced from 5)
nt. <- 1
# Number of iterations to discard as burn-in
nb. <- 500 # reduced from 5000 for example speed
# Number of MCMC chains
nc. <- 1 # reduced from 3 for example speed

### Create empty matrix to store model output ---
sim_in <- vector("list", nsims)
sim_out <- vector("list", nsims)


# Simulation ---
\donttest{
system.time(for(s in 1:nsims){
  cat("Simulation",s,"of",nsims,"\n")

  # Simulate data
  sim_data <- simCountData(nStrata = nStrata., routesPerStrata = routePerStrat.,
                           nYears = nYears., alphaStrat = alphaStrat.,
                           sdRoute = sdRoute., sdYear = sdYear.)
  sim_in[[s]] <- sim_data


  # Estimate population-level abundance
  out_mcmc <- modelCountDataJAGS(count_data = sim_data, ni = ni., nt = nt.,
                                 nb = nb., nc = nc.)

  # Store model output
  sim_out[[s]] <- out_mcmc
  remove(out_mcmc)

})


### Check that relative abundance is, on average, equal for each population
prop.table(sapply(sim_in, function(x) return(rowsum(colSums(x$C), x$strat))), 2)

rel_names <- paste0('relN[', 1:nStrata., ']')
rel_abund1 <- data.frame(sim=1:nsims,
                         ra1.mean=NA, ra2.mean=NA, ra3.mean=NA, ra4.mean=NA,
                         ra1.low=NA, ra2.low=NA, ra3.low=NA, ra4.low=NA,
                         ra1.high=NA, ra2.high=NA, ra3.high=NA, ra4.high=NA,
                         ra1.cover=0, ra2.cover=0, ra3.cover=0, ra4.cover=0)
for (s in 1:nsims) {
  rel_abund1[s, 2:5] <- summary(sim_out[[s]])$statistics[rel_names, "Mean"]
  rel_abund1[s, 6:9] <- summary(sim_out[[s]])$quantiles[rel_names, 1]
  rel_abund1[s, 10:13] <- summary(sim_out[[s]])$quantiles[rel_names, 5]
}
rel_abund1 <- transform(rel_abund1,
                        ra1.cover = (ra1.low<=0.25 & ra1.high>=0.25),
                        ra2.cover = (ra2.low<=0.25 & ra2.high>=0.25),
                        ra3.cover = (ra3.low<=0.25 & ra3.high>=0.25),
                        ra4.cover = (ra4.low<=0.25 & ra4.high>=0.25))

summary(rel_abund1)
}

