nSimulationsCMR <- 100   #
originPos13 <- matrix(c(rep(seq(-99, -81, 2), each = 10),
                          rep(seq(49, 31, -2), 10)), 100, 2)
targetPos13 <- matrix(c(rep(seq(-79, -61, 2), each = 10),
                           rep(seq(9, -9, -2), 10)), 100, 2)
originPosCMR <- rowsum(originPos13, c(rep(1:2, 5, each = 5),
                                         rep(3:4, 5, each = 5))) / 25
originPosCMR
targetPosCMR <- rowsum(targetPos13, c(rep(1:2, 5, each = 5),
                                           rep(3:4, 5, each = 5))) / 25
targetPosCMR

originDCMR <- distFromPos(originPosCMR, 'ellipsoid')
targetDCMR <- distFromPos(targetPosCMR, 'ellipsoid')
originRelAbundCMR <- rep(0.25, 4)
# the second intermediate psi scenario, the "low" level
psiMatCMR <- matrix(c(0.55, 0.2, 0.15, 0.1, 0.1, 0.55, 0.2, 0.15, 0.15, 0.1,
                     0.55, 0.2, 0.2, 0.15, 0.1, 0.55), 4, 4, byrow=T) #samplePsis[[5]]
MCCMR <- calcMC(originDCMR, targetDCMR, psiMatCMR, originRelNCMR)
MCCMR
# Storage list for examples
cmrExamples <- vector('list', nSimulationsCMR)
for (r in 1:nSimulationsCMR) {
  mod.name <- paste('psiB.enc2.band100', r, sep='.')
  in.file <- paste('out_', mod.name, '.gzip', sep='')   #inputs/files are from Cohen et al 2014 Ecology and Evolution, Simulation
  load(in.file)
  cmrExamples[[r]] <- get(mod.name)
  rm(list=mod.name)
}
devtools::use_data(cmrExamples, overwrite = T)

load("rel_abun_sim2.RData")   #Modeled seperately, see above
abundExamples <- sim_out
devtools::use_data(abundExamples, overwrite = T)
