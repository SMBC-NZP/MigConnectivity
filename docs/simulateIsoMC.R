## ----echo = FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, 
                      comment = "#>", 
                      tidy = FALSE,
                      message = FALSE,
                      warning = FALSE)

## ----eval = FALSE----------------------------------------------------------------------------
#  library(terra)
#  library(MigConnectivity)
#  library(sf)

## ----eval = FALSE----------------------------------------------------------------------------
#  # read in raster layer
#  download.file(
#  "https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity/master/data-raw/Spatial_Layers/bbsoven.txt",
#                destfile = "bbsoven.txt")
#  
#  OVENabund <- terra::rast("bbsoven.txt")
#  
#  OVENdist <- OVENabund
#  OVENdist[OVENdist>0]<-1
#  OVENdist[OVENdist==0]<-NA
#  
#  OVEN_single_poly <- terra::as.polygons(OVENdist)
#  terra::crs(OVEN_single_poly) <- MigConnectivity::projections$WGS84
#  terra::crs(OVENabund) <- MigConnectivity::projections$WGS84

## ----eval = FALSE----------------------------------------------------------------------------
#  ### get isotope raster for template
#  rasterTemplate <- getIsoMap()

## ----eval = FALSE----------------------------------------------------------------------------
#  terra::crs(rasterTemplate) <- MigConnectivity::projections$WGS84
#  rasterTemplate <- terra::crop(terra::mask(rasterTemplate,
#                                            OVEN_single_poly),
#                                OVEN_single_poly)
#  
#  # rasterize the distribution for relative abundance so that raster
#  # dimensions and resolution match the isotope layer
#  relativeAbund <- terra::project(OVENabund,rasterTemplate)
#  relativeAbund  <- relativeAbund /
#    unlist(c(terra::global(relativeAbund, fun = "sum", na.rm = T)))

## ----eval = FALSE----------------------------------------------------------------------------
#  # generate target sites
#  targetRanges <- vector('list',5)
#  # 3' latitude
#  targetRanges[[1]] <- terra::as.polygons(terra::rast(xmin = -180,
#                                           xmax = -40,
#                                           ymin = 25,
#                                           ymax = 85,
#                                           resolution = c(140,3)))
#  
#  # 5' latitude
#  targetRanges[[2]] <- terra::as.polygons(terra::rast(xmin = -180,
#                                               xmax = -40,
#                                               ymin = 25,
#                                               ymax = 85,
#                                               resolution = c(140,5)))
#  
#  # 10' latitude
#  targetRanges[[3]] <- terra::as.polygons(terra::rast(xmin = -180,
#                                           xmax = -40,
#                                           ymin = 25,
#                                           ymax = 85,
#                                           resolution = c(140,10)))
#  
#  # 12 isotope units
#  featherIso <- (0.95*getIsoMap(period = "GrowingSeason")+(-17.57))
#  iso <- terra::crop(featherIso, terra::ext(c(-180,-40,25,85)))
#  isocut <- terra::classify(iso,
#                            rcl = seq(terra::global(iso, fun = "min",
#                                                    na.rm = TRUE)$min,
#                                      terra::global(iso, fun = "max",
#                                                    na.rm = TRUE)$max,12))
#  targetRanges[[4]] <- terra::as.polygons(isocut)
#  
#  
#  # 12*2 isotope units
#  isocut <- terra::classify(iso,
#                            rcl = seq(terra::global(iso, fun = "min",
#                                                    na.rm = TRUE)$min,
#                                      terra::global(iso,fun = "max",
#                                                    na.rm = TRUE)$max, 24))
#  targetRanges[[5]] <- terra::as.polygons(isocut)
#  
#  
#  TR <- lapply(targetRanges, sf::st_as_sf)
#  vTR <- lapply(TR, sf::st_make_valid)
#  vTR <- lapply(vTR, FUN = function(x){sf::st_set_crs(x, 4326)})
#  
#  sf_use_s2(FALSE)
#  OVEN_single_poly <- sf::st_as_sf(OVEN_single_poly) %>% sf::st_make_valid()
#  
#  
#  # Keep only the targetSites that intersect with the OVEN polygon
#  targetRanges <- lapply(vTR, FUN = function(x){sf::st_intersection(x,OVEN_single_poly)})
#  
#  targetRanges <- lapply(targetRanges, FUN = function(x){y <- sf::st_as_sf(x)
#                                                         z <- sf::st_transform(y,4326)
#                                                         return(z)})
#  
#  for(i in 1:5){
#    targetRanges[[i]]$Target <- 1:nrow(targetRanges[[i]])
#  }
#  targetRanges <- lapply(targetRanges, sf::st_make_valid)
#  sf_use_s2(TRUE)

## ----eval = FALSE----------------------------------------------------------------------------
#  #Generate random breeding locations using the 10' target sites
#  Site1 <- st_as_sf(sf::st_sample(targetRanges[[3]][1:2,], size =100, type = "random"))
#  Site2 <- st_as_sf(sf::st_sample(targetRanges[[3]][2:3,], size =100, type = "random"))
#  Site3 <- st_as_sf(sf::st_sample(targetRanges[[3]][2:4,], size =100, type = "random"))
#  
#  # Capture coordinates
#  capCoords <- array(NA,c(3,2))
#  capCoords[1,] <- cbind(-98.17,28.76)
#  capCoords[2,] <- cbind(-93.70,29.77)
#  capCoords[3,] <- cbind(-85.000,29.836)
#  
#  featherIso <- (0.95*getIsoMap(period = "GrowingSeason")+(-17.57))
#  
#  # Extract simulated data
#  iso_dat <- data.frame(Site = rep(1:3, each = 100),
#                        xcoords = c(sf::st_coordinates(Site1)[,1],
#                                    sf::st_coordinates(Site2)[,1],
#                                    sf::st_coordinates(Site3)[,1]),
#                        ycoords = c(sf::st_coordinates(Site1)[,2],
#                                    sf::st_coordinates(Site2)[,2],
#                                    sf::st_coordinates(Site3)[,2]),
#                        targetSite = unlist(unclass(sf::st_intersects(rbind(Site1,
#                                                                            Site2,
#                                                                            Site3),
#                                          targetRanges[[3]]))),
#                        featherIso = terra::extract(featherIso,
#                                                    terra::vect(rbind(Site1,Site2,
#                                                                      Site3))))
#  
#  iso_dat <- iso_dat[complete.cases(iso_dat),]
#  
#  # generate transition data from simulation
#  sim_psi <- table(iso_dat$Site, iso_dat$targetSite)
#  
#  for(i in 1:nrow(sim_psi)){
#    sim_psi[i,]<-sim_psi[i,]/sum(sim_psi[i,])
#  }
#  
#  states <- rnaturalearth::ne_states(country = "United States of America",
#                                     returnclass = "sf")
#  originSites <- states[(states$woe_name %in% c("Texas","Louisiana","Florida")),]
#  originSites <- sf::st_transform(originSites, 4326)
#  
#  originDist <- distFromPos(st_coordinates(sf::st_centroid(originSites,
#                                                              byid = TRUE)))
#  targetDistMC <- distFromPos(sf::st_coordinates(sf::st_centroid(targetRanges[[3]],
#                                                                byid = TRUE)))

## ----eval = FALSE----------------------------------------------------------------------------
#  originRelAbund <- rep(1/3, 3)
#  nTargetSetups <- 5
#  nSims <- 2            #SET LOW FOR EXAMPLE
#  # nSims <- 200
#  nOriginSites = 3
#  
#  targetPoints0 <- matrix(NA, 300, 2, dimnames = list(NULL, c("x", "y")))
#  
#  
#  a <- Sys.time()
#  output.sims <- vector("list", nTargetSetups)
#  for(target in 1:nTargetSetups){
#    sim.output <- data.frame(targetSetup = rep(NA,nSims),
#                             sim = rep(NA,nSims),
#                             MC.generated = rep(NA,nSims),
#                             MC.realized = rep(NA,nSims),
#                             MC.est = rep(NA,nSims),
#                             MC.low = rep(NA,nSims),
#                             MC.high = rep(NA,nSims),
#                             rM.realized = rep(NA,nSims),
#                             rM.est = rep(NA,nSims),
#                             rM.low = rep(NA,nSims),
#                             rM.high = rep(NA,nSims))
#    set.seed(9001)
#    targetSites <- targetRanges[[target]]
#    sf_use_s2(FALSE)
#    targetSites <- sf::st_make_valid(targetSites)
#    targetDist <- distFromPos(sf::st_coordinates(sf::st_centroid(targetSites,
#                                                                 byid = TRUE)))
#  
#    # Extract simulated data
#    iso_dat <- data.frame(Site = rep(1:3,each = 100),
#                          xcoords = c(sf::st_coordinates(Site1)[,1],
#                                      sf::st_coordinates(Site2)[,1],
#                                      sf::st_coordinates(Site3)[,1]),
#                          ycoords = c(sf::st_coordinates(Site1)[,2],
#                                      sf::st_coordinates(Site2)[,2],
#                                      sf::st_coordinates(Site3)[,2]),
#                          targetSite = unlist(unclass(sf::st_intersects(rbind(Site1,
#                                                                              Site2,
#                                                                              Site3),
#                                          targetSites))),
#                          featherIso = terra::extract(featherIso,
#                                                      terra::vect(rbind(Site1,
#                                                                        Site2,
#                                                                        Site3))))
#  
#    iso_dat <- iso_dat[complete.cases(iso_dat),]
#  
#    # generate transition data from simulation
#    sim_psi <- prop.table(table(iso_dat$Site,factor(iso_dat$targetSite,
#                                                    1:nrow(targetSites))), 1)
#  
#    sf_use_s2(TRUE)
#    MC.generated <- calcMC(originDist = originDist,
#                            targetDist = targetDist,
#                            originRelAbund = originRelAbund,
#                            psi = sim_psi)
#  
#    for (sim in 1:nSims) {
#      cat("Simulation Run", sim, "of", nSims, "for target",target,"at", date(), "\n")
#      sim_move <- simMove(rep(100, nOriginSites), originDist, targetDist, sim_psi, 1, 1)
#      originAssignment <- sim_move$animalLoc[,1,1,1]
#      targetAssignment <- sim_move$animalLoc[,2,1,1]
#      sf_use_s2(FALSE)
#      for (i in 1:300) {
#        targetPoint1 <- sf::st_sample(targetSites[targetAssignment[i],],
#                                 size = 1, type = "random", iter = 25)
#        targetPoints0[i,] <- sf::st_coordinates(targetPoint1)
#      }
#      targetPoints <- sf::st_as_sf(as.data.frame(targetPoints0), crs = 4326,
#                                   coords = c("x", "y"))
#  
#      # Extract simulated data
#      iso_dat <- data.frame(Site = originAssignment,
#                            xcoords = targetPoints0[,1],
#                            ycoords = targetPoints0[,2],
#                            targetSite = unlist(unclass(sf::st_intersects(targetPoints,
#                                                                          targetSites))),
#                            featherIso = terra::extract(featherIso,
#                                                        terra::vect(targetPoints)))
#  
#  
#      iso_dat <- iso_dat[complete.cases(iso_dat),]
#  
#      # generate transition data from simulation
#      sim_psi0 <- table(iso_dat$Site,
#                        factor(iso_dat$targetSite,
#                               min(targetSites$Target):max(targetSites$Target)))
#      sim_psi_realized <- prop.table(sim_psi0, 1)
#  
#      # get points ready for analysis
#      nSite1 <- table(iso_dat$Site)[1]
#      nSite2 <- table(iso_dat$Site)[2]
#      nSite3 <- table(iso_dat$Site)[3]
#  
#      nTotal <- nSite1+nSite2+nSite3
#  
#      originCap <- array(NA, c(nTotal,2))
#  
#      wherecaught <- rep(paste0("Site", 1:3), c(nSite1, nSite2, nSite3))
#  
#      originCap[which(pmatch(wherecaught,"Site1",dup = TRUE)==1),1] <- capCoords[1,1]
#      originCap[which(pmatch(wherecaught,"Site1",dup = TRUE)==1),2] <- capCoords[1,2]
#  
#      originCap[which(pmatch(wherecaught,"Site2",dup = TRUE)==1),1] <- capCoords[2,1]
#      originCap[which(pmatch(wherecaught,"Site2",dup = TRUE)==1),2] <- capCoords[2,2]
#  
#      originCap[which(pmatch(wherecaught,"Site3",dup = TRUE)==1),1] <- capCoords[3,1]
#      originCap[which(pmatch(wherecaught,"Site3",dup = TRUE)==1),2] <- capCoords[3,2]
#  
#      originPoints <- sf::st_as_sf(as.data.frame(originCap), crs = 4326, coords = 1:2)
#      sf_use_s2(TRUE)
#  
#      MC.realized <- calcMC(originDist = originDist,
#                            targetDist = targetDist,
#                            originRelAbund = originRelAbund,
#                            psi = sim_psi_realized,
#                            sampleSize=nTotal)
#  
#  
#  
#      originPointDists <- distFromPos(originCap)
#      targetPointDists <- distFromPos(cbind(iso_dat$xcoords, iso_dat$ycoords))
#  
#  
#      simAssign <- isoAssign(isovalues = iso_dat$featherIso.GrowingSeasonD,
#                             isoSTD = 12,
#                             intercept = -17.57,
#                             slope = 0.95,
#                             odds = 0.67,
#                             restrict2Likely = FALSE,
#                             nSamples = 500,
#                             sppShapefile = OVEN_single_poly,
#                             relAbund = relativeAbund,
#                             verbose = 2,
#                             isoWeight = -0.7,
#                             abundWeight = 0,
#                             assignExtent = c(-179,-60,15,89),
#                             element = "Hydrogen",
#                             period = "GrowingSeason")
#      psi <- estTransition(targetRaster = simAssign,
#                          targetSites = targetSites,
#                          originPoints = originPoints,
#                          originSites = originSites, isRaster = TRUE,
#                          nSamples = 100, nSim = 5, verbose = 0)
#      rM <- estMantel(targetRaster = simAssign,
#                      targetSites = targetSites,
#                      originPoints = originPoints,
#                      isGL = FALSE, isTelemetry = FALSE,
#                      originSites = originSites, isRaster = TRUE,
#                      nBoot = 100, nSim = 5, verbose = 0, captured = "origin")
#      simEst <- estStrength(originRelAbund = originRelAbund,
#                            targetDist = targetDist,
#                            psi = psi,
#                            originDist = originDist,
#                            nSamples = 100,
#                            verbose = 0,
#                            alpha = 0.05)
#  
#      sim.output$targetSetup[sim] <- target
#      sim.output$sim[sim]<-sim
#      sim.output$MC.generated[sim] <- MC.generated
#      sim.output$MC.realized[sim] <- MC.realized
#      sim.output$MC.est[sim] <- simEst$MC$mean
#      sim.output$MC.low[sim] <- simEst$MC$bcCI[1]
#      sim.output$MC.high[sim] <- simEst$MC$bcCI[2]
#      sim.output$rM.realized[sim] <- calcMantel(originDist = originPointDists,
#                                                targetDist = targetPointDists)$pointCorr
#      sim.output$rM.est[sim] <- rM$corr$mean
#      sim.output$rM.low[sim] <- rM$corr$bcCI[1]
#      sim.output$rM.high[sim] <- rM$corr$bcCI[2]
#  
#    }
#    output.sims[[target]] <- sim.output
#  }
#  Sys.time()-a
#  
#  output.sims <- do.call("rbind", output.sims)
#  #

## ----eval = FALSE----------------------------------------------------------------------------
#  sim.output <- transform(output.sims,
#                          MC.generated.cover = as.integer((MC.low <= MC.generated &
#                                                             MC.high >= MC.generated)),
#                          MC.realized.cover = as.integer((MC.low <= MC.realized &
#                                                            MC.high >= MC.realized)),
#                          MC.generated.error = MC.est - MC.generated,
#                          MC.realized.error = MC.est - MC.realized,
#                          rM.cover = as.integer((rM.low <= rM.realized &
#                                                   rM.high >= rM.realized)),
#                          rM.error = rM.est - rM.realized)
#  
#  summary(sim.output)
#  # Examine results
#  aggregate(MC.generated.error ~ targetSetup, sim.output, mean)
#  aggregate(MC.generated.error ~ targetSetup, sim.output, function (x) mean(abs(x)))
#  aggregate(MC.generated.cover ~ targetSetup, sim.output, mean)
#  aggregate(MC.realized.error ~ targetSetup, sim.output, mean)
#  aggregate(MC.realized.error ~ targetSetup, sim.output, function (x) mean(abs(x)))
#  aggregate(MC.realized.cover ~ targetSetup, sim.output, mean)
#  aggregate(I(MC.realized - MC.generated) ~ targetSetup, sim.output, mean)
#  aggregate(rM.error ~ targetSetup, sim.output, mean)
#  aggregate(rM.error ~ targetSetup, sim.output, function (x) mean(abs(x)))
#  aggregate(rM.cover ~ targetSetup, sim.output, mean)
#  aggregate(MC.est ~ targetSetup, sim.output, mean)

