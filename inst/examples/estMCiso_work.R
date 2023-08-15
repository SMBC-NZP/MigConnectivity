library(raster); library(MigConnectivity)

# read in a shapefile with sf package
OVENdist <- sf::st_read("data-raw/Spatial_Layers/OVENdist.shp")

# keep ony the breeding info
OVENdist <- OVENdist[OVENdist$ORIGIN==2,]

# set the crs to WGS84
sf::st_crs(OVENdist) <- 4326

# read in data
OVENvals <- read.csv("data-raw/deltaDvalues.csv")

# save only those from outside NH
OVENvals <- OVENvals[grep(OVENvals[,1],pattern = "NH", invert = TRUE),]

nSamples <- 1000
nAnimals <- nrow(OVENvals)

set.seed(122)
a <- Sys.time()
b <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               #oddsRatio = FALSE,
               #odds = NULL,
               #SingleCellAssign = TRUE,
               nSamples = nSamples,
               #dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               #return = "sim.cell",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
Sys.time()-a

summary(b)

# OVENdist <- rgeos::gUnaryUnion(OVENdist, id = OVENdist$ORIGIN)
OVENdist <- aggregate(OVENdist, list(OVENdist$ORIGIN), head, n=1)

r <- raster(extent(OVENdist)+10,res = c(40,10))

r1 <- raster::rasterize(OVENdist,r,mask = TRUE)

r1[] <- 1:ncell(r1)

r2 <- sf::st_as_sf(rasterToPolygons(r1))
sf::st_crs(r2) <- 4326

targetSites <- raster::intersect(OVENdist,r2)
targetSites$layer <- 1:nrow(targetSites)


targetDist <- distFromPos(sf::st_coordinates(sf::st_centroid(targetSites)))

plot(sf::st_geometry(targetSites))
plot(sp::SpatialPoints(t(b$SingleCell[1,,])),add = TRUE, pch = 19)

results <- array(NA, c(nAnimals, nSamples))
for(i in 1:nSamples) {
  point.samples <- sf::st_as_sf(data.frame(t(b$SingleCell[i,,])),
                                coords = c("Longitude","Latitude"),
                                crs = 4326)
  results[, i] <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = point.samples,
                                                      y = targetSites,
                                                      sparse = TRUE))))
}
# What is proportion of simulated points not in targetSites, with restrict2Likely TRUE?
sum(is.na(results)) / length(results)
# 0

set.seed(122)
a <- Sys.time()
b <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               #oddsRatio = FALSE,
               #odds = NULL,
               #SingleCellAssign = TRUE,
               nSamples = nSamples,
               #dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               #return = "sim.cell",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = FALSE)
Sys.time()-a

results <- array(NA, c(nAnimals, nSamples))
for(i in 1:nSamples) {
  point.samples <- sf::st_as_sf(data.frame(t(b$SingleCell[i,,])),
                                coords = c("Longitude","Latitude"),
                                crs = 4326)
  results[, i] <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = point.samples,
                                                                        y = targetSites,
                                                                        sparse = TRUE))))
}

# What is proportion of simulated points not in targetSites, with restrict2Likely FALSE?
sum(is.na(results)) / length(results)
# Same:
# 0

# Plot points that aren't in targetSites only
# notIn <- which(is.na(t(results)), arr.ind = T)
# plot(targetSites)
# for (i in 1:nrow(notIn))
#  plot(SpatialPoints(t(b[notIn[i,1],,notIn[i,2]])),add = TRUE, pch = 19)

Countries <- sf::st_read("data-raw/Spatial_Layers/TM_WORLD_BORDERS-0.3.shp")
JAM <- Countries[Countries$NAME == "Jamaica",]
# JAM <- as(JAM,"SpatialPolygons")

States <- sf::st_read("data-raw/Spatial_Layers/st99_d00.shp")
FL <- States[States$NAME == "Florida",]
FL <- aggregate(FL, list(FL$STATE), head, n=1)


originSites <- rbind(JAM[,c("NAME","geometry")],
                     FL[,c("NAME","geometry")])


originDist <- distFromPos(sf::st_coordinates(sf::st_centroid(originSites)))

ever <- length(grep(OVENvals[,1],pattern = "EVER")) / nAnimals
jam <- length(grep(OVENvals[,1],pattern = "JAM")) / nAnimals


originRelAbundance <- as.matrix(c(jam, ever))



# # # # # # # # # # # # # # # # # # # #
# NEEDS CONVERSION from sp to sf BELOW



targetSites <- sf::st_transform(targetSites, srid = 9820)#crs = CRS(MigConnectivity::projections$Lambert)
originSites <- sf::st_transform(originSites, 9820)#CRS(MigConnectivity::projections$Lambert))

set.seed(12)
originLongLat <- array(NA,c(nrow(OVENvals),2))
sample1 <- sp::spsample(originSites[1,],length(grep(OVENvals[,1],pattern = "JAM")),type = "random")
sample2 <- sp::spsample(originSites[2,],length(grep(OVENvals[,1],pattern = "EVER")),type = "random")
originLongLat[grep(OVENvals[,1],pattern = "JAM"),]<-sample1@coords
originLongLat[grep(OVENvals[,1],pattern = "EVER"),]<-sample2@coords

originPoints <- sp::SpatialPoints(originLongLat, proj4string = sp::CRS(originSites@proj4string@projargs))

targetSites <- sp::spTransform(targetSites,sp::CRS(MigConnectivity::projections$WGS84))
originSites <- sp::spTransform(originSites,sp::CRS(MigConnectivity::projections$WGS84))
originPoints <- sp::spTransform(originPoints,sp::CRS(MigConnectivity::projections$WGS84))

#library(MigConnectivity)

system.time(MC <- estMC(targetDist = targetDist,
            originRelAbund = originRelAbundance,
            targetIntrinsic = b,
            targetSites = targetSites,
            targetAssignment = NULL,
            originDist = originDist,
            originPoints = originPoints,
            originSites = originSites,
            originAssignment = NULL,
            originNames = NULL,
            targetNames = NULL,
            nSamples = 100,
            verbose = 2,
            nSim = 5,
            calcCorr = TRUE,
            alpha = 0.05,
            approxSigTest = FALSE,
            sigConst = 0,
            resampleProjection = MigConnectivity::projections$EquidistConic,
            maxTries = 300,
            isIntrinsic = TRUE))

str(MC)
summary(MC, digits = 2)

(mcDiff <- diffMC(list(OVEN1 = Combined, OVEN2 = MC)))
(rMDiff <- diffMantel(list(OVEN1 = Combined, OVEN2 = MC)))

object.size(b)
bp <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = FALSE,
               odds = NULL,
               SingleCellAssign = FALSE,
               nSamples = nSamples,
               dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "probability",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
object.size(bp)
bpd <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = FALSE,
               odds = NULL,
               SingleCellAssign = FALSE,
               nSamples = nSamples,
               dataFrame = TRUE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "probability",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
object.size(bpd)
bo <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = TRUE,
               odds = NULL,
               SingleCellAssign = FALSE,
               nSamples = nSamples,
               dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "odds",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
object.size(bo)
bod <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = TRUE,
               odds = NULL,
               SingleCellAssign = FALSE,
               nSamples = nSamples,
               dataFrame = TRUE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "odds",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
object.size(bod)
pop <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = TRUE,
               odds = NULL,
               SingleCellAssign = FALSE,
               nSamples = nSamples,
               dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "population",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
object.size(pop)
popd <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = TRUE,
               odds = NULL,
               SingleCellAssign = FALSE,
               nSamples = nSamples,
               dataFrame = TRUE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "population",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
object.size(popd)
print(object.size(b), units = "Mb")
print(object.size(bp), units = "Mb")
print(object.size(pop), units = "Mb")
print(object.size(bpd), units = "Mb")
print(object.size(popd), units = "Mb")
print(object.size(b) + object.size(bp) + object.size(bo) + object.size(pop), units = "Mb")
print(object.size(b) + object.size(bpd) + object.size(bod) + object.size(popd), units = "Mb")
print(object.size(b) + object.size(bp) + object.size(bo) + object.size(pop) +
  object.size(bpd) + object.size(bod) + object.size(popd), units = "Mb")
