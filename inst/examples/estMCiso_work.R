library(raster); library(MigConnectivity)

OVENdist <- shapefile("data-raw/Spatial_Layers/OVENdist.shp")
OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

OVENvals <- read.csv("data-raw/deltaDvalues.csv")

OVENvals <- OVENvals[grep(OVENvals[,1],pattern = "NH", invert = TRUE),]

nSamples <- 1000
nAnimals <- nrow(OVENvals)

set.seed(122)
a <- Sys.time()
b <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = FALSE,
               odds = NULL,
               SingleCellAssign = TRUE,
               nSamples = nSamples,
               dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "sim.cell",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = TRUE)
Sys.time()-a

str(b)

OVENdist <- rgeos::gUnaryUnion(OVENdist, id = OVENdist$ORIGIN)

r <- raster(extent(OVENdist)+10,res = c(40,10))

r1 <- raster::rasterize(OVENdist,r,mask = TRUE)

r1[] <- 1:ncell(r1)

r2 <- rasterToPolygons(r1)

targetSites <- raster::intersect(OVENdist,r2)
targetSites$layer <- 1:length(targetSites)
targetSites <- as(targetSites,"SpatialPolygons")

targetDist <- distFromPos(rgeos::gCentroid(targetSites,byid = TRUE)@coords)

plot(targetSites)
plot(SpatialPoints(t(b[1,,])),add = TRUE, pch = 19)

results <- array(NA, c(nAnimals, nSamples))
for(i in 1:nSamples) {
  results[, i] <- over(SpatialPoints(t(b[i, , ]),
                                     proj4string = CRS(targetSites@proj4string@projargs)),
                       targetSites)
}

# What is proportion of simulated points not in targetSites, with restrict2Likely TRUE?
sum(is.na(results)) / length(results)
# [1] 0.07018543

set.seed(122)
a <- Sys.time()
b <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = FALSE,
               odds = NULL,
               SingleCellAssign = TRUE,
               nSamples = nSamples,
               dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "sim.cell",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual",
               restrict2Likely = FALSE)
Sys.time()-a

results <- array(NA, c(nAnimals, nSamples))
for(i in 1:nSamples) {
  results[, i] <- over(SpatialPoints(t(b[i, , ]),
                                     proj4string = CRS(targetSites@proj4string@projargs)),
                       targetSites)
}

# What is proportion of simulated points not in targetSites, with restrict2Likely FALSE?
sum(is.na(results)) / length(results)
# Same:
# [1] 0.07018543

# Plot points that aren't in targetSites only
notIn <- which(is.na(t(results)), arr.ind = T)
plot(targetSites)
for (i in 1:nrow(notIn))
  plot(SpatialPoints(t(b[notIn[i,1],,notIn[i,2]])),add = TRUE, pch = 19)

Countries <- shapefile("data-raw/Spatial_Layers/TM_WORLD_BORDERS-0.3.shp")
JAM <- Countries[Countries$NAME == "Jamaica",]
JAM <- as(JAM,"SpatialPolygons")

States <- shapefile("data-raw/Spatial_Layers/st99_d00.shp")
FL <- States[States$NAME == "Florida",]
FL <- rgeos::gUnaryUnion(FL,id = FL$STATE)

originSites <- rbind(JAM,FL)

originDist <- distFromPos(rgeos::gCentroid(originSites,byid = TRUE)@coords)

ever <- length(grep(OVENvals[,1],pattern = "EVER")) / nAnimals
jam <- length(grep(OVENvals[,1],pattern = "JAM")) / nAnimals


originRelAbundance <- as.matrix(c(jam, ever))

targetSites <- sp::spTransform(targetSites, CRS(MigConnectivity::projections$Lambert))
originSites <- sp::spTransform(originSites, CRS(MigConnectivity::projections$Lambert))

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
            targetPoints = b,
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
            approxSigTest = F,
            sigConst = 0,
            resampleProjection = MigConnectivity::projections$EquidistConic,
            maxTries = 300,
            isIntrinsic = TRUE))

str(MC)

(mcDiff <- diffMC(list(OVEN1 = Combined, OVEN2 = MC)))
(rMDiff <- diffMantel(list(OVEN1 = Combined, OVEN2 = MC)))

object.size(b)
bp <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = FALSE,
               odds = NULL,
               SingleCellAssign = F,
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
               SingleCellAssign = F,
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
               SingleCellAssign = F,
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
               SingleCellAssign = F,
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
               SingleCellAssign = F,
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
               SingleCellAssign = F,
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
