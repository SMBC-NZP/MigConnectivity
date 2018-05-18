source("R/getIsoMap.R")
source("R/IsoAssign.R")
source("R/estConnectivity.R")
source("R/utilityFunctions.R")
library(raster); library(MigConnectivity)

OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

OVENvals <- read.csv("data-raw/deltaDvalues.csv")

OVENvals <- OVENvals[grep(OVENvals[,1],pattern = "NH", invert = TRUE),]

a <- Sys.time()
b <- isoAssign(isovalues = OVENvals[,2],
               isoSTD = 12,
               intercept = -10,
               slope = 0.8,
               oddsRatio = FALSE,
               odds = NULL,
               SingleCellAssign = TRUE,
               nSim = 1000,
               dataFrame = FALSE,
               sppShapefile = OVENdist,
               assignExtent = NULL,
               return = "sim.cell",
               element = "Hydrogen",
               surface = FALSE,
               period = "Annual")
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

plot(targetSites, add=T, col = "red")
plot(SpatialPoints(t(b[1,,])),add = TRUE, pch = 19)
plot(point.sample,add = TRUE, pch = 19)

results <- array(NA,c(151,1000))
for(i in 1:1000){
results[,i]<-over(SpatialPoints(t(b[i,,]),proj4string= CRS(targetSites@proj4string@projargs)),targetSites)
}

Countries <- shapefile("data-raw/Spatial_Layers/TM_WORLD_BORDERS-0.3.shp")
JAM <- Countries[Countries$NAME == "Jamaica",]
JAM <- as(JAM,"SpatialPolygons")

States <- shapefile("data-raw/Spatial_Layers/st99_d00.shp")
FL <- States[States$NAME == "Florida",]
FL <- rgeos::gUnaryUnion(FL,id = FL$STATE)

originSites <- rbind(JAM,FL)

originDist <- distFromPos(rgeos::gCentroid(originSites,byid = TRUE)@coords)

ever <- length(grep(OVENvals[,1],pattern = "EVER"))/nrow(OVENvals)
jam <- length(grep(OVENvals[,1],pattern = "JAM"))/nrow(OVENvals)


originRelAbundance <- as.matrix(c(jam,ever))

targetSites <- sp::spTransform(targetSites,CRS(MigConnectivity::projections$Lambert))
originSites <- sp::spTransform(originSites,CRS(MigConnectivity::projections$Lambert))

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

MC <- estMC(targetDist = targetDist,
            originRelAbund = originRelAbundance,
            targetPoints = b,
            targetSites = targetSites,
            #sampleSize = dim(b)[3],
            targetAssignment=NULL,
            originDist = originDist,
            originPoints=originPoints,
            originSites=originSites,
            originAssignment=NULL,
            originNames=NULL,
            targetNames=NULL,
            nSamples = 100,
            verbose = 2,
            nSim = 4,
            calcCorr=TRUE,
            alpha = 0.05,
            approxSigTest = F,
            sigConst = 0,
            resampleProjection = MigConnectivity::projections$EquidistConic,
            maxTries = 300,
            isIntrinsic = TRUE)

str(MC)

(mcDiff <- diffMC(list(OVEN1 = Combined, OVEN2 = MC)))
(rMDiff <- diffMantel(list(OVEN1 = Combined, OVEN2 = MC)))
