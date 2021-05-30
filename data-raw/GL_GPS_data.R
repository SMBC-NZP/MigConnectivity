################################################################################
#
#      Migratory Connectivity Metric - Cohen et al.
#
#      Geolocator data and GPS data from OVENBIRDS
#      Geolocator data - Hallworth et al. 2015 - Ecological Applications
#      GPS data - Hallworth and Marra 2015 Scientific Reports
#
#      Script written by M.T.Hallworth & J.A.Hostetler
################################################################################
# load required packages

library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(SpatialTools)
library(geosphere)
library(maptools)
library(shape)
library(ade4)
library(sf)

###################################################################
#
# geoVcov & geoBias
#
###################################################################

# WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# Lambert<-"+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# EquidistConic <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"

WGS84 <- "+init=epsg:4326"
Lambert <- 'PROJCS["North_America_Albers_Equal_Area_Conic",
            GEOGCS["GCS_North_American_1983",
            DATUM["North_American_Datum_1983",
            SPHEROID["GRS_1980",6378137,298.257222101]],
            PRIMEM["Greenwich",0],
            UNIT["Degree",0.017453292519943295]],
            PROJECTION["Albers_Conic_Equal_Area"],
            PARAMETER["False_Easting",0],
            PARAMETER["False_Northing",0],
            PARAMETER["longitude_of_center",-96],
            PARAMETER["Standard_Parallel_1",20],
            PARAMETER["Standard_Parallel_2",60],
            PARAMETER["latitude_of_center",40],
            UNIT["Meter",1],
            AUTHORITY["EPSG","102008"]]'

EquidistConic <- 'PROJCS["North_America_Equidistant_Conic",
                  GEOGCS["GCS_North_American_1983",
				  DATUM["North_American_Datum_1983",
				  SPHEROID["GRS_1980",6378137,298.257222101]],
				  PRIMEM["Greenwich",0],
				  UNIT["Degree",0.017453292519943295]],
				  PROJECTION["Equidistant_Conic"],
				  PARAMETER["False_Easting",0],
				  PARAMETER["False_Northing",0],
				  PARAMETER["Longitude_Of_Center",-96],
				  PARAMETER["Standard_Parallel_1",20],
				  PARAMETER["Standard_Parallel_2",60],
				  PARAMETER["Latitude_Of_Center",40],
				  UNIT["Meter",1],
				  AUTHORITY["EPSG","102010"]]'

# Define capture locations in the winter #

captureLocations<-matrix(c(-77.93,18.04,  # Jamaica
                           -80.94,25.13,  # Florida
                           -66.86,17.97,  # Puerto Rico
                           -71.72,43.95), # New Hampshire
                          nrow=4,ncol=2,byrow=TRUE)

# Convert capture locations into sf #
colnames(captureLocations) <- c("Longitude","Latitude")

CapLocs <- sf::st_as_sf(data.frame(captureLocations),
                        coords = c("Longitude","Latitude"),
                        crs = 4326)

# Project Capture locations #

CapLocsM<-sf::st_transform(CapLocs, 'ESRI:102010')

# Retrieve raw non-breeding locations from github #
# First grab the identity of the bird so we can loop through the files #
# For this example we are only interested in the error around non-breeding locations #
# here we grab only the birds captured during the non-breeding season #

winterBirds <- dget("https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity/master/data-raw/GL_NonBreedingFiles/winterBirds.txt")

# create empty list to store the location data #
Non_breeding_files <- vector('list',length(winterBirds))

# Get raw location data from Github #
for(i in 1:length(winterBirds)){
  Non_breeding_files[[i]] <- dget(paste0("https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity/master/data-raw/GL_NonBreedingFiles/NonBreeding_",winterBirds[i],".txt"))
}

# Remove locations around spring Equinox and potential migration points - same NB time frame as Hallworth et al. 2015 #
# two steps because subset on shapefile doesn't like it in a single step

Non_breeding_files <- lapply(Non_breeding_files,FUN = function(x){month <- as.numeric(format(x$Date,format = "%m"))
x[which(month != 3 & month != 4),]})


Jam <- c(1:9)   # locations within the list of winterBirds captured in Jamaica
Fla <- c(10:12) # locations within the list of winterBirds in Florida
PR <- c(13:16)  # locations within the list of winterBirds in Puerto Rico

# Turn the locations into shapefiles #

NB_GL <- lapply(Non_breeding_files,
                FUN = function(x){
                  sf::st_as_sf(x,
                               coords = c("Longitude","Latitude"),
                               crs =  4326)})

# Project into UTM projection #

NB_GLmeters <- lapply(NB_GL,
                      FUN = function(x){sf::st_transform(x,'ESRI:102010')})

# Process to determine geolocator bias and variance-covariance in meters #

# generate empty vector to store data #
# 16 birds were recovered during the non-breeding season
LongError <- rep(NA,length(winterBirds))
LatError <- rep(NA,length(winterBirds))

# Calculate the error in longitude derived from geolocators from the true capture location #
LongError[Jam] <- unlist(lapply(NB_GLmeters[Jam],
                                FUN = function(x){mean(sf::st_coordinates(x)[,1]-sf::st_coordinates(CapLocsM)[1,1])}))
LongError[Fla] <- unlist(lapply(NB_GLmeters[Fla],
                                FUN = function(x){mean(sf::st_coordinates(x)[,1]-sf::st_coordinates(CapLocsM)[2,1])}))
LongError[PR] <- unlist(lapply(NB_GLmeters[PR],
                               FUN = function(x){mean(sf::st_coordinates(x)[,1]-sf::st_coordinates(CapLocsM)[3,1])}))

# Calculate the error in latitude derived from geolocators from the true capture location #
LatError[Jam] <- unlist(lapply(NB_GLmeters[Jam],
                                FUN = function(x){mean(sf::st_coordinates(x)[,2]-sf::st_coordinates(CapLocsM)[1,2])}))
LatError[Fla] <- unlist(lapply(NB_GLmeters[Fla],
                                FUN = function(x){mean(sf::st_coordinates(x)[,2]-sf::st_coordinates(CapLocsM)[2,2])}))
LatError[PR] <- unlist(lapply(NB_GLmeters[PR],
                               FUN = function(x){mean(sf::st_coordinates(x)[,2]-sf::st_coordinates(CapLocsM)[3,2])}))
# Get co-variance matrix for error of known non-breeding deployment sites #

geo.error.model <- lm(cbind(LongError,LatError) ~ 1) # lm does multivariate normal models if you give it a matrix dependent variable!

geo.bias <- coef(geo.error.model)
geo.vcov <- vcov(geo.error.model)

###################################################################
#
#   Winter Locations - targetPoints
#     length = n animals tracked
#
###################################################################

#########################################################################################
#
# Here instead of using the raw points - use the KDE to estimate location mean locations
#
#########################################################################################
# Non-breeding #

NB_KDE_names<-list.files("data-raw/NonBreeding_Clipped_KDE", pattern="*_clip.txt",full.names=TRUE)

NB_KDE<-lapply(NB_KDE_names,raster)

nGL <- length(NB_KDE_names)

# Get weighted means from KDE #
kdelist<-vector('list',nGL)
NB_kde_long<-NB_kde_lat<-rep(NA,nGL)

for(i in 1:nGL){
kdelist[[i]]<-rasterToPoints(NB_KDE[[i]])
NB_kde_long[i]<-weighted.mean(x=kdelist[[i]][,1],w=kdelist[[i]][,3])
NB_kde_lat[i]<-weighted.mean(x=kdelist[[i]][,2],w=kdelist[[i]][,3])
}

# Replace estimated locations with TRUE capture locations - Jam, Fla, PR birds #
NB_kde_long[c(1,9,21,22,24,25,29,31,36)] <- -77.94 # Jamaica
NB_kde_lat[c(1,9,21,22,24,25,29,31,36)] <- 18.04 # Jamaica

NB_kde_long[c(17,18,23)] <- -80.94 # FLA
NB_kde_lat[c(17,18,23)] <- 25.13 # FLA

NB_kde_long[c(32,33,34,35)] <- -66.86 # PR
NB_kde_lat[c(32,33,34,35)] <- 17.97 # PR

weightedNB <- sf::st_as_sf(data.frame(longitude = NB_kde_long,
                                      latitude = NB_kde_lat),
                           coords = c("longitude","latitude"),
                           crs = 4326)

weightedNBm <- sf::st_transform(weightedNB, 'ESRI:102010')

# USE ONLY BIRDS CAPTURED DURING BREEDING SEASON - GEOLOCATORS #
summerDeploy<-c(2,3,4,5,6,7,8,10,11,12,13,14,15,16,19,20,26,27,28,30)
nB_GL <- length(summerDeploy)

#######################################################################################################################################################
#
# Add the GPS data into the mix
#
#######################################################################################################################################################
GPSdata<-read.csv("data-raw/Ovenbird_GPS_HallworthMT_FirstLast.csv")

nGPS <- nrow(GPSdata)/2

GPSpts <- sf::st_as_sf(GPSdata,
                       coords = c("Longitude","Latitude"),
                       crs = 4326)

GPSptsm <- sf::st_transform(GPSpts, 'ESRI:102010')

# First add GPS locations to both breeding and non-breeding data sets #
cap<-seq(1,2*nGPS,2)
wint<-seq(2,2*nGPS,2)

# Using the weighted locations #

weightedNB_breeDeployOnly <- sf::st_as_sf(data.frame(Longitude = c(NB_kde_long[summerDeploy],GPSdata[wint,2]),
                                                     Latitude = c(NB_kde_lat[summerDeploy],GPSdata[wint,1])),
                                          coords = c("Longitude","Latitude"),
                                          crs = 4326)

NB_breedDeploy<-sf::st_transform(weightedNB_breeDeployOnly, 'ESRI:102010')

isGL <- c(rep(TRUE,20),rep(FALSE,19))
targetPoints <- NB_breedDeploy


###################################################################
#
#  Capture Locations - OriginPoints
#     length = n animals tracked
#
###################################################################


OriginData <- data.frame(Longitude = c(rep(captureLocations[4,1],20),GPSdata[cap,2]),
                         Latitude = c(rep(captureLocations[4,2],20),GPSdata[cap,1]))

Origin <- sf::st_as_sf(OriginData,
                       coords = c("Longitude","Latitude"),
                       crs = 4326)

originPoints<-sf::st_transform(Origin, 'ESRI:102010')

###################################################################
#
#  Origin & Target sites
#
###################################################################
World<-sf::st_read("data-raw/Spatial_Layers/TM_WORLD_BORDERS-0.3.shp")

World<-sf::st_transform(World,'ESRI:102010')

States<-sf::st_read("data-raw/Spatial_Layers/st99_d00.shp")
States<-sf::st_transform(States,'ESRI:102010')

# Non-breeding - Target sites #
Florida<-subset(States,subset=NAME=="Florida")
Florida <- aggregate(Florida, list(Florida$NAME), head, n=1)

Cuba<-subset(World,subset=NAME=="Cuba")

Hisp <- rbind(subset(World,subset=NAME=="Haiti"),
                     subset(World,subset=NAME=="Dominican Republic"))
Hisp <- aggregate(Hisp, list(Hisp$REGION), head, n=1)
Hisp$NAME <- "Hispaniola"

# Change IDs to merge files together

targetSites <- rbind(Cuba[,c("NAME","geometry")],
                     Florida[,c("NAME","geometry")],
                     Hisp[,c("NAME","geometry")])

# Make polygons -
# Breeding - Make square region around capture location - equal size around NH and MD.

# Polygon around MD #
mdvertx<-c((1569680-(536837/2)),(1569680-(536837/2)),(1569680+(536837/2)),(1569680+(536837/2)))
mdverty<-c(-212648,324189,324189,-212648)
mdp<-Polygon(cbind(mdvertx,mdverty))
MDbreedPoly<-SpatialPolygons(list(Polygons(list(mdp),ID=1)))

# Polygon around NH #
nhvertx<-c((1810737-(536837/2)),(1810737-(536837/2)),(1810737+(536837/2)),(1810737+(536837/2)))
nhverty<-c(324189,861026,861026,324189)
nhbp<-Polygon(cbind(nhvertx,nhverty))
NHbreedPoly<-SpatialPolygons(list(Polygons(list(nhbp),ID=1)))

NHbreedPoly<-spChFIDs(NHbreedPoly,"NH")
MDbreedPoly<-spChFIDs(MDbreedPoly,"MD")

crs(NHbreedPoly) <- crs(MDbreedPoly) <- sp::CRS(Lambert)

originSites<-suppressWarnings(spRbind(NHbreedPoly,MDbreedPoly))
crs(originSites)<-Lambert

originSites <- spTransform(originSites,CRS(EquidistConic))

originSites <- sf::st_as_sf(originSites)
###################################################################
#
#  Get relative abundance within breeding "population" polygons #
#
###################################################################

# Breeding Bird Survey Abundance Data #
BBSoven<-raster("data-raw/Spatial_Layers/bbsoven.txt")
crs(BBSoven)<-sf::st_crs(4326)$proj4string

BBSovenMeters<-projectRaster(BBSoven,crs=EquidistConic)

NHbreedPoly <- spTransform(NHbreedPoly,CRS(EquidistConic))
MDbreedPoly <- spTransform(MDbreedPoly,CRS(EquidistConic))

NHabund<-extract(BBSovenMeters,NHbreedPoly)
MDabund<-extract(BBSovenMeters,MDbreedPoly)
TotalOvenAbund<-sum(NHabund[[1]],na.rm=TRUE)+sum(MDabund[[1]],na.rm=TRUE)


BreedRelAbund<-array(NA,c(2,1))
BreedRelAbund[1,1]<-sum(NHabund[[1]],na.rm=TRUE)/TotalOvenAbund
BreedRelAbund[2,1]<-sum(MDabund[[1]],na.rm=TRUE)/TotalOvenAbund

originRelAbund<-BreedRelAbund

###################################################################
#
#  Generate Distance matrices
#
###################################################################
# First need to project from meters to Lat/Long -WGS84
# define current projection #

# project to WGS84
NHbreedPolyWGS<-spTransform(NHbreedPoly,sf::st_crs(4326)$proj4string)
MDbreedPolyWGS<-spTransform(MDbreedPoly,sf::st_crs(4326)$proj4string)

BreedDistMat<-array(NA,c(2,2))
rownames(BreedDistMat)<-colnames(BreedDistMat)<-c(1,2)
diag(BreedDistMat)<-0
BreedDistMat[1,2]<-BreedDistMat[2,1]<-distVincentyEllipsoid(gCentroid(MDbreedPolyWGS, byid=TRUE, id = MDbreedPolyWGS@polygons[[1]]@ID)@coords,
                                                         gCentroid(NHbreedPolyWGS, byid=TRUE, id = NHbreedPolyWGS@polygons[[1]]@ID)@coords)


# Project to WGS84 #
FloridaWGS<-sf::st_transform(Florida,4326)
CubaWGS<-sf::st_transform(Cuba,4326)
HispWGS<-sf::st_transform(Hisp,4326)


NBreedDistMat<-array(NA,c(3,3))
rownames(NBreedDistMat)<-colnames(NBreedDistMat)<-c(3,4,5)

diag(NBreedDistMat)<-0
NBreedDistMat[2,1]<-NBreedDistMat[1,2]<-distVincentyEllipsoid(st_coordinates(st_centroid(FloridaWGS)),
                                                              st_coordinates(st_centroid(CubaWGS)))
NBreedDistMat[3,1]<-NBreedDistMat[1,3]<-distVincentyEllipsoid(st_coordinates(st_centroid(FloridaWGS)),
                                                              st_coordinates(st_centroid(HispWGS)))
NBreedDistMat[3,2]<-NBreedDistMat[2,3]<-distVincentyEllipsoid(st_coordinates(st_centroid(CubaWGS)),
                                                              st_coordinates(st_centroid(HispWGS)))


originDist<-BreedDistMat
targetDist<-NBreedDistMat

###################################################################
#Convert
###################################################################
targetPoints <- sf::st_as_sf(targetPoints)
targetSites <- sf::st_as_sf(targetSites)
originPoints <- sf::st_as_sf(originPoints)
originSites <- sf::st_as_sf(originSites)

###################################################################
#
#  Write required data to the data folder
#
###################################################################

# Put all components of the OVEN Geolocator and GPS data into a named list
OVENdata<-vector('list',12)
names(OVENdata)<-c("geo.bias","geo.vcov","isGL","targetPoints","originPoints",
                   "targetSites","originSites","originRelAbund","originDist",
                   "targetDist","originNames","targetNames")

OVENdata[[1]]<-geo.bias
OVENdata[[2]]<-geo.vcov
OVENdata[[3]]<-isGL
OVENdata[[4]]<-targetPoints
OVENdata[[5]]<-originPoints
OVENdata[[6]]<-targetSites
OVENdata[[7]]<-originSites
OVENdata[[8]]<-originRelAbund
OVENdata[[9]]<-originDist
OVENdata[[10]]<-targetDist
OVENdata[[11]]<-c("NH", "MD")
OVENdata[[12]]<-c("FL", "Cuba", "Hisp")

# Save to data folder
usethis::use_data(OVENdata, overwrite = T)





