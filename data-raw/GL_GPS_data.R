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

###################################################################
#
# geoVcov & geoBias
#
###################################################################

WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Lambert<-"+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"


# Define capture locations in the winter #

captureLocations<-matrix(c(-77.93,18.04,  # Jamaica
                           -80.94,25.13,  # Florida
                           -66.86,17.97,  # Puerto Rico
                           -71.72,43.95), # New Hampshire
                            nrow=4,ncol=2,byrow=TRUE)

CapLocs<-SpatialPoints(captureLocations,CRS(WGS84))

# Project Capture locations #

CapLocsM<-spTransform(CapLocs, CRS(Lambert))

#  Non-breeding GL location data
#  Locations derived during the non-breeding season - ALL birds

# non-breeding file names #

NB_shapefile_names<-list.files("data-raw/GL_NonBreedingFiles", pattern="*.shp",full.names=TRUE)

# Convert into shapefile #

NB_GL<-lapply(NB_shapefile_names, shapefile)

# Better to use variables for frequently used numbers - clearer, more portable code
nGL <- length(NB_GL)

# Remove locations around spring Equinox and potential migration points - same NB time frame as Hallworth et al. 2014 #
# two steps because subset on shapefile doesn't like it in a single step

for(i in 1:nGL){
  NB_GL[[i]]<-subset(NB_GL[[i]],subset=(Month!=3))
}
for(i in 1:nGL){
  NB_GL[[i]]<-subset(NB_GL[[i]],subset=(Month!=4))
}

# Project data from WGS84 into equal area conic to get true error in m/km #

NB_GLm<-vector('list',nGL)
for(i in 1:nGL){
NB_GLm[[i]]<-spTransform(NB_GL[[i]],CRS(Lambert))
}

# Calc known location error in winter #

# Geolocators deployed during non-breeding season - using Raw Location Points #

Jam<-c(30,1,15,16,19,23,25,2,18) # these are the locations within the list of birds captured in Jam
Fla<-c(11,12,17)                 # these are the locations within the list of birds captured in Fla
PR<-c(26,27,28,29)               # these are the locations within the list of birds captured in PR

# Better to use variables for frequently used numbers - clearer, more portable code
nJam <- length(Jam)
nFla <- length(Fla)
nPR <- length(PR)
nNB_GL <- nJam + nFla + nPR

#####################################################################
# Determine the error associated with Geolocators captured at known #
# Non-breeding locations in the Caribbean                           #
#####################################################################

LongError<-rep(NA,nNB_GL) # 16 birds were recovered during the non-breeding season
LatError<-rep(NA,nNB_GL)  # 16 birds were recovered during the non-breeding season

for(i in 1:nJam){   # loop through the 9 Jam birds
#error in
LongError[i]<-mean(NB_GLm[[Jam[i]]]@coords[,1]-CapLocsM@coords[1,1])
LatError[i]<-mean(NB_GLm[[Jam[i]]]@coords[,2]-CapLocsM@coords[1,2])
}
for(i in 1:nFla){   # loop through the 3 Fla birds
LongError[nJam+i]<-mean(NB_GLm[[Fla[i]]]@coords[,1]-CapLocsM@coords[2,1])
LatError[nJam+i]<-mean(NB_GLm[[Fla[i]]]@coords[,2]-CapLocsM@coords[2,2])
}
for(i in 1:nPR){   # loop through the 4 PR birds
LongError[nJam+nFla+i]<-mean(NB_GLm[[PR[i]]]@coords[,1]-CapLocsM@coords[3,1])
LatError[nJam+nFla+i]<-mean(NB_GLm[[PR[i]]]@coords[,2]-CapLocsM@coords[3,2])
}

# Get co-variance matrix for error of known non-breeding deployment sites #

geo.error.model <- lm(cbind(LongError,LatError) ~ 1) # lm does multivariate normal models if you give it a matrix dependent variable!
summary(geo.error.model)

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
NB_KDE_names<-NB_KDE_names[c(1:28,30:37)]
NB_KDE<-lapply(NB_KDE_names,raster)

# Get weighted means from KDE #
kdelist<-vector('list',nGL)
NB_kde_long<-NB_kde_lat<-rep(NA,nGL)

for(i in 1:nGL){
kdelist[[i]]<-rasterToPoints(NB_KDE[[i]])
NB_kde_long[i]<-weighted.mean(x=kdelist[[i]][,1],w=kdelist[[i]][,3])
NB_kde_lat[i]<-weighted.mean(x=kdelist[[i]][,2],w=kdelist[[i]][,3])
}

# Replace estimated locations with TRUE capture locations - Jam, Fla, PR birds #
NB_kde_long[c(1,9,21,22,24,25,29,31,36)]<--77.94 # Jamaica
NB_kde_lat[c(1,9,21,22,24,25,29,31,36)]<-18.04 # Jamaica

NB_kde_long[c(17,18,23)]<--80.94 # FLA
NB_kde_lat[c(17,18,23)]<-25.13 # FLA

NB_kde_long[c(32,33,34,35)]<--66.86 # PR
NB_kde_lat[c(32,33,34,35)]<-17.97 # PR

weightedNB<-SpatialPoints(as.matrix(cbind(NB_kde_long,NB_kde_lat)))
crs(weightedNB)<-WGS84
weightedNBm<-spTransform(weightedNB,CRS(Lambert))

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
GPSpts<-SpatialPoints(as.matrix(cbind(GPSdata[,2],GPSdata[,1]),nrow=nGPS,ncol=2,byrow=TRUE),CRS(WGS84))
GPSptsm<-spTransform(GPSpts,CRS(Lambert))

# First add GPS locations to both breeding and non-breeding data sets #
cap<-seq(1,2*nGPS,2)
wint<-seq(2,2*nGPS,2)

# Using the weighted locations #

weightedNB_breeDeployOnly<-SpatialPoints(as.matrix(cbind(c(NB_kde_long[summerDeploy],GPSdata[wint,2]),c(NB_kde_lat[summerDeploy],GPSdata[wint,1]))))
crs(weightedNB_breeDeployOnly)<-WGS84
NB_breedDeploy<-spTransform(weightedNB_breeDeployOnly,CRS(Lambert))

isGL<-c(rep(TRUE,20),rep(FALSE,19))
targetPoints<-NB_breedDeploy


###################################################################
#
#  Capture Locations - OriginPoints
#     length = n animals tracked
#
###################################################################

Origin<-SpatialPoints(cbind(c(rep(captureLocations[4,1],20),GPSdata[cap,2]),c(rep(captureLocations[4,2],20),GPSdata[cap,1])))
crs(Origin)<-WGS84

originPoints<-spTransform(Origin,CRS(Lambert))


###################################################################
#
#  Origin & Target sites
#
###################################################################
World<-shapefile("data-raw/Spatial_Layers/TM_WORLD_BORDERS-0.3.shp")
World<-spTransform(World,CRS(Lambert))
States<-shapefile("data-raw/Spatial_Layers/st99_d00.shp")
States<-spTransform(States,CRS(Lambert))

# Non-breeding - Target sites #
Florida<-subset(States,subset=NAME=="Florida")
Florida<-gUnaryUnion(Florida)
Cuba<-subset(World,subset=NAME=="Cuba")
Hisp<-gUnion(subset(World,subset=NAME=="Haiti"),subset(World,subset=NAME=="Dominican Republic"))

# Change IDs to merge files together
Cuba<-spChFIDs(Cuba,"Cuba")
Florida<-spChFIDs(Florida,"Florida")
Hisp<-spChFIDs(Hisp,"Hisp")

#Combine into a single SpatialPolygon
WinterRegion1 <- spRbind(Florida,Cuba)
WinterRegions<-spRbind(WinterRegion1,Hisp)

targetSites<-WinterRegions

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

originSites<-spRbind(NHbreedPoly,MDbreedPoly)
crs(originSites)<-Lambert

###################################################################
#
#  Get relative abundance within breeding "population" polygons #
#
###################################################################

# Breeding Bird Survey Abundance Data #
BBSoven<-raster("data-raw/Spatial_Layers/bbsoven.txt")
crs(BBSoven)<-WGS84
BBSovenMeters<-projectRaster(BBSoven,crs=Lambert)


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

crs(NHbreedPoly)<-crs(MDbreedPoly)<-Lambert

# project to WGS84
NHbreedPolyWGS<-spTransform(NHbreedPoly,CRS(WGS84))
MDbreedPolyWGS<-spTransform(MDbreedPoly,CRS(WGS84))

BreedDistMat<-array(NA,c(2,2))
rownames(BreedDistMat)<-colnames(BreedDistMat)<-c(1,2)
diag(BreedDistMat)<-0
BreedDistMat[1,2]<-BreedDistMat[2,1]<-distVincentyEllipsoid(gCentroid(MDbreedPolyWGS, byid=TRUE, id = MDbreedPolyWGS@polygons[[1]]@ID)@coords,
                                                         gCentroid(NHbreedPolyWGS, byid=TRUE, id = NHbreedPolyWGS@polygons[[1]]@ID)@coords)


# Project to WGS84 #
FloridaWGS<-spTransform(Florida,CRS(WGS84))
CubaWGS<-spTransform(Cuba,CRS(WGS84))
HispWGS<-spTransform(Hisp,CRS(WGS84))


NBreedDistMat<-array(NA,c(3,3))
rownames(NBreedDistMat)<-colnames(NBreedDistMat)<-c(3,4,5)

diag(NBreedDistMat)<-0
NBreedDistMat[2,1]<-NBreedDistMat[1,2]<-distVincentyEllipsoid(gCentroid(FloridaWGS, byid=FALSE)@coords,
                                                           gCentroid(CubaWGS, byid=FALSE)@coords)
NBreedDistMat[3,1]<-NBreedDistMat[1,3]<-distVincentyEllipsoid(gCentroid(FloridaWGS, byid=FALSE)@coords,
                                                           gCentroid(HispWGS, byid=FALSE)@coords)
NBreedDistMat[3,2]<-NBreedDistMat[2,3]<-distVincentyEllipsoid(gCentroid(CubaWGS, byid=FALSE)@coords,
                                                         gCentroid(HispWGS, byid=FALSE)@coords)


originDist<-BreedDistMat
targetDist<-NBreedDistMat

###################################################################
#
#  Write required data to the data folder
#
###################################################################

# Put all components of the OVEN Geolocator and GPS data into a named list
OVENdata<-vector('list',10)
names(OVENdata)<-c("geo.bias","geo.vcov","isGL","targetPoints","originPoints",
                   "targetSites","originSites","originRelAbund","originDist","targetDist")

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

# Save to data folder
devtools::use_data(OVENdata, overwrite = T)





