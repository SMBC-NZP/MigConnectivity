## ----echo = FALSE, message = FALSE, warning = FALSE------------------------------------------
oo <- getOption("rmarkdown.html_vignette.check_title")
on.exit(options(rmarkdown.html_vignette.check_title = oo))
options(rmarkdown.html_vignette.check_title = FALSE)

library(MigConnectivity)

## --------------------------------------------------------------------------------------------
# Ovenbird data included with the package
data(OVENdata, package = "MigConnectivity") 

names(OVENdata)

## ----message = FALSE, warning = FALSE, error=FALSE, fig.cap="Figure 1. Origin and Target sites used to estimate Ovenbird migratory connectivity using light-level geolocator and GPS tags deployed in the eastern portion of their distribution"----
# Load packages used below
library(sf)
library(terra)
library(maps)

# Set the crs to WGS84
originSitesWGS84 <- st_transform(OVENdata$originSites, 4326)
targetSitesWGS84 <- st_transform(OVENdata$targetSites, 4326)
# Create a simple plot of the origin and and target sites 
op <- graphics::par(no.readonly = TRUE)
on.exit(graphics::par(op))
par(mar=c(0,0,0,0))
plot(originSitesWGS84,
     xlim=c(terra::ext(targetSitesWGS84)[1],
            terra::ext(targetSitesWGS84)[2]),

     ylim=c(terra::ext(targetSitesWGS84)[3],
            terra::ext(originSitesWGS84)[4]))

plot(targetSitesWGS84,
     add = TRUE,
     col=c("gray70","gray35","gray10"))

map("world", add=TRUE)

## ----message=FALSE, warning = FALSE, error=FALSE---------------------------------------------
(OVENpsi <- estTransition(isGL=OVENdata$isGL, # Logical vector: light-level geolocator(T)/GPS(F)
                    isTelemetry = !OVENdata$isGL, # Location vector: light-level geolocator(F)/GPS(T)
                 geoBias = OVENdata$geo.bias, # Light-level geolocator location bias
                 geoVCov = OVENdata$geo.vcov, #Light-level geolocator covariance matrix
                 targetSites = OVENdata$targetSites, # Non-breeding / target sites
                 originSites = OVENdata$originSites, # Breeding / origin sites 
                 targetNames = OVENdata$targetNames, # Names of nonbreeding/target sites
                 originNames = OVENdata$originNames, # Names of breeding/origin sites
                 originPoints = OVENdata$originPoints, # Capture locations 
                 targetPoints = OVENdata$targetPoints, # Target locations from devices
                 resampleProjection = sf::st_crs(OVENdata$targetPoints), 
                 maxTries = 300,
                 verbose = 0,   # output options - see help(estTransition)
                 nSamples = 10)) # This is set low for example 


## ----echo = FALSE----------------------------------------------------------------------------
#saveRDS(OVENpsi,"OVENpsi.rds")

## ----eval = FALSE----------------------------------------------------------------------------
#  # not run
#  str(OVENpsi, max.level = 2)

## ----fig.width=5, fig.height=5,fig.cap="Figure 2. Transition probablities of Ovenbirds from breeding origin sites to non-breeding target sites"----
# THE FIGS ARE TOO SMALL IN VINGETTE HTML SO NEED TO ADD ADDITIONAL PARAMETERS TO MAKE IT LOOK HALFWAY DECENT

plot(OVENpsi, legend = "top", cex = 0.5, las = 2)

## --------------------------------------------------------------------------------------------
# Read in the processed Yellow Warbler data 
# see the Worked Examples Vignette for details on how data structured/created
set.seed(1)
newDir <- tempdir()
baseURL <- 'https://github.com/SMBC-NZP/MigConnectivity/blob/devpsi2/data-raw/YEWA/'
file.name <- "yewa_estTrans_data.rds"
url1 <- paste0(baseURL, file.name, '?raw=true')
temp <- paste(newDir, file.name, sep = '/')
utils::download.file(url1, temp, mode = 'wb')
YEWA_target_sites <- readRDS(temp)
YEWAdata <- readRDS(temp)
unlink(temp)

# take a quick look at the data # 
str(YEWAdata, 1)

## --------------------------------------------------------------------------------------------
## Run analysis for psi
psiYEWA <- estTransition(originSites = YEWAdata$originSites,
                         targetSites = YEWAdata$targetSites,
                         originPoints = YEWAdata$originPoints,
                         targetPoints = YEWAdata$targetPoints,
                         originAssignment = YEWAdata$originAssignment,
                         originNames = YEWAdata$originNames,
                         targetNames = YEWAdata$targetNames,
                         nSamples = 10,
                         isGL = YEWAdata$isGL,
                         isTelemetry = YEWAdata$isTelemetry,
                         isRaster = YEWAdata$isRaster,
                         isProb = YEWAdata$isProb,
                         captured = YEWAdata$captured,
                         geoBias = YEWAdata$geoBias,
                         geoVCov = YEWAdata$geoVCov,
                         originRaster = YEWAdata$originRaster,
                         verbose = 2, 
                         maxTries = 400,
                         resampleProjection = YEWAdata$resampleProjection,
                         nSim = 10,
                         dataOverlapSetting = "none",
                         targetRelAbund = YEWAdata$targetRelAbund,
                         returnAllInput = FALSE)





# Take a look at the results # 
psiYEWA

## ----echo = FALSE----------------------------------------------------------------------------
#saveRDS(psiYEWA,"psiYEWA.rds")

## ----fig.width=7, fig.height=5, fig.cap="Estimated transition probabilities of the Yellow Warbler between breeding and nonbreeding periods"----
par(mar = c(2,3,0,0))
plot(psiYEWA, legend = "top", cex = 0.4)

