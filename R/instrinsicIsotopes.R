#' Generate probabilitistic isotope assignments
#'
#' The \code{isoAssign} function generates origin assignments using stable-hydrogen isotopes in tissue. The function generates
#' a probability surface of origin assignment from a vector of stable-isotope values for each animal/sample of interest.
#' Probabilistic assignments are constructed by first converting observed stable-isotope ratios (isoscape) in either precipitation or surface
#' waters into a 'tissuescape' using a user-provided intercept, slope and standard deviation.
#' See \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0035137}{Hobson et. al. (2012)}.
#'
#'
#' @param isovalues vector of tissue isotope values
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param odds odds ratio to use to set likely and unlikely locations defaults to 0.67
#' @param restrict2Likely if \code{TRUE} restricts locations to fall within the 'likely' assignment
#'        locations.
#' @param nSamples integer specifying how many random samples to draw from a multinomial distribution.
#' @param sppShapefile SpatialPolygon layer defining species range. Assignments are restricted to these
#'        areas.
#' @param relAbund raster with relative abundance (must match extent of isotope assignment)
#' @param isoWeight weighting value to apply to isotope assignment
#' @param abundWeight weighting value to apply to relative abundance prior
#' @param population vector identifying location where animal was captured. Same order as \code{isovalues}
#' @param assignExtent definition for the extent of the assignment. Can be used in place of \code{sppShapefile} to
#'        limit assignment. Input should follow \code{c(xmin,xmax,ymin,ymax)} in degrees longitude and latitude.
#' @param element The elemental isotope of interest. Currently the only
#'     elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface if "TRUE" returns surface water values. Defaults is 'FALSE'
#'     which returns the isotopes ratio found in precipitation.
#' @param period The time period of interest. If 'Annual' returns a raster
#'     of mean annual values in precipitation for the \code{element}. If
#'     'GrowingSeason' returns growing season values in precipitation for
#'      \code{element} of interest.
#' @param seed numeric value fed to \code{set.seed} for random number generation. Default = NULL.
#' @param verbose takes values 0, 1 (default) or 2. 0 prints no output during run. 1 prints
#'  a message detailing where in the process the function is. 2 prints the animal currently being sampled.
#' @return returns an \code{isoAssign} object containing the following:
#'  \describe{
#'   \item{\code{probassign}}{raster stack of individual probabilistic assignments}
#'   \item{\code{oddsassign}}{raster stack that includes likely vs unlikely origin for each animal}
#'   \item{\code{popassign}}{a raster for population level assignment (sum of \code{oodsassign} if \code{population} = NULL).
#'   If \code{population} is a vector then returns a raster stack for each unique \code{population} provided}
#'   \item{\code{probDF}}{data.frame of individual probability surfaces}
#'   \item{\code{oddsDF}}{data.frame of likely vs unlikley surfaces}
#'   \item{\code{popDF}}{data.frame of population level assignment}
#'   \item{\code{SingeCell}}{array of coordinates (longitude,latitude) for single cell assignment}
#'   \item{\code{targetSites}}{\code{SpatialPolygons} layer representing isotope bands equivilant to \code{isoSTD}}
#'   \item{\code{RandomSeed}}{the RNG seed used when generating locations from the multinomial distribution}
#'   }
#'
#' @seealso{\code{\link{weightAssign}}}
#' @export
#'
#' @examples
#' \dontrun{
#' OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
#' OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#' raster::crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#'
#' OVENvals <- read.csv("data-raw/deltaDvalues.csv")
#'
#'a <- Sys.time()
#'b <- isoAssign(isovalues = OVENvals[,2],
#'               isoSTD = 12,
#'               intercept = -10,
#'               slope = 0.8,
#'               odds = NULL,
#'               restrict2Likely = TRUE,
#'               nSamples = 1000,
#'               sppShapefile = OVENdist,
#'               assignExtent = c(-179,-60,15,89),
#'               element = "Hydrogen",
#'               surface = FALSE,
#'               period = "Annual")
#'Sys.time()-a
#'}
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of
#' Mexico. Ecography 42: 658-669.
#'
#' Hobson, K. A., S. L. Van Wilgenburg, L. I. Wassenaar, and K. Larson. 2012.
#' Linking hydrogen isotopes in feathers and precipitation: sources of
#' variance and consequences for assignment to isoscapes. PLoS ONE 7: e35137.

isoAssign <- function(isovalues,
                      isoSTD,
                      intercept,
                      slope,
                      odds = 0.67,
                      restrict2Likely = TRUE,
                      nSamples = NULL,
                      sppShapefile = NULL,
                      relAbund = NULL,
                      isoWeight = NULL,
                      abundWeight = NULL,
                      population = NULL,
                      assignExtent = c(-179,-60,15,89),
                      element = "Hydrogen",
                      surface = FALSE,
                      period = "Annual",
                      seed = NULL,
                      verbose=1){
# force verbose to default when outside specified range.
if(!(verbose %in% c(0,1,2))){verbose = 1}

# download isoscape map
isomap <- getIsoMap(element = element, surface = surface, period = period)

# 1. if sppShapefile == NULL - use extent option
if(is.null(sppShapefile)){
  isomap <- raster::crop(isomap,raster::extent(assignExtent))
}

# Series of checks for a species range map inputs
if(!is.null(sppShapefile)){
# 2. if sppShapefile provided check that it has a projection defined
#    if not stop - if so, mask the isoscape to range
  if(is.na(raster::crs(sppShapefile))){
    stop("coordinate system needed for sppShapefile")}
# 3. if the projections don't match - project into same as isomap then mask
  # quick check
  if(class(sppShapefile) %in% c("SpatialPolygons","SpatialPolygonsDataFrame")){
  if(sppShapefile@proj4string@projargs != isomap@crs@projargs){
    sppShapefile <- sp::spTransform(sppShapefile, sp::CRS(isomap@crs@projargs))
  }
  }
  if(class(sppShapefile)[1] %in% "sf"){
    if(!identical(sf::st_crs(sppShapefile),sf::st_crs(4326))){
      sppShapefile <- sf::st_transform(sppShapefile, 4326)
    }
    }


# convert sp file to sf
if(class(sppShapefile) %in% c("SpatialPolygon","SpatialPolygonDataFrame")){
  sppShapefile <- sf::st_as_sf(sppShapefile)
}

if(verbose>0){cat("\n Restricting possible assignments to species distribution \n")}

sppShapefile$INOUT<-1

# mask the isomap to sppShapefile
isomap <- raster::crop(isomap, sppShapefile)
isomap <- raster::mask(isomap, sppShapefile)
}
if(!is.null(relAbund) && !inherits(relAbund,"RasterLayer")){stop("relAbund should be a raster layer")}
if(!is.null(relAbund) && inherits(relAbund,"RasterLayer")){
# if isomap and relAbund don't have the same resolution and/or extent
# change to relAbund to match isomap
if(!raster::compareRaster(isomap,relAbund, values = FALSE)){
# project to match isomap
relAbund <- raster::projectRaster(relAbund,isomap)
# re-scale to ensure sums to 1
relAbund <- relAbund/raster::cellStats(relAbund,sum)
}
# if relAbund isn't a probability surface - generate probability surface #
if(raster::cellStats(relAbund,sum)!=1){relAbund <- relAbund/raster::cellStats(relAbund,sum)}
}
# generate a 'feather'/animal isoscape
animap <- raster::calc(isomap, function(x){y <- slope*x+intercept})

# generate targetSites - seq from min to max values by isoSTD
isocut <- raster::cut(animap, breaks= seq(from = raster::cellStats(animap,min),
                                  to = raster::cellStats(animap,max),
                                  by = isoSTD))

# use those cuts to make polygons
targetSites <- raster::rasterToPolygons(isocut, dissolve = TRUE)
raster::crs(targetSites) <- raster::crs(isomap)

#if sppShapefile !NULL then clip targetSites to distribution
if(!is.null(sppShapefile)){
targetSites <- raster::intersect(targetSites,sppShapefile)
}
# rename the targetSites to simplify output
targetSites<-targetSites[,1]
names(targetSites) <- c("targetSite")
#targetSites <- rgeos::gUnaryUnion(targetSites, id=targetSites$targetSite)
# spatially explicit assignment
assign <- function(x,y) {((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))}

# apply the assignment function to all input values
if(verbose>0){cat("\n Generating probabilistic assignments \n")}

assignments <- lapply(isovalues, FUN = function(x){assign(x, y = animap)})

# stack the assignment probabilities into a single raster stack
assignments <- raster::stack(assignments)

# Transform the assignments into a true probability surface #
assign2prob <- assignments/raster::cellStats(assignments, sum)

# Weighted Assignments
if(inherits(relAbund,"RasterLayer") && is.null(isoWeight) && is.null(abundWeight)){
  if(verbose>0){cat("\n Creating posterior assignments where isotope & abundance have equal weight \n")}
assign2prob <- assign2prob*relAbund
assign2prob <- assign2prob/raster::cellStats(assign2prob,sum)
}
if(inherits(relAbund,"RasterLayer") && !is.null(isoWeight) && !is.null(abundWeight)){
  if(verbose>0){cat("\n Creating weighted posterior assignments \n")}
  isoWeight <- 10^isoWeight
  abundWeight <- 10^abundWeight
assign2prob <- (assign2prob^isoWeight)*(relAbund^abundWeight)
assign2prob <- assign2prob/raster::cellStats(assign2prob,sum)
}

if(inherits(relAbund,"RasterLayer") && is.null(isoWeight) && !is.null(abundWeight)){
  if(verbose>0){cat("\n Creating posterior abundance weighted assignments \n")}
  abundWeight <- 10^abundWeight
assign2prob <- assign2prob*(relAbund^abundWeight)
assign2prob <- assign2prob/raster::cellStats(assign2prob,sum)
}
if(inherits(relAbund,"RasterLayer") && !is.null(isoWeight) && is.null(abundWeight)){
  if(verbose>0){cat("\n Creating posterior isotope weighted assignments \n")}
  isoWeight <- 10^isoWeight
assign2prob <- (assign2prob^isoWeight)*relAbund
assign2prob <- assign2prob/raster::cellStats(assign2prob,sum)
}
# Create a dataframe with XY coords and probabilites for each animal
assign2probDF <- data.frame(raster::rasterToPoints(assign2prob))

# function to make an odds ratio (likely vs unlikely) assignment
oddsFun <- function(x,odds = odds){
  predict(smooth.spline(x = cumsum(sort(x)),
                        sort(x),
                        spar = 0.1),(1-odds))$y
}

# if odds is left null - use default of 0.33
if (is.null(odds)){odds <- 0.67}

# extract values from the probability assignment
matvals <- raster::rasterToPoints(assign2prob)
# XY coords of raster
matvalsXY <- matvals[,1:2]

# drop XY from matvals
matvals <- matvals[,-(1:2)]

if(verbose>0){cat("\n Generating likely vs unlikely assignments \n")}
# apply the odds function

cuts <- apply(matvals,2,FUN = oddsFun,odds = odds)

# reclassify the rasters based on likely v unlikely

step1 <- mapply(FUN = function(x,y){raster::reclassify(assign2prob[[x]],cbind(0,y,0))},
                x = 1:raster::nlayers(assign2prob),
                y = cuts)
step1 <- raster::stack(step1)
step2 <- mapply(FUN = function(x,y){raster::reclassify(step1[[x]],cbind(y,1,1))},
                x = 1:raster::nlayers(step1),
                y = cuts)
step2 <- raster::stack(step2)

#step1 <- raster::reclassify(assign2prob,cbind(0,cuts,0))
#step2 <- raster::reclassify(step1,cbind(cuts,1,1))

# convert to dataframe
LikelyUnlikely <- raster::rasterToPoints(step2)
step2DF <- data.frame(LikelyUnlikely)

if(is.null(population)){
# Return the population level odds assignment - i.e, how many animals
SamplePop <- sum(step2)
# convert to dataframe
SamplePopDF <- data.frame(raster::rasterToPoints(SamplePop))
} else {
nPop <- length(unique(population))
Pops <- unique(population)
POPs <- POPSdf <- vector("list",nPop)
for(p in 1:nPop){
aniPop <- which(population == Pops[p])
pop_p <- step2[[aniPop]]
POPs[[p]] <- sum(pop_p)
POPSdf[[p]] <- data.frame(raster::rasterToPoints(POPs[[p]]))
}
SamplePop <- raster::stack(POPs)
names(SamplePop)<-Pops
SamplePopDF <- do.call('cbind',POPSdf)
}

# SINGLE CELL PROBABILITY ASSIGNMENTS - this makes MC possible with isotopes
if (is.null(nSamples)){ nSamples <- 1000}
if(is.null(seed)){seed <- as.numeric(Sys.time())}

# Set random number generator to replicate results
set.seed(seed)

# generate empty array to fill with locations
# make a simulated array twice the size to weed out locations
# that fall outside of distribution
xysimulation <- array(NA,c(nSamples+floor((nSamples/2)),2,raster::nlayers(assign2prob)))
# give names for sf to convert down the line
dimnames(xysimulation)[[2]] <- c("Longitude","Latitude")

xysim <- array(NA, c(nSamples, 2, raster::nlayers(assign2prob)))
# name the array
  #dimnames(xysim)[[1]] <- 1:nSamples
  dimnames(xysim)[[2]] <- c("Longitude","Latitude")
  dimnames(xysim)[[3]] <- names(assign2prob)

# converts raster to matrix of XY then probs
#matvals <- raster::rasterToPoints(assign2prob)

# Restrict random point estimates to 'likely' origin area #
if(restrict2Likely){
    matvals<- matvals* LikelyUnlikely[,3:ncol(LikelyUnlikely)]
}

  if(verbose>0){cat("\n Generating single cell assignments \n")}
  # This draws samples nSamples per animal (faster than looping over nSamples) and fills the xysim with x,y coords
 if(verbose>1){
    cat("\b\b\b\b\b\b");
    cat("\n      animal # ",sprintf("%3d",1),"\n");
    utils::flush.console()
  }
for(i in 1:ncol(matvals)) {
  if(verbose>1){
    cat("\b\b\b\b\b\b");
    cat(sprintf("%3d",i),"\n");
    utils::flush.console()
  }

    multidraw <- rmultinom(n = nSamples+floor((nSamples/2)), size = 1, prob = matvals[,i])
    xysimulation[,1,i] <- matvalsXY[which(multidraw == 1, arr.ind = TRUE)[,1],1]
    xysimulation[,2,i] <- matvalsXY[which(multidraw == 1, arr.ind = TRUE)[,1],2]
    # check to see which are in the distribution and which fall outside
    if(!is.null(sppShapefile)){
      sppShapefile <- sf::st_as_sf(sppShapefile)
   # randpoints <- sp::SpatialPoints(cbind(xysimulation[,1,i],xysimulation[,2,i]),
   #                                 sp::CRS(sppShapefile@proj4string@projargs))
     randpoints <- sf::st_as_sf(data.frame(xysimulation[,,i]),
                                coords = c("Longitude","Latitude"),
                                crs = 4326)

   # inout <- sp::over(randpoints,sppShapefile)
    inout <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = randpoints,
                                                                   y = sppShapefile,
                                                          sparse = TRUE))))
    # How many are in
   # InDist <- randpoints[which(inout$INOUT == 1),]
    InDist <- randpoints[which(inout > 0),]
    samplecoords <- sample(1:nrow(InDist),
                           size = nSamples,
                           replace = FALSE)
    #xysim[,1,i] <- InDist@coords[samplecoords,1]
    #xysim[,2,i] <- InDist@coords[samplecoords,2]
    xysim[,1,i] <- sf::st_coordinates(InDist)[samplecoords,1]
    xysim[,2,i] <- sf::st_coordinates(InDist)[samplecoords,2]
    }else{
    randsamples <- sample(1:nrow(xysimulation),size = nSamples,replace = FALSE)
    xysim[,1,i] <- xysimulation[randsamples,1,i]
    xysim[,2,i] <- xysimulation[randsamples,2,i]
    }
    # while numInDist is less than nSamples - redraw and fill NA with multinomial draw
  #samplenum <- 1
  #while(numInDist<nSamples){
   # if(verbose){cat(samplenum," resampling events to ensure points fall within distribution \n")}
   # identify which rows need new points
   # NeedsNewPoints <- which(is.na(inout$INOUT),arr.ind = TRUE)
   # how many are needed
   # NumNewDraws <- length(NeedsNewPoints)
   # sample 1000 new draws
   # NewDraws <- rmultinom(n = nSamples, size = 1, prob = matvals[,i])
   # take a subsample of the new draws of size NumNewDraws
   # subsample <- sample(x = which(NewDraws == 1,arr.ind = TRUE)[,1],
   #                     size = NumNewDraws,
   #                     replace = FALSE)

   # fill with the subsample
   # xysim[NeedsNewPoints,1,i-2] <- matvals[subsample,1]
   # xysim[NeedsNewPoints,2,i-2] <- matvals[subsample,2]

   # check to see which are in the distribution and which fall outside
   # inout <- sp::over(sp::SpatialPoints(cbind(xysim[,1,i-2],xysim[,2,i-2]),
   #                                     sp::CRS(sppShapefile@proj4string@projargs)),sppShapefile)
   # How many are in
   # numInDist <- sum(inout$INOUT,na.rm = TRUE)
   # while numInDist is less than nSamples - redraw and fill NA with multinomial draw
   # samplenum <- samplenum+1
  #}
}


isoAssignReturn <- structure(list(probassign = assign2prob,
                        oddsassign = step2,
                        popassign = SamplePop,
                        probDF = assign2probDF,
                        oddsDF = step2DF,
                        popDF = SamplePopDF,
                        SingleCell = xysim,
                        targetSites = targetSites,
                        RandomSeed = seed),
                        class = c("isoAssign", "intrinsicAssign"))


if(verbose>0){cat("\n Done \n"); cat("\n Random number seed used = ", seed,"\n")}
return(isoAssignReturn)
}



#' Get Isoscape map
#' getIsoMap
#'
#' The \code{getIsoMap} function downloads predicted isoscape maps from
#'  \url{http://wateriso.utah.edu/waterisotopes/}. The function first checks
#' whether the isoscapes are located within the current working directory
#' \code{getwd()}. If a local copy of the isoscape is found, it's read into
#' the environment. If not, the isoscape is downloaded and imported
#' as a raster.
#'
#'
#' @param element The elemental isotope of interest. Currently the only
#'     elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface if "TRUE" returns surface water values. Default is 'FALSE'
#'     which returns the isotopes ratio found in precipitation.
#' @param period The time period of interest. If 'Annual' (default) returns a raster
#'     of mean annual values in precipitation for the \code{element}. If
#'     'GrowingSeason' returns growing season values in precipitation for
#'      \code{element} of interest.
#'
#' @return returns a global \code{RasterLayer} (resolution = 0.333'x0.3333') object for the \code{element} and \code{period} of interest
#'
#' @export
#' @examples
#' \dontrun{
#' map <- getIsoMap(element = "Hydrogen", period = "Annual")
#' }

getIsoMap<-function(element = "Hydrogen", surface = FALSE, period = "Annual"){

  # and read into R as raster - otherwise read MAD into R
  if(!(element %in% c("Hydrogen","Oxygen")))
    stop("element must be either Hydrogen or Oxygen")

  if(!(period %in% c("Annual","GrowingSeason")))
    stop("period must be either Annual or GrowingSeason")

  # if "/AnnualD" isn't in working directory somewhere download MAD from website
  # Download Mean Annual Deuterium values from wateriso.utah.edu #
  haveIsoMap <- list.dirs(path = getwd(), recursive = TRUE)
  if(surface == TRUE){
    if(element == "Hydrogen"){
      if(!(paste0(getwd(),"/mwswh_fin") %in% haveIsoMap)){

        # Create temporary file #
        tf <- tempfile(pattern = "file", tmpdir = getwd(), fileext = "")

        # Download the file #
        utils::download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/Surface_H.zip",
                      destfile = tf,
                      quiet = TRUE,
                      extra = getOption("download.file.extra"))

        # unzip the downloaded file #
        utils::unzip(zipfile = tf,
              files = NULL, list = FALSE, overwrite = TRUE,
              junkpaths = FALSE, exdir = getwd(), unzip = "internal",
              setTimes = FALSE)

        # Delete zipped folder #
        file.remove(tf)

        m_s_d <- raster::raster(paste0(getwd(),"/mwswh_fin/w001001.adf"))
      }else{
        m_s_d <- raster::raster(paste0(haveIsoMap[grep(haveIsoMap,pattern = "/mwswh_fin$")],"/w001001.adf"))
      }
      names(m_s_d)<-"MeanSurfaceD"
      return(m_s_d)
    }
    if(element == "Oxygen"){
      if(!(paste0(getwd(),"/mwswh_fin") %in% haveIsoMap)){

        # Create temporary file #
        tf <- tempfile(pattern = "file", tmpdir = getwd(), fileext = "")

        # Download the file #
        utils::download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/Surface_O.zip",
                      destfile = tf,
                      quiet = TRUE,
                      extra = getOption("download.file.extra"))

        # unzip the downloaded file #
        utils::unzip(zipfile = tf,
              files = NULL, list = FALSE, overwrite = TRUE,
              junkpaths = FALSE, exdir = getwd(), unzip = "internal",
              setTimes = FALSE)

        # Delete zipped folder #
        file.remove(tf)

        m_s_o <- raster::raster(paste0(getwd(),"/mwswo_fin/w001001.adf"))
      }else{
        m_s_o <- raster::raster(paste0(haveIsoMap[grep(haveIsoMap,pattern = "/mwswh_fin$")],"/w001001.adf"))
      }
      names(m_s_o)<-"MeanSurfaceO"
      return(m_s_o)
    }

  }
  if(element == "Hydrogen" & period == "Annual"){
    if(!(paste0(getwd(),"/AnnualD") %in% haveIsoMap)){

      # Create temporary file #
      tf <- tempfile(pattern = "file", tmpdir = getwd(), fileext = "")

      # Download the file #
      utils::download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/AnnualD.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
            files = NULL, list = FALSE, overwrite = TRUE,
            junkpaths = FALSE, exdir = getwd(), unzip = "internal",
            setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      m_a_d <- raster::raster(paste0(getwd(),"/AnnualD/mad/w001001.adf"))
    }else{
      m_a_d <- raster::raster(paste0(haveIsoMap[grep(haveIsoMap,pattern = "/AnnualD/mad$")],"/w001001.adf"))
    }
    names(m_a_d)<-"MeanAnnualD"
    return(m_a_d)
  }

  if(element == "Hydrogen" & period == "GrowingSeason"){
    if(!(paste0(getwd(),"/GSD") %in% haveIsoMap)){

      # Create temporary file #
      tf <- tempfile(pattern = "file", tmpdir = getwd(), fileext = "")

      # Download the file #
      utils::download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/GSD.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
            files = NULL, list = FALSE, overwrite = TRUE,
            junkpaths = FALSE, exdir = getwd(), unzip = "internal",
            setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      g_s_d <- raster::raster(paste0(getwd(),"/GSD/gsd/w001001.adf"))
    }else{
      g_s_d <- raster::raster(paste0(haveIsoMap[grep(haveIsoMap,pattern = "/GSD/gsd$")],"/w001001.adf"))
    }
    g_s_d[g_s_d == -999]<-NA
    names(g_s_d)<-"GrowingSeasonD"
    return(g_s_d)
  }
  if(element == "Oxygen" & period == "Annual"){
    if(!(paste0(getwd(),"/GSO") %in% haveIsoMap)){

      # Create temporary file #
      tf <- tempfile(pattern = "file", tmpdir = getwd(), fileext = "")

      # Download the file #
      utils::download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/AnnualO.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
            files = NULL, list = FALSE, overwrite = TRUE,
            junkpaths = FALSE, exdir = getwd(), unzip = "internal",
            setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      m_a_o <- raster::raster(paste0(getwd(),"/AnnualO/mao/w001001.adf"))
    }else{
      m_a_o <- raster::raster(paste0(haveIsoMap[grep(haveIsoMap,pattern = "/AnnualO/mao$")],"/w001001.adf"))
    }
    names(m_a_o)<-"MeanAnnualO"
    return(m_a_o)
  }

  if(element == "Oxygen" & period == "GrowingSeason"){
    if(!(paste0(getwd(),"/GSO") %in% haveIsoMap)){

      # Create temporary file #
      tf <- tempfile(pattern = "file", tmpdir = getwd(), fileext = "")

      # Download the file #
      utils::download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/GSO.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
            files = NULL, list = FALSE, overwrite = TRUE,
            junkpaths = FALSE, exdir = getwd(), unzip = "internal",
            setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      g_s_o <- raster::raster(paste0(getwd(),"/GSO/gso/w001001.adf"))
    }else{
      g_s_o <- raster::raster(paste0(haveIsoMap[grep(haveIsoMap,pattern = "/GSO/gso$")],"/w001001.adf"))
    }
    g_s_o[g_s_o == -999]<-NA
    names(g_s_o)<-"GrowingSeasonO"
    return(g_s_o)
  }

}

#' Calculate Weights for Isotope Assignments
#' weightAssign
#'
#' The primary purpose of this function is to determine whether weighting likelihood based isotope assignments
#' and prior information, such as relative abundance can improve the model performance compared to the
#' isotope-only model. To do this, we raise the likelihood and prior values to powers from 0.1
#' to 10 and measure model performance using the assignment error rate and assignment area. Weights < 1 flatten
#' the likelihood/prior distributions (giving relatively more weight to smaller values) and weights > 1
#' sharpen the distributions (giving relatively less weight to smaller values. The \code{weightAssign} function
#' generates origin assignments using stable-hydrogen isotopes in tissue. If first generates
#' a probability surface of origin assignment from a vector of stable-isotope values for each animal/sample
#' captured at a known location. Probabilistic assignments are constructed by first converting observed
#' stable-isotope ratios (isoscape) in either precipitation or surface waters into a 'tissuescape' using
#' a user-provided intercept, slope and standard deviation. See
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0035137}{Hobson et. al. (2012)}.
#'
#'  See \href{https://onlinelibrary.wiley.com/doi/10.1002/ece3.2605}{Rushing et al. (2017)} for more information.
#'
#' @param knownLocs matrix of capture locations of the same length as \code{isovalues}
#' @param isovalues vector of tissue isotope values from known locations
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param odds odds ratio to use to set likely and unlikely locations defaults to 0.67
#' @param relAbund raster layer of relative abundance that sums to 1.
#' @param weightRange vector of length 2 within minimum and maximum values to weight isotope and relative abundance.
#'        Default = c(-1,1)
#' @param sppShapefile SpatialPolygon layer defining species range. Assignments are restricted to these
#'        areas.
#' @param assignExtent definition for the extent of the assignment. Can be used in place of \code{sppShapefile} to
#'        limit assignment. Input should follow \code{c(xmin,xmax,ymin,ymax)} in degrees longitude and latitude.
#' @param element The elemental isotope of interest. Currently the only
#'     elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface if "TRUE" returns surface water values. Defaults is 'FALSE'
#'     which returns the isotopes ratio found in precipitation.
#' @param period The time period of interest. If 'Annual' returns a raster
#'     of mean annual values in precipitation for the \code{element}. If
#'     'GrowingSeason' returns growing season values in precipitation for
#'      \code{element} of interest.
#'
#' @return returns an \code{weightAssign} object containing the following:
#'   \describe{
#'    \item{\code{top}}{data.frame with the optimal weightings}
#'    \item{\code{frontier}}{data.frame with values that fall along the Pareto frontier}
#'    \item{\code{performance}}{data.frame with error rate and assignment area for each weight combination}
#' }
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of migratory
#' connectivity for birds en route to breeding through the Gulf of Mexico.
#' Ecography 42: 658-669.
#'
#' Rushing, C. S., P. P. Marra and C. E. Studds. 2017. Incorporating breeding
#' abundance into spatial assignments on continuous surfaces. Ecology and
#' Evolution 3: 3847-3855.
#' @export
#'
#' @examples
#' \dontrun{
#' OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
#' OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#' raster::crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#'
#' OVENvals <- read.csv("data-raw/deltaDvalues.csv")
#'
#' HBEFbirds <- OVENvals[grep("NH",OVENvals[,1]),]
#' knownLocs <- cbind(rep(-73,nrow(HBEFbirds)),rep(43,nrow(HBEFbirds)))
#'
#' utils::download.file("https://www.mbr-pwrc.usgs.gov/bbs/ra15/ra06740.zip", destfile = "oven.zip")
#' utils::unzip("oven.zip")
#' oven_dist <- raster::shapefile("ra06740.shp")
#'
#' # Empty raster with the same dimensions as isoscape and Ovenbird distribution
#' r <- raster::raster(nrow = 83, ncol = 217, res = c(0.333333, 0.333333),
#'                     xmn = -125.0001, xmx = -52.66679, ymn = 33.33321, ymx = 60.99985,
#'                     crs = MigConnectivity::projections$WGS84)
#'
#' relativeAbund <- raster::rasterize(sp::spTransform(oven_dist, sp::CRS(r@crs@projargs)),r)
#' relativeAbund <- relativeAbund /raster::cellStats(relativeAbund ,sum)
#'
#' BE <- weightAssign(knownLocs = knownLocs,
#'                  isovalues = HBEFbirds[,2],
#'                  isoSTD = 12,
#'                  intercept = -10,
#'                  slope = 0.8,
#'                  odds = 0.67,
#'                  relAbund = relativeAbund,
#'                  weightRange = c(-1,1),
#'                  sppShapefile = OVENdist,
#'                  assignExtent = c(-179,-60,15,89),
#'                  element = "Hydrogen",
#'                  surface = FALSE,
#'                  period = "Annual")
#'}
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. In revision. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of Mexico.
#'
#' Hobson, K. A., S. L. Van Wilgenburg, L. I. Wassenaar, and K. Larson. 2012.
#' Linking hydrogen isotopes in feathers and precipitation: sources of
#' variance and consequences for assignment to isoscapes. PLoS ONE 7: e35137.
#'
#' Rushing, C. S., P. P. Marra, and C. E. Studds. 2017. Incorporating breeding
#' abundance into spatial assignments on continuous surfaces. Ecology and
#' Evolution 7: 3847-3855.
weightAssign <- function(knownLocs,
                         isovalues,
                         isoSTD,
                         intercept,
                         slope,
                         odds = 0.67,
                         relAbund,
                         weightRange = c(-1,1),
                         sppShapefile = NULL,
                         assignExtent = c(-179,-60,15,89),
                         element = "Hydrogen",
                         surface = FALSE,
                         period = "Annual"){
a <- Sys.time()
  # Check to make sure knownLocs are added by user
  if(is.null(knownLocs)){
    stop("Known locations are needed to determine weighted assignemnts")}
  if(nrow(knownLocs)!=length(isovalues)){
    stop("A known location (knownLocs) is needed for each isotope value (isovalues)")
  }
  # download isoscape map
  isomap <- getIsoMap(element = element, surface = surface, period = period)

  # 1. if sppShapefile == NULL - use extent option
  if(is.null(sppShapefile)){
    isomap <- raster::crop(isomap,raster::extent(assignExtent))
  }

  # Series of checks for a species range map inputs
  if(!is.null(sppShapefile)){
    # 2. if sppShapefile provided check that it has a projection defined
    #    if not stop - if so, mask the isoscape to range
    if(is.na(raster::crs(sppShapefile))){
      stop("coordinate system needed for sppShapefile")}
    # 3. if the projections don't match - project into same as isomap then mask
    if(sppShapefile@proj4string@projargs != isomap@crs@projargs){
      sppShapefile <- sp::spTransform(sppShapefile, sp::CRS(isomap@crs@projargs))
    }

    cat("\n Restricting possible assignments to species distribution \n")

    sppShapefile$INOUT<-1

    # mask the isomap to sppShapefile
    isomap <- raster::mask(isomap, sppShapefile)
    isomap <- raster::crop(isomap, sppShapefile)
  }

  # generate a 'feather'/animal isoscape
  animap <- raster::calc(isomap, function(x){y <- slope*x+intercept})

  # spatially explicit assignment
  assign <- function(x,y) {((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))}

  # apply the assignment function to all input values
  cat("\n Generating probabilistic assignments \n")

  assignments <- lapply(isovalues, FUN = function(x){assign(x, y = animap)})

  # stack the assignment probabilities into a single raster stack
  assignments <- raster::stack(assignments)
    # Transform the assignments into a true probability surface #
  assignIsoprob <- assignments/raster::cellStats(assignments, sum)

  # function to make an odds ratio (likely vs unlikely) assignment
  oddsFun <- function(x,odds = odds){
    predict(smooth.spline(x = cumsum(sort(x)),
                          sort(x),
                          spar = 0.1),(1-odds))$y
  }
    # if odds is left null - use default of 0.67
  if (is.null(odds)){odds <- 0.67}

  # HERE IS WHERE TO IMPLEMENT THE WEIGHTED ASSIGNMENT
  # The weighting scheme
  weight_range <-seq(from = weightRange[1], to = weightRange[2], by = .1)
  weights <-expand.grid(x = 10 ^ weight_range, y = 10 ^ weight_range)
  weights <- rbind(data.frame(x=1,y=0),weights)
  names(weights) <-c("iso_weight", "abun_weight")

sum_weights <- vector('list',nrow(weights))
cat("\n Interating through possible weighted assignments \n")
pb <- utils::txtProgressBar(min = 0, max = nrow(weights), style = 3)
for(i in 1:nrow(weights)){
  utils::setTxtProgressBar(pb, i)
    tempAssign <- (assignIsoprob^weights$iso_weight[i])*
                  (relAbund^weights$abun_weight[i])

    tempAssign <- tempAssign/raster::cellStats(tempAssign,sum)

    matvalsWeight <- raster::rasterToPoints(tempAssign)

    cuts <- apply(matvalsWeight[,3:ncol(matvalsWeight)],2,FUN = oddsFun,odds = odds)

    step1 <- raster::stack(raster::reclassify(tempAssign,cbind(0,cuts,0)))
    step2 <- raster::stack(raster::reclassify(step1,cbind(cuts,1,1)))

    correctAssign <- diag(raster::extract(step2,knownLocs))
    errorRate <- 1-mean(correctAssign)
    areaAssign <- raster::cellStats(raster::area(step2)*step2,sum)

sum_weights[[i]] <- data.frame(isoWeight=log(weights$iso_weight[i],base = 10),
                              abundWeight=log(weights$abun_weight[i],base = 10),
                              error = errorRate,
                              area = mean(areaAssign))
}
close(pb)
bind_weights <- do.call('rbind',sum_weights)

matvalsWeight <- raster::rasterToPoints(assignIsoprob)

cuts <- apply(matvalsWeight[,3:ncol(matvalsWeight)],2,FUN = oddsFun,odds = odds)

step1 <- raster::stack(raster::reclassify(assignIsoprob,cbind(0,cuts,0)))
step2 <- raster::stack(raster::reclassify(step1,cbind(cuts,1,1)))

correctAssign <- diag(raster::extract(step2,knownLocs))
errorRate <- 1-mean(correctAssign)
areaAssign <- raster::cellStats(raster::area(step2)*step2,sum)

sum_weight <- data.frame(isoWeight=1,
                         abundWeight=NA,
                         error = errorRate,
                         area = mean(areaAssign))

iso_performance <- rbind(sum_weight,bind_weights)
iso_performance$area_percent <- iso_performance$area/max(iso_performance$area)

pareto <- function(df){
  df.sorted <- df[with(df,order(df$error, df$area_percent)),]
  front <- df.sorted[which(!duplicated(cummin(df.sorted$area_percent))),]
  return(front)
}
top_assign <- function(front, df){
  baseline <- df[is.na(df$abundWeight),]
  error_base <- baseline$error
  area_base <- baseline$area_percent
  top <- front[(front$error < error_base  & front$area_percent < area_base),]
  return(top)
}

cat("\n finding optimal assignment weights \n")
frontier <-pareto(iso_performance)

top <-top_assign(front = frontier, df = iso_performance)

b <- Sys.time()-a
cat("\n Done: calculation took", b,attributes(b)$units,"\n")
return(structure(list(top = top,
                      frontier = frontier,
                      performance = iso_performance), class = "weightAssign"))

}

