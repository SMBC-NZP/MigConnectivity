#' IsoAssign
#' The \code{IsoAssign} function does ....
#'
#' @param isovalues vector of tissue isotope values
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param oddsRatio logical to convert probability into odds ratio
#' @param odds if \code{oddsRatio} the odds ratio to use to set likely and unlikely locations
#' @param singleCellAssign if TRUE generates \code{nSim} single location assignments using a multinomial
#'        distribution using the assignment probability.
#' @param nSim integer specifying how many random samples to draw from a multinomial distribution
#' @param dataFrame logical defining whether to return results as a data.frame. If FALSE (default) returns a
#'        raster layer.
#' @param sppShapefile SpatialPolygon layer defining species range. Assignments are restricted to these
#'        areas.
#' @param assignExtent definition for the extent of the assignment. Can be used in place of \code{sppShapefile} to
#'        limit assignment. Input should follow \code{c(xmin,xmax,ymin,ymax)} in degrees longitude and latitude.
#' @param return which assignment map to return. 'probability' returns assignment probability.
#'     'odds' returns a single map for each isovalue as likely or unlikely assignment. 'population'
#'     the population level assignment (sum of all 'odds').
#' @param element The elemental isotope of interest. Currently the only
#'     elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface if "TRUE" returns surface water values. Defaults is 'FALSE'
#'     which returns the isotopes ratio found in precipitation.
#' @param period The time period of interest. If 'Annual' returns a raster
#'     of mean annual values in precipitation for the \code{element}. If
#'     'GrowingSeason' returns growing season values in precipitation for
#'      \code{element} of interest.
#'
#' @return raster stack, raster layer or data.frame if \code{dataFrame = TRUE} of assignment probabilites
#'
#' @export
#'
#' @example
#' \dontrun{
#' OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
#' OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#' crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#'
#' OVENvals <- read.csv("data-raw/deltaDvalues.csv")
#'
#' a <- Sys.time()
#' b <- isoAssign(isovalues = OVENvals[,2],
#'               isoSTD = 12,
#'               intercept = -10,
#'               slope = 0.8,
#'               oddsRatio = FALSE,
#'               odds = NULL,
#'               SingleCellAssign = TRUE,
#'               nSim = 1000,
#'               dataFrame = FALSE,
#'               sppShapefile = OVENdist,
#'               assignExtent = c(-179,-60,15,89),
#'               return = "sim.cell",
#'               element = "Hydrogen",
#'               surface = FALSE,
#'               period = "Annual")
#' Sys.time()-a}

isoAssign <- function(isovalues,
                      isoSTD,
                      intercept,
                      slope,
                      oddsRatio = FALSE,
                      odds = NULL,
                      SingleCellAssign = FALSE,
                      nSim = NULL,
                      dataFrame = FALSE,
                      sppShapefile = NULL,
                      assignExtent = c(-179,-60,15,89),
                      return = "probability",
                      element = "Hydrogen",
                      surface = FALSE,
                      period = "Annual"){

  # quick input check #
if(!return %in% c("probability","population","odds","sim.cell")){
  stop("return must be either probability,population,odds,sim.cell")}

# download isoscape map
isomap <- getIsoMap(element = element, surface = surface, period = period)

# Series of checks for a species range map inputs
# 1. if sppShapefile == NULL - use extent option
if(is.null(sppShapefile)){
isomap <- raster::crop(isomap,raster::extent(assignExtent))
}
# 2. if sppShapefile provided check that it has a projection defined
#    if not stop - if so, mask the isoscape to range
if(!is.null(sppShapefile) & is.na(raster::crs(sppShapefile))){
  stop("coordinate system needed for sppShapefile")}
if(!is.null(sppShapefile) & (sppShapefile@proj4string@projargs == isomap@crs@projargs)){
isomap <- raster::crop(isomap, sppShapefile)
isomap <- raster::mask(isomap, sppShapefile)
}
# 3. if the projections don't match - project into same as isomap then mask
if(!is.null(sppShapefile) & !(sppShapefile@proj4string@projargs == isomap@crs@projargs)){
sppShapefile <- sp::spTransform(sppShapefile, sp::CRS(isomap@crs@projargs))
isomap <- raster::crop(isomap, sppShapefile)
isomap <- raster::mask(isomap, sppShapefile)
}

# generate a 'feather'/animal isoscape
animap <- raster::calc(isomap, function(x){y <- slope*x+intercept})

# spatially explicit assignment
assign <- function(x,y) {((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))}

assignments <- lapply(isovalues, FUN = function(x){assign(x, y = animap)})

assignments <- raster::stack(assignments)

# Transform the assignments into a true probability surface #
assign2prob <- assignments/raster::cellStats(assignments,sum)

if(dataFrame == TRUE){
assing2probDF <- dataFrame(raster::rasterToPoints(assign2prob))
}

oddsFun <- function(x,odds = odds){
  predict(smooth.spline(x = cumsum(sort(x)),
                        sort(x),
                        spar = 0.1),odds)$y
}

if(oddsRatio == TRUE){
 if(is.null(odds)) odds <- 0.33;

matvals <- raster::values(assign2prob)

cuts <- apply(matvals,2,FUN = oddsFun,odds = odds)

step1 <- raster::reclassify(assign2prob,cbind(0,cuts,0))
step2 <- raster::reclassify(step1,cbind(cuts,1,1))

SamplePop <- sum(step2)

if(dataFrame == TRUE){
  step2DF <- dataFrame(raster::rasterToPoints(step2))
  SamplePopDF <- dataFrame(raster::rasterToPoints(SamplePop))
}
}

if(SingleCellAssign == TRUE){
  if(is.null(nSim)) nSim<-1000;
xysim <- array(NA,c(nSim,2,raster::nlayers(assign2prob)))
#dimnames(xysim)[[1]] <- 1:nSim
dimnames(xysim)[[2]] <- c("Longitude","Latitude")
dimnames(xysim)[[3]] <- names(assign2prob)

matvals <- raster::rasterToPoints(assign2prob)

for(i in 1:nSim){
  multidraw <- apply(matvals[,3:ncol(matvals)],2,FUN = function(x){rmultinom(n = 1, size = 1, prob = x)})
  xysim[i,1,] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],1]
  xysim[i,2,] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],2]
}
}

if(return == "probability" & dataFrame == FALSE){return(assign2prob)}
if(return == "probability" & dataFrame == TRUE){return(assing2probDF)}
if(return == "odds" & dataFrame == FALSE){return(step2)}
if(return == "odds" & dataFrame == TRUE){return(step2DF)}
if(return == "population" & dataFrame == FALSE){return(SamplePop)}
if(return == "population" & dataFrame == TRUE){return(SamplePopDF)}
if(return == "sim.cell"){return(xysim)}
}
