#' isoAssign
#' The \code{isoAssign} function generates origin assignments using stable-hydrogen isotopes in tissue. The function generates
#' a probability surface of origin assignment from a vector of stable-isotope values for each animal/sample of interest.
#' Probabilistic assignments are constructed by first converting observed stable-isotope ratios (isoscape) in either precipitation or surface
#' waters into a 'tissuescape' using a user-provided intercept, slope and standard deviation. See \href{http://journals.plos.org.mutex.gmu.edu/plosone/article?id=10.1371/journal.pone.0035137}{Hobson et. al. 2012}
#'
#'
#' @param isovalues vector of tissue isotope values
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param oddsRatio logical to convert probability into odds ratio
#' @param odds if \code{oddsRatio} the odds ratio to use to set likely and unlikely locations
#' @param singleCellAssign if TRUE generates \code{nSim} single location assignments using a multinomial
#'        distribution using the assignment probability.
#' @param restrict2Likely if \code{TRUE} restricts locations to fall within the most probable assignment
#'        locations if \code{singleCellAssign = TRUE}. Only used when \code{oddsRatio = TRUE}.
#' @param nSamples integer specifying how many random samples to draw from a multinomial distribution.
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
#'               restrict2Likely = TRUE,
#'               nSamples = 1000,
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
                      restrict2Likely = TRUE,
                      nSamples = NULL,
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
  if(is.null(nSamples)) nSamples<-1000;
xysim <- array(NA,c(nSamples,2,raster::nlayers(assign2prob)))
#dimnames(xysim)[[1]] <- 1:nSamples
dimnames(xysim)[[2]] <- c("Longitude","Latitude")
dimnames(xysim)[[3]] <- names(assign2prob)

# converts raster to matrix of XY then probs
matvals <- raster::rasterToPoints(assign2prob)

# Restrict random point estimates to 'likely' origin area #
if(restrict2Likely && oddsRatio){
binaryLikely <- raster::rasterToPoints(step2)
matvals[,3:ncol(matvals)] <- matvals[,3:ncol(matvals)]*binaryLikely[,3:ncol(matvals)]
matvals[matvals==0]<-NA
}

# This draws samples 1 per animal nSamples number of times and fills the xysim with x,y coords
for(i in 1:nSamples){
  multidraw <- apply(matvals[,3:ncol(matvals)],2,FUN = function(x){rmultinom(n = 1, size = 1, prob = x)})
  xysim[i,1,] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],1]
  xysim[i,2,] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],2]
}
}
#What should be returned
if(return == "probability" & dataFrame == FALSE){return(assign2prob)}
if(return == "probability" & dataFrame == TRUE){return(assing2probDF)}
if(return == "odds" & dataFrame == FALSE){return(step2)}
if(return == "odds" & dataFrame == TRUE){return(step2DF)}
if(return == "population" & dataFrame == FALSE){return(SamplePop)}
if(return == "population" & dataFrame == TRUE){return(SamplePopDF)}
if(return == "sim.cell"){return(xysim)}
}



#' Get Isoscape map
#' The \code{getIsoMap} function downloads predicted isoscape maps from
#'  \url{http://wateriso.utah.edu/waterisotopes/}. The function first checks
#' whether the isoscapes are located within the current working directory
#' \code{getwd()}. If a local copy of the isoscape is found, it's read into
#' the environment. If not, the isoscape is downloaded and imported
#' as a raster.
#'
#' @param element The elemental isotope of interest. Currently the only
#'     elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface if "TRUE" returns surface water values. Defaults is 'FALSE'
#'     which returns the isotopes ratio found in precipitation.
#' @param period The time period of interest. If 'Annual' returns a raster
#'     of mean annual values in precipitation for the \code{element}. If
#'     'GrowingSeason' returns growing season values in precipitation for
#'      \code{element} of interest.
#'
#' @return raster layer
#'
#' @export
#'
#' @examples
#' \dontrun{
#' getIsoMap(element = "Hydrogen", period = "Annual")
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
        download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/Surface_H.zip",
                      destfile = tf,
                      quiet = TRUE,
                      extra = getOption("download.file.extra"))

        # unzip the downloaded file #
        unzip(zipfile = tf,
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
        download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/Surface_O.zip",
                      destfile = tf,
                      quiet = TRUE,
                      extra = getOption("download.file.extra"))

        # unzip the downloaded file #
        unzip(zipfile = tf,
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
      download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/AnnualD.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      unzip(zipfile = tf,
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
      download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/GSD.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      unzip(zipfile = tf,
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
      download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/AnnualO.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      unzip(zipfile = tf,
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
      download.file(url = "http://wateriso.utah.edu/waterisotopes/media/ArcGrids/GSO.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      unzip(zipfile = tf,
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
