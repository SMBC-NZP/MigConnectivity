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

