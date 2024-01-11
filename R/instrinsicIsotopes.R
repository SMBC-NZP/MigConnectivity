#' Generate probabilistic isotope assignments
#'
#' The \code{isoAssign} function generates origin assignments using
#' stable-hydrogen isotopes in tissue. The function generates a probability
#' surface of origin assignment from a vector of stable-isotope values for each
#' animal/sample of interest. Probabilistic assignments are constructed by first
#' converting observed stable-isotope ratios (isoscape) in either precipitation
#' or surface waters into a 'tissuescape' using a user-provided intercept, slope
#' and standard deviation. See
#' \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0035137}{Hobson et. al. (2012)}.
#'
#'
#' @param isovalues vector of tissue isotope values
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param odds odds ratio to use to set likely and unlikely locations defaults
#'        to 0.67
#' @param restrict2Likely if \code{TRUE} restricts locations to fall within the
#'        'likely' assignment locations
#' @param nSamples integer specifying how many random samples to draw from a
#'        multinomial distribution.
#' @param sppShapefile A polygon spatial layer (sf - MULTIPOLYGON) defining
#'        species range. Assignments are restricted to these areas.
#' @param relAbund raster (\code{SpatRast}) with relative abundance (must match
#'        extent of isotope assignment)
#' @param isoWeight weighting value to apply to isotope assignment
#' @param abundWeight weighting value to apply to relative abundance prior
#' @param population vector identifying location where animal was captured.
#'        Same order as \code{isovalues}
#' @param assignExtent definition for the extent of the assignment. Can be used
#'        in place of \code{sppShapefile} to limit assignment. Input should
#'        follow \code{c(xmin,xmax,ymin,ymax)} in degrees longitude and
#'        latitude
#' @param element The elemental isotope of interest. Currently the only
#'        elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface DEPRECATED function no longer returns surface water values.
#'     Default is 'FALSE' which returns the precipitation isotopes ratio.
#' @param period The time period of interest. If 'Annual' returns a raster
#'        of mean annual values in precipitation for the \code{element}. If
#'        'GrowingSeason' returns growing season values in precipitation for
#'        \code{element} of interest
#' @param seed numeric value fed to \code{set.seed} for random number
#'        generation. Default = NULL
#' @param verbose takes values 0, 1 (default) or 2. 0 prints no output during
#'        run. 1 prints a message detailing where in the process the function
#'        is. 2 prints the animal currently being sampled.
#' @param generateSingleCell if 'TRUE' generates a single origin location using
#'        the posterior assignment distribution - this takes a while to run.
#'        If 'FALSE' (default), no coordinates are generated.
#' @param mapDirectory Directory to save/read isotope map from. Can use relative
#'     or absolute addressing. The default value (NULL) downloads to a temporary
#'     directory, so we strongly recommend changing this from the default unless
#'     you're sure you're not going to need these data more than once.
#'
#' @return returns an \code{isoAssign} object containing the following:
#'  \describe{
#'   \item{\code{probassign}}{SpatRast stack of individual probabilistic assignments}
#'   \item{\code{oddsassign}}{SpatRast stack that includes likely vs unlikely origin for each animal}
#'   \item{\code{popassign}}{a SpatRast for population level assignment (sum of \code{oodsassign} if \code{population} = NULL).
#'   If \code{population} is a vector then returns a raster stack for each unique \code{population} provided}
#'   \item{\code{probDF}}{data.frame of individual probability surfaces}
#'   \item{\code{oddsDF}}{data.frame of likely vs unlikely surfaces}
#'   \item{\code{popDF}}{data.frame of population level assignment}
#'   \item{\code{SingeCell}}{array of coordinates (longitude,latitude) for
#'   single cell assignment}
#'   \item{\code{targetSites}}{\code{sf - MULTIPOLYGON} layer representing
#'   isotope bands equivalent to \code{isoSTD}}
#'   \item{\code{RandomSeed}}{the RNG seed used when generating locations from
#'   the multinomial distribution}
#'   }
#'
#' @seealso{\code{\link{weightAssign}}}
#' @export
#'
#' @examples
#' \donttest{
#' extensions <- c("shp", "shx", "dbf", "sbn", "sbx")
#' tmp <- tempdir()
#' for (ext in extensions) {
#' download.file(paste0(
#'               "https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity",
#'                      "/master/data-raw/Spatial_Layers/OVENdist.",
#'                      ext),
#'               destfile = paste0(tmp, "/OVENdist.", ext), mode = "wb")
#' }
#' OVENdist <- sf::st_read(paste0(tmp, "/OVENdist.shp"))
#' OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#' sf::st_crs(OVENdist) <- sf::st_crs(4326)
#'
#' download.file(paste0(
#'   "https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity",
#'                      "/master/data-raw/deltaDvalues.csv"),
#'               destfile = paste0(tmp, "/deltaDvalues.csv"))
#' OVENvals <- read.csv(paste0(tmp, "/deltaDvalues.csv"))
#'
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
#'               period = "GrowingSeason") # this setting for demonstration only
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
                      verbose=1,
                      generateSingleCell = FALSE,
                      mapDirectory = NULL) {
  # force verbose to default when outside specified range.
  if(!(verbose %in% c(0,1,2))){
    verbose = 1
  }

  # download isoscape map
  isomap <- getIsoMap(element = element, surface = surface, period = period,
                      mapDirectory = mapDirectory)

  # 1. if sppShapefile == NULL - use extent option
  if(is.null(sppShapefile)){
    isomap <- terra::crop(isomap,terra::ext(assignExtent))
  } else {
  # Series of checks for a species range map inputs
  # 2. if sppShapefile provided check that it has a projection defined
  #    if not stop - if so, mask the isoscape to range
    if(inherits(sppShapefile,"SpatVector") & is.na(terra::crs(sppShapefile))){stop("coordinate system needed for sppShapefile")}
    if(inherits(sppShapefile,"sf") & is.na(sf::st_crs(sppShapefile))){stop("coorindate system needed for sppShapfile")}
    #if(is.na(terra::crs(sppShapefile))){
    #  stop("coordinate system needed for sppShapefile")
    #}
  # 3. if the projections don't match - project into same as isomap then mask
    # quick check
    if(inherits(sppShapefile, "SpatialPolygons") ||
       inherits(sppShapefile, "SpatialPolygonsDataFrame")){
      warning("The sp package is no longer supported or maintained - consider using the sf package instead")
      #if sp object convert to sf
      sppShapefile <- sf::st_as_sf(sppShapefile)
      if(!identical(sf::st_crs(sppShapefile),sf::st_crs(4326))){
        sppShapefile <- sf::st_transform(sppShapefile, 4326)
      }
    }
    if(class(sppShapefile)[1] %in% "sf"){
      if(!identical(sf::st_crs(sppShapefile),sf::st_crs(4326))){
        sppShapefile <- sf::st_transform(sppShapefile, 4326)
      }
    }

    if(verbose>0){
      cat("\n Restricting possible assignments to species distribution \n")
    }

    sppShapefile$INOUT<-1
    if (!inherits(sppShapefile, "SpatVector")) {
      sppShapefile <- terra::vect(sppShapefile)
    }
    isomap <- terra::project(isomap, terra::crs(sppShapefile))
    # mask the isomap to sppShapefile
    isomap <- terra::crop(isomap, sppShapefile)
    isomap <- terra::mask(isomap, sppShapefile)
  }
  if(!is.null(relAbund) && inherits(relAbund,"RasterLayer")){
    warning("The raster package relies on rgdal which is no longer supported - consider using the terra package instead \n")
  }
  if(!is.null(relAbund) && !inherits(relAbund,"SpatRaster")){
    #convert from raster to terra
    relAbund <- terra::rast(relAbund)
  }
  # if isomap and relAbund don't have the same resolution and/or extent
  # change to relAbund to match isomap
  if(!is.null(relAbund) && !terra::compareGeom(isomap,relAbund)){
    # project to match isomap
    relAbund <- terra::project(relAbund,isomap)
    # re-scale to ensure sums to 1
    relAbund <- propSpatRaster(relAbund)
  }
  # generate a 'feather'/animal isoscape
  animap <- terra::app(isomap, fun = function(x){y <- slope*x+intercept})

  # generate targetSites - seq from min to max values by isoSTD
  isocut <- terra::classify(animap,
                            rcl= seq(from = unlist(terra::global(animap,
                                                                 fun = "min",
                                                                 na.rm = TRUE)),
                                     to = unlist(terra::global(animap,
                                                               fun = "max",
                                                               na.rm = TRUE)),
                                     by = isoSTD))

  # use those cuts to make polygons
  # added makeValid - just in case but prevents errors -
  # we could add an on error clause
  targetSites_iso <- terra::makeValid(terra::as.polygons(isocut))

  # rename the targetSites to simplify output
  targetSites_iso <- targetSites_iso[,1]
  names(targetSites_iso) <- c("targetSite")

  # ensure that the targetSites are aggregated
  targetSites <- terra::aggregate(targetSites_iso, by = "targetSite", dissolve = TRUE)

  terra::crs(targetSites) <- terra::crs(isomap)

  #if sppShapefile !NULL then clip targetSites to distribution
  if(!is.null(sppShapefile)){
    targetSites <- terra::intersect(targetSites,
                                    sppShapefile)

    targetSites <- terra::aggregate(targetSites_iso, by = "targetSite", dissolve = TRUE)
  }

  #targetSites <- rgeos::gUnaryUnion(targetSites, id=targetSites$targetSite)
  # spatially explicit assignment
  assign <- function(x,y) {
    ((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))
  }

  # apply the assignment function to all input values
  if(verbose>0){
    cat("\n Generating probabilistic assignments \n")
  }

  assignments <- lapply(isovalues, FUN = function(x){assign(x, y = animap)})

  # stack the assignment probabilities into a single raster stack
  assignments <- terra::rast(assignments)
  terra::set.names(assignments, paste0("lyr.", 1:terra::nlyr(assignments)))
  # Transform the assignments into a true probability surface #
  assign2prob <- propSpatRaster(assignments)
  # Weighted Assignments
  if(inherits(relAbund,"spatRast") && is.null(isoWeight) && is.null(abundWeight)){
    if(verbose>0){
      cat("\n Creating posterior assignments where isotope & abundance have equal weight \n")
    }
    test <- assign2prob*relAbund
    assign2prob <- propSpatRaster(test)
  }
  else if(inherits(relAbund,"spatRast") && !is.null(isoWeight) && !is.null(abundWeight)){
    if(verbose>0){cat("\n Creating weighted posterior assignments \n")}
    if (is.null(isoWeight))
      isoWeight <- 0
    if (is.null(abundWeight))
      abundWeight <- 0
    isoWeight <- 10^isoWeight
    abundWeight <- 10^abundWeight
    assign2prob <- (assign2prob^isoWeight)*(relAbund^abundWeight)
    assign2prob <- propSpatRaster(assign2prob)
  }

  # Create a dataframe with XY coords and probabilites for each animal
  assign2probDF <- data.frame(terra::as.points(assign2prob))


  # if odds is left null - use default of 0.33
  if (is.null(odds)){odds <- 0.67}

  # extract values from the probability assignment
  matvals <- terra::as.data.frame(assign2prob, xy = TRUE)
  # XY coords of raster
  matvalsXY <- matvals[,1:2]

  # drop XY from matvals
  matvals <- matvals[,-(1:2)]

  if(verbose>0){
    cat("\n Generating likely vs unlikely assignments \n")
  }
  # apply the odds function

  cuts <- apply(matvals,2,FUN = oddsFun,odds = odds)
  # reclassify the rasters based on likely v unlikely

  step1 <- mapply(FUN = function(x,y){terra::classify(assign2prob[[x]],rcl=cbind(0,y,0))},
                  x = 1:terra::nlyr(assign2prob),
                  y = cuts)
  step1 <- terra::rast(step1)
  step2 <- mapply(FUN = function(x,y){terra::classify(step1[[x]], rcl=cbind(y,1,1))},
                  x = 1:terra::nlyr(step1),
                  y = cuts)
  step2 <- terra::rast(step2)


  # convert to data frame
  step2DF <- terra::as.data.frame(step2, xy = TRUE)
  LikelyUnlikely <- as.matrix(step2DF)


  if(is.null(population)){
    # Return the population level odds assignment - i.e, how many animals
    SamplePop <- sum(step2)
    # convert to data frame
    SamplePopDF <- data.frame(terra::as.points(SamplePop))
  }
  else {
    nPop <- length(unique(population))
    Pops <- unique(population)
    POPs <- POPSdf <- vector("list",nPop)
    for(p in 1:nPop){
      aniPop <- which(population == Pops[p])
      pop_p <- step2[[aniPop]]
      POPs[[p]] <- sum(pop_p)
      POPSdf[[p]] <- data.frame(terra::as.points(POPs[[p]]))
    }
    SamplePop <- terra::rast(POPs)
    names(SamplePop)<-Pops
    SamplePopDF <- do.call('cbind', POPSdf)
  }

  # SINGLE CELL PROBABILITY ASSIGNMENTS - this makes MC possible with isotopes
  if (is.null(nSamples)){ nSamples <- 1000}
  if(is.null(seed)){seed <- as.numeric(Sys.time())}

  # Set random number generator to replicate results
  set.seed(seed)

  # generate empty array to fill with locations
  # make a simulated array twice the size to weed out locations
  # that fall outside of distribution
  if(generateSingleCell){

  xysimulation <- array(NA,c(nSamples+floor((nSamples/2)),2,terra::nlyr(assign2prob)))

  # give names for sf to convert down the line
  dimnames(xysimulation)[[2]] <- c("Longitude","Latitude")

  xysim <- array(NA, c(nSamples, 2, terra::nlyr(assign2prob)))
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
} # end generateSingleCell

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
  else {
    xysim <- NULL
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
#'  \url{https://wateriso.utah.edu/waterisotopes/}. The function first checks
#' whether the isoscapes are located within the directory
#' \code{mapDirectory}. If a local copy of the isoscape is found, it's read into
#' the environment. If not, the isoscape is downloaded and imported
#' as a raster.
#'
#'
#' @param element The elemental isotope of interest. Currently the only
#'     elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface DEPRECATED function no longer returns surface water values.
#'     Default is 'FALSE' which returns the precipitation isotopes ratio.
#' @param period The time period of interest. If 'Annual' (default) returns a
#'     raster of mean annual values in precipitation for the \code{element}. If
#'     'GrowingSeason' returns growing season values in precipitation for
#'      \code{element} of interest.
#' @param mapDirectory Directory to save/read isotope map from. Can use relative
#'     or absolute addressing. The default value (NULL) downloads to a temporary
#'     directory, so we strongly recommend changing this from the default unless
#'     you're sure you're not going to need these data more than once.
#'
#' @return returns a global \code{RasterLayer} (resolution = 0.333'x0.3333')
#'     object for the \code{element} and \code{period} of interest
#'
#' @export
#' @examples
#' \donttest{
#'   map <- getIsoMap(element = "Hydrogen", period = "GrowingSeason")
#' }

getIsoMap<-function(element = "Hydrogen", surface = FALSE, period = "Annual",
                    mapDirectory = NULL){

  # and read into R as raster - otherwise read MAD into R
  if(!(element %in% c("Hydrogen", "Oxygen")))
    stop("element must be either Hydrogen or Oxygen")

  if(!(period %in% c("Annual", "GrowingSeason")))
    stop("period must be either Annual or GrowingSeason")

  # if "/AnnualD" isn't in working directory somewhere download MAD from website
  # Download Mean Annual Deuterium values from wateriso.utah.edu #
  if (is.null(mapDirectory)) {
    haveIsoMap <- ""
    destDirectory <- tempdir()
  }
  else {
    mapDirectory <- path.expand(mapDirectory)
    haveIsoMap <- list.dirs(path = mapDirectory, recursive = TRUE)
    destDirectory <- ifelse(surface || period == "GrowingSeason",
                            paste0(mapDirectory, "/GlobalPrecipGS"),
                            paste0(mapDirectory, "/GlobalPrecip"))
  }

  # Create temporary file location #
  tf <- tempfile(pattern = "file", fileext = "")

  if(surface){
    cat(paste("This function is no longer set up to retrieve surface data;",
              "providing rainfall isotope data\n"))
  }
  if(element == "Hydrogen" & period == "Annual"){

    if(!(destDirectory %in% haveIsoMap)){

      cat("NOTE: The Global H2 precipitation file is large (>600 MB) and may take a while to download\n")
      oldTO <- getOption("timeout")
      options(timeout = max(700, oldTO))


      # Download the file #
      utils::download.file(url = "https://wateriso.utah.edu/waterisotopes/media/ArcGrids/GlobalPrecip.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      options(timeout = oldTO)
      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
                  files = NULL,
				  list = FALSE,
				  overwrite = TRUE,
                  junkpaths = FALSE,
				  exdir = destDirectory,
				  unzip = "internal",
                  setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      m_a_d <- terra::rast(paste0(destDirectory, "/d2h_MA.tif"))
    }else{
      dir <- grep(haveIsoMap, pattern = "/GlobalPrecip$", value = TRUE)
      m_a_d <- terra::rast(paste0(dir, "/d2h_MA.tif"))
    }
    names(m_a_d)<-"MeanAnnualD"
    return(m_a_d)
  }

  if(element == "Hydrogen" & period == "GrowingSeason"){
    if(!(destDirectory %in% haveIsoMap)){

      # Download the file #
      utils::download.file(
        url = "https://wateriso.utah.edu/waterisotopes/media/ArcGrids/GlobalPrecipGS.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
                   files = NULL,
				   list = FALSE,
				   overwrite = TRUE,
                   junkpaths = FALSE,
				   exdir = destDirectory,
				   unzip = "internal",
				   setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      g_s_d <- terra::rast(paste0(destDirectory, "/d2h_GS.tif"))
    }else{
      dir <- grep(haveIsoMap, pattern = "/GlobalPrecipGS$", value = TRUE)
      g_s_d <- terra::rast(paste0(dir, "/d2h_GS.tif"))
    }
    g_s_d[g_s_d == -999]<-NA
    names(g_s_d)<-"GrowingSeasonD"
    return(g_s_d)
  }
  if(element == "Oxygen" & period == "Annual"){
    if(!(destDirectory %in% haveIsoMap)){

      cat("NOTE: The Global precipitation file is large (>600 MB) and may take a while to download \n")
      oldTO <- getOption("timeout")
      options(timeout = max(700, oldTO))

      # Download the file #
      utils::download.file(
        url = "https://wateriso.utah.edu/waterisotopes/media/ArcGrids/GlobalPrecip.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      options(timeout = oldTO)

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
	               files = NULL,
				   list = FALSE,
				   overwrite = TRUE,
                   junkpaths = FALSE,
				   exdir = destDirectory,
				   unzip = "internal",
                   setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      m_a_o <- terra::rast(paste0(destDirectory, "/d18o_MA.tif"))
    }else{
      dir <- grep(haveIsoMap, pattern = "/GlobalPrecip$", value = TRUE)
      m_a_o <- terra::rast(paste0(dir, "/d18o_MA.tif"))
    }
    names(m_a_o)<-"MeanAnnualO"
    return(m_a_o)
  }

  if(element == "Oxygen" & period == "GrowingSeason"){
    if(!(destDirectory %in% haveIsoMap)){

      # Download the file #
      utils::download.file(url = "https://wateriso.utah.edu/waterisotopes/media/ArcGrids/GlobalPrecipGS.zip",
                    destfile = tf,
                    quiet = TRUE,
                    extra = getOption("download.file.extra"))

      # unzip the downloaded file #
      utils::unzip(zipfile = tf,
	               files = NULL,
 				   list = FALSE,
 				   overwrite = TRUE,
            junkpaths = FALSE,
			exdir = destDirectory,
			unzip = "internal",
            setTimes = FALSE)

      # Delete zipped folder #
      file.remove(tf)

      g_s_o <- terra::rast(paste0(destDirectory,"/d18o_GS.tif"))
    }else{
      dir <- grep(haveIsoMap, pattern = "/GlobalPrecipGS$", value = TRUE)
      g_s_o <- terra::rast(paste0(dir, "/d18o_GS.tif"))
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
#' \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0035137}{Hobson et. al. (2012)}.
#'
#' @param knownLocs matrix of capture locations of the same length as
#'        \code{isovalues}
#' @param isovalues vector of tissue isotope values from known locations
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param odds odds ratio to use to set likely and unlikely locations defaults
#'        to 0.67
#' @param relAbund raster layer of relative abundance that sums to 1.
#' @param weightRange vector of length 2 within minimum and maximum values to
#'        weight isotope and relative abundance. Default = c(-1,1)
#' @param sppShapefile A polygon spatial layer (sf - MULTIPOLYGON) defining
#'        species range. Assignments are restricted to these areas.
#' @param assignExtent definition for the extent of the assignment. Can be used
#'        in place of \code{sppShapefile} to limit assignment. Input should
#'        follow \code{c(xmin,xmax,ymin,ymax)} in degrees longitude and latitude.
#' @param element The elemental isotope of interest. Currently the only
#'        elements that are implemented are 'Hydrogen' (default) and 'Oxygen'
#' @param surface DEPRECATED function no longer returns surface water values.
#'     Default is 'FALSE' which returns the precipitation isotopes ratio.
#' @param period The time period of interest. If 'Annual' returns a raster
#'        of mean annual values in precipitation for the \code{element}. If
#'        'GrowingSeason' returns growing season values in precipitation for
#'        \code{element} of interest.
#' @param mapDirectory Directory to save/read isotope map from. Can use relative
#'     or absolute addressing. The default value (NULL) downloads to a temporary
#'     directory, so we strongly recommend changing this from the default unless
#'     you're sure you're not going to need these data more than once.
#'
#' @return returns an \code{weightAssign} object containing the following:
#'   \describe{
#'    \item{\code{top}}{data.frame with the optimal weightings}
#'    \item{\code{frontier}}{data.frame with values that fall along the Pareto
#'      frontier}
#'    \item{\code{performance}}{data.frame with error rate and assignment area
#'      for each weight combination}
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
#' Evolution 3: 3847-3855. \doi{10.1002/ece3.2605}
#' @export
#'
#' @examples
#' \donttest{
#' extensions <- c("shp", "shx", "dbf", "sbn", "sbx")
#' tmp <- tempdir()
#' for (ext in extensions) {
#' download.file(paste0(
#'               "https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity",
#'                      "/master/data-raw/Spatial_Layers/OVENdist.",
#'                      ext),
#'               destfile = paste0(tmp, "/OVENdist.", ext), mode = "wb")
#' }
#' OVENdist <- sf::st_read(paste0(tmp, "/OVENdist.shp"))
#' OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#' sf::st_crs(OVENdist) <- sf::st_crs(4326)
#'
#' download.file(paste0("https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity",
#'                      "/master/data-raw/deltaDvalues.csv"),
#'               destfile = paste0(tmp, "/deltaDvalues.csv"))
#' OVENvals <- read.csv(paste0(tmp, "/deltaDvalues.csv"))
#'
#' HBEFbirds <- OVENvals[grep("NH",OVENvals[,1]),]
#'
#' # Create a spatial object of known capture sites
#' knownLocs <- sf::st_as_sf(data.frame(Long = rep(-73,nrow(HBEFbirds)),
#'                                     Lat = rep(43,nrow(HBEFbirds))),
#'                          coords = c("Long","Lat"),
#'                          crs = 4326)
#'
#' #Get OVEN abundance from BBS estimates and read into R #
#' utils::download.file("https://www.mbr-pwrc.usgs.gov/bbs/ra15/ra06740.zip",
#'                      destfile = paste0(tmp, "/oven.zip"))
#' utils::unzip(paste0(tmp, "/oven.zip"), exdir = tmp)
#' oven_dist <- sf::st_read(paste0(tmp, "/ra06740.shp"))
#'
#' # Empty raster with the same dimensions as isoscape and Ovenbird distribution
#'
#' # We do this manually here but the weightedAssign function has been updated
#' # to ensure the isoscape and abundance rasts have the same extent using
#' # resampling to match  relAbund to the isoscape.
#' r <- terra::rast(nrow = 331, ncol = 870,
#'                  res = c(0.0833333, 0.0833333),
#'                  xmin = -125.1667, xmax = -52.66672,
#'                  ymin = 33.49995, ymax = 61.08327,
#'                  crs = sf::st_crs(4326)$wkt)
#'
#' # rasterize the polygons from BBS - this is not needed if working with a
#' # rasterized surface
#' relativeAbun<-terra::rasterize(terra::vect(sf::st_transform(oven_dist,4326)),
#'                                r,
#'                                field = "RASTAT")
#'
#' relativeAbund <- relativeAbun/terra::global(relativeAbun, sum,
#'                                             na.rm = TRUE)$sum
#'
#'
#' BE <- weightAssign(knownLocs = knownLocs,
#'                    isovalues = HBEFbirds[,2],
#'                    isoSTD = 12,
#'                    intercept = -10,
#'                    slope = 0.8,
#'                    odds = 0.67,
#'                    relAbund = relativeAbund,
#'                    weightRange = c(-1, 1),
#'                    sppShapefile = OVENdist,
#'                    assignExtent = c(-179,-60,15,89),
#'                    element = "Hydrogen",
#'                    period = "Annual")
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
                         period = "Annual",
                         mapDirectory = NULL){
a <- Sys.time()
  # Check to make sure knownLocs are added by user
  if(is.null(knownLocs)){
    stop("Known locations are needed to determine weighted assignemnts")}
  if(nrow(knownLocs)!=length(isovalues)){
    stop("A known location (knownLocs) is needed for each isotope value (isovalues)")
  }
  # download isoscape map
  isomap <- getIsoMap(element = element, surface = surface, period = period,
                      mapDirectory = mapDirectory)

  # 1. if sppShapefile == NULL - use extent option
  if(is.null(sppShapefile)){
    isomap <- terra::crop(isomap,terra::ext(assignExtent))
  }

  # Series of checks for a species range map inputs
  if(!is.null(sppShapefile)){
    # 2. if sppShapefile provided check that it has a projection defined
    #    if not stop - if so, mask the isoscape to range
    if(is.na(sf::st_crs(sppShapefile))){
      stop("coordinate system needed for sppShapefile")}
    # 3. if the projections don't match - project into same as isomap then mask
    if(sf::st_crs(sppShapefile)!= terra::crs(isomap)){
      sppShapefile <- sf::st_transform(sppShapefile, terra::crs(isomap))
    }

    cat("\n Restricting possible assignments to species distribution \n")

    sppShapefile$INOUT<-1

    # mask the isomap to sppShapefile
    isomap <- terra::mask(isomap, terra::vect(sppShapefile))
    isomap <- terra::crop(isomap, terra::vect(sppShapefile))
  }

  # generate a 'feather'/animal isoscape
  animap <- slope*isomap+intercept


  # spatially explicit assignment
  assign <- function(x,y) {((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))}

  # apply the assignment function to all input values
  cat("\n Generating probabilistic assignments \n")

  assignments <- do.call(c,lapply(isovalues, FUN = function(x){assign(x, y = animap)}))

  # stack the assignment probabilities into a single raster stack
  # assignments <- terra::rast(assignments)

  # Transform the assignments into a true probability surface #
  assignIsoprob <- propSpatRaster(assignments)
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

# good pratice to check that the extents match -
# if not resample the relaAbun fed from the user

if(!terra::compareGeom(assignIsoprob,relAbund, stopOnError = FALSE)){
  relAbund <- terra::resample(relAbund, assignIsoprob)
}

cat("\n Iterating through possible weighted assignments \n")
pb <- utils::txtProgressBar(min = 0, max = nrow(weights), style = 3)
for(i in 1:nrow(weights)){
  utils::setTxtProgressBar(pb, i)

  tempAssign <- (assignIsoprob^weights$iso_weight[i])*
                  (relAbund^weights$abun_weight[i])

    tempAssign <- propSpatRaster(tempAssign)

    matvalsWeight <- terra::as.data.frame(tempAssign)

    cuts <- apply(matvalsWeight, 2, FUN = oddsFun,odds = odds)

    step1 <- c(terra::classify(tempAssign,rcl=cbind(0,cuts,0)))
    step2 <- c(terra::classify(step1,rcl=cbind(cuts,1,1)))

    correctAssign <- diag(as.matrix(terra::extract(step2,terra::vect(knownLocs), ID = FALSE)))

    errorRate <- 1-mean(correctAssign)
    areaAssign <- terra::global(terra::cellSize(step2[[1]])*step2,fun = "sum",na.rm = TRUE)

sum_weights[[i]] <- data.frame(isoWeight=log(weights$iso_weight[i],base = 10),
                              abundWeight=log(weights$abun_weight[i],base = 10),
                              error = errorRate,
                              area = mean(areaAssign$sum))
}
close(pb)
bind_weights <- do.call('rbind',sum_weights)

# replacing here because kept getting [readstart] warning for some reason
# I think its because the rast was in memory and accessed in the function above
# and it doesn't like accessing it again in a different environment.
# I'm not exactly sure
matvalsWeight2 <- terra::as.data.frame(propSpatRaster(assignments))

cuts <- apply(matvalsWeight2,2,FUN = oddsFun,odds = odds)

step1 <- c(terra::classify(assignIsoprob,rcl=cbind(0,cuts,0)))
step2 <- c(terra::classify(step1,rcl=cbind(cuts,1,1)))

correctAssign <- diag(as.matrix(terra::extract(step2,terra::vect(knownLocs),ID = FALSE)))
errorRate <- 1-mean(correctAssign)
areaAssign <- terra::global(terra::cellSize(step2[[1]])*step2,fun = "sum",na.rm = TRUE)

sum_weight <- data.frame(isoWeight=1,
                         abundWeight=NA,
                         error = errorRate,
                         area = mean(areaAssign$sum,na.rm = TRUE))

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

