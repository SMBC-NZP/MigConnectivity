#' IsoAssign
#' The \code{IsoAssign} function does ....
#'
#' @param isovalues vector of tissue isotope values
#' @param isoSTD standard deviation from calibration
#' @param intercept intercept value from calibration
#' @param slope value from calibration
#' @param oddsRatio logical to convert probability into odds ratio
#' @param odds if \code{oddsRatio} the odds ratio to use to set likely and unlikely locations
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
#' @return raster stack or raster layer
#'
#' @export
#'
#' @example
#' \dontrun{
#' z <- IsoAssign(isovalues = runif(n = 10,-150,-60),
#'               isoSTD = 12,
#'               intercept = -2,
#'               slope = 1,
#'               oddsRatio = TRUE,
#'               odds = 0.33,
#'               return = "population",
#'               element = "Hydrogen",
#'               surface = FALSE,
#'               period = "Annual")}
#'
IsoAssign <- function(isovalues,
                      isoSTD,
                      intercept,
                      slope,
                      oddsRatio = FALSE,
                      odds = NULL,
                      SingleCellAssign = FALSE,
                      nSim = NULL,
                      Data.Frame = FALSE,
                      assignExtent = c(-179,-60,15,89),
                      return = "probability",
                      element = "Hydrogen",
                      surface = FALSE,
                      period = "Annual"){
  # quick input check #
if(!return %in% c("probability","population","odds","sim.cell")){
  stop("return must be either probability,population,odds,sim.cell")}

isomap <- getIsoMap(element = element, surface = surface, period = period)

isomap <- raster::crop(isomap,raster::extent(assignExtent))

animap <- raster::calc(isomap, function(x){y <- slope*x+intercept})

# spatially explicit assignment
assign <- function(x,y) {((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))}

assignments <- lapply(isovalues, FUN = function(x){assign(x, y = animap)})

assignments <- raster::stack(assignments)

# Transform the assignments into a true probability surface #
assign2prob <- assignments/raster::cellStats(assignments,sum)

if(Data.Frame == TRUE){
assing2probDF <- data.frame(raster::rasterToPoints(assign2prob))
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

if(Data.Frame == TRUE){
  step2DF <- data.frame(raster::rasterToPoints(step2))
  SamplePopDF <- data.frame(raster::rasterToPoints(SamplePop))
}
}

if(SingleCellAssign == TRUE){
  if(!is.null(nSim)) nSim<-1000;
xysim <- array(NA,c(nSim,2,raster::nlayers(assign2prob)))
dimnames(xysim)[[1]] <- 1:nSim
dimnames(xysim)[[2]] <- c("Longitude","Latitude")
dimnames(xysim)[[3]] <- names(assign2prob)

matvals <- raster::rasterToPoints(assign2prob)

for(i in 1:nSim){
  multidraw <- apply(matvals[,3:ncol(matvals)],2,FUN = function(x){rmultinom(n = 1, size = 1, prob = x)})
  xysim[i,1,] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],1]
  xysim[i,2,] <- matvals[which(multidraw == 1, arr.ind = TRUE)[,1],2]
}
}

if(return == "probability" & Data.Frame == FALSE){return(assign2prob)}
if(return == "probability" & Data.Frame == TRUE){return(assing2probDF)}
if(return == "odds" & Data.Frame == FALSE){return(step2)}
if(return == "odds" & Data.Frame == TRUE){return(step2DF)}
if(return == "population" & Data.Frame == FALSE){return(SamplePop)}
if(return == "population" & Data.Frame == TRUE){return(SamplePopDF)}
if(return == "sim.cell"){return(xysim)}
}
