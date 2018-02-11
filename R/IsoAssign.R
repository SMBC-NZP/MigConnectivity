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
                      return = "probability",
                      element = "Hydrogen",
                      surface = FALSE,
                      period = "Annual"){

isomap <- getIsoMap(element = element, surface = surface, period = "Annual")

animap <- raster::calc(isomap, function(x){y <- slope*x+intercept})

# spatially explicit assignment
assign <- function(x,y) {((1/(sqrt(2 * 3.14 * isoSTD))) * exp((-1/(2 * isoSTD^2)) * ((x) - y)^2))}

assignments <- lapply(isovalues, FUN = function(x){assign(x, y = animap)})

assignments <- raster::stack(assignments)

# Transform the assignments into a true probability surface #
assign2prob <- assignments/raster::cellStats(assignments,sum)

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
}

if(return == "probability"){return(assign2prob)}
if(return == "odds"){return(step2)}
if(return == "population"){return(SamplePop)}

}
