# @export
is.isoAssign <- function(x) inherits(x, "isoAssign")
# @export
is.estMC <- function(x) inherits(x, "estMC")

# @export
#print <- function(x,...) UseMethod("print")
#' @export
print.estMigConnectivity <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Migratory Connectivity Estimate(s)\n")
  if (inherits(x, "estMC")) {
    if (is.null(x$psi)) {
      x$psi <- list(mean = apply(x$samplePsi, 2:3, mean),
                    se = apply(x$samplePsi, 2:3, sd),
                    simpleCI = apply(x$samplePsi, 2:3, quantile,
                                     probs = c(alpha/2, 1-alpha/2),
                                     na.rm=TRUE, type = 8, names = F))
      x$MC <- list(mean = x$meanMC, se = x$seMC, simpleCI = x$simpleCI)
      x$input <- list(alpha = x$alpha,
                      originNames = ifelse(is.null(dimnames(x$samplePsi)[2]),
                                           LETTERS[1:dim(x$samplePsi)[2]],
                                           dimnames(x$samplePsi)[2]),
                      targetNames = ifelse(is.null(dimnames(x$samplePsi)[3]),
                                           1:dim(x$samplePsi)[3],
                                           dimnames(x$samplePsi)[3]))
    }
    cat("Transition probability (psi) estimates (mean):\n")
    print(x$psi$mean)
    cat("   se:\n")
    print(x$psi$se)
    cat("   ", ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile):\n")
    print(array(paste(format(x$psi$simpleCI[1,,],digits = digits, trim = TRUE),
                      format(x$psi$simpleCI[2,,],digits = digits, trim = TRUE),
                      sep = ' - '), dim = dim(x$psi$mean),
                dimnames = list(x$input$originNames, x$input$targetNames)),
          quote = FALSE)
    cat("MC estimate (mean):", format(x$MC$mean, digits = digits), "+/- (SE)",
        format(x$MC$se, digits = digits), '\n')
    cat("   ", ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile): ",
        paste(format(x$MC$simpleCI, digits = digits, trim = TRUE),
              collapse = ' - '), '\n', sep = "")
    # cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
    #     "% confidence interval (bias-corrected): ",
    #     paste(format(x$bcCI, digits = digits, trim = TRUE), collapse = ' - '),
    #     '\n', sep = "")
    # cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
    #     "% credible interval (highest posterior density): ",
    #     paste(format(x$hpdCI, digits = digits, trim = TRUE), collapse = ' - '),
    #     '\n', sep = "")
    # cat("   median:", format(x$medianMC, digits = digits), '\n')
    # if (!is.na(x$pointMC))
    #   cat("   point calculation (not considering error):",
    #       format(x$pointMC, digits = digits), '\n')
  }
  if (inherits(x, "estMantel")) {
    if (is.null(x$corr)) {
      x$corr <- list(mean = x$meanCorr, se = x$seCorr,
                     simpleCI = x$simpleCICorr)
      x$input <- list(alpha = x$alpha)
    }
    cat("rM estimate (mean):", format(x$corr$mean, digits = digits), "+/- (SE)",
        format(x$corr$se, digits = digits), '\n')
    cat("   ", ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile): ",
        paste(format(x$corr$simpleCI, digits = digits, trim = TRUE),
              collapse = ' - '), '\n', sep = "")
    # cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
    #     "% confidence interval (bias-corrected): ",
    #     paste(format(x$bcCICorr, digits = digits, trim = TRUE), collapse = ' - '),
    #     '\n', sep = "")
    # cat("   median:", format(x$medianCorr, digits = digits), '\n')
    # if (!is.na(x$pointCorr))
    #   cat("   point calculation (not considering error):",
    #       format(x$pointCorr, digits = digits), '\n')
  }
  cat("This is a subset of what's available inside MigConnectivity outputs.\n")
  cat("For more info, try ?estMC or str(obj_name, max.levels = 2).\n")
}

#' @export
print.intrinsicAssign <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
if(inherits(x,"isoAssign")){
  cat("Individual Probability Surfaces \n")
  print(x$probassign,...)
  cat("\n Individual likely/unlikely Surfaces \n")
  print(x$oddsassign,...)
  cat("\n Population-level assignment Surface \n")
  print(x$popassign,...)
  cat("\n Individual Probability data frame* \n")
  utils::str(x$probDF[,1:5])
  cat("\n Individual likely/unlikely data frame* \n")
  utils::str(x$oddsDF[,1:5])
  cat("\n Individual single cell assignment \n")
  utils::str(x$SingleCell)
  cat("\n Target sites spatial layer \n")
  print(x$targetSites)
  cat("\n Random number seed set to: \n")
  print(x$RandomSeed)
  cat("\n * only first few columns are printed")
}
}
# @export
#summary <- function(x,...) UseMethod("summary")
#' @export
summary.estMigConnectivity <- function(object, ...)
{
  print.estMigConnectivity(object, ...)
}

#' @export
summary.intrinsicAssign<-function(object, ...){
  print.intrinsicAssign(object,...)
}

#' Basic plot function for the different isoAssign outputs
#' @param x an isoAssign object
#' @param map which \code{isoAssign} output to plot either 'probability', 'population' or 'odds'
#' @param ... additional arguments passed to plot function
#' @seealso{\code{isoAssign}}
#' @return A basic plot of the isotope assignments. If \code{map = 'population'} returns a single map.
#' If \code{map = 'probability' or map = 'odds'} a map for each individual is returned. User is asked for input before each individual is drawn.
#' @examples
#' \dontrun{
#' OVENdist <- raster::shapefile("data-raw/Spatial_Layers/OVENdist.shp")
#' OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#' raster::crs(OVENdist) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#'
#' OVENvals <- read.csv("data-raw/deltaDvalues.csv")
#'
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
#'
#'plot(b, map = "population")
#'}
#'
#' @export
plot.intrinsicAssign <- function(x,map,...){
  if(inherits(x,"isoAssign")){
plot.isoAssign <- function(x,map,...){
    if(!(map %in% c("probability","population","odds"))){
      stop("map must be either probability, population, or odds")}
    op <- graphics::par(no.readonly = TRUE)
    if(map == "population"){
      raster::plot(x$popassign,horiz = TRUE,...)
      }
    if(map == "probability"){
      for(i in 1:raster::nlayers(x$probassign)){
        raster::plot(x$probassign[[i]],horiz = TRUE,...)
        graphics::par(ask = TRUE)
      }
    graphics::par(op)
    }
    if(map == "odds"){
      for(i in 1:raster::nlayers(x$probassign)){
        raster::plot(x$oddsassign[[i]],horiz = TRUE,...)
        graphics::par(ask = TRUE)
      }
    graphics::par(op)
  }
  on.exit(graphics::par(op))
  }
  }
 if(inherits(x,"weightAssign")){
   graphics::par(bty = "L")
   graphics::plot((x$performance$area/max(x$performance$area))~x$performance$error,
        las = 1, ylab = "Assignment Area",
        xlab = "Error",pch = 19, cex = 1.25, col = "gray")
   graphics::points((x$frontier$area/max(x$performance$area))~x$frontier$error, col = "red", pch = 19, cex = 1.25)
   graphics::points((x$top$area/max(x$performance$area))~x$top$error, col = "blue", pch = 19, cex = 1.25)
 }

}

