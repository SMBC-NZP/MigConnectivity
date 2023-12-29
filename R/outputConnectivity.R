# @export
is.isoAssign <- function(x) inherits(x, "isoAssign")
# @export
is.estMC <- function(x) inherits(x, "estMC")

# @export
#print <- function(x,...) UseMethod("print")
#' @export
print.estMigConnectivity <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Migratory Connectivity Estimates\n")
  if (inherits(x, "estPsi")) {
    dimnames(x$psi$mean) <- dimnames(x$psi$se) <- list(x$input$originNames,
                                                       x$input$targetNames)
    cat("\nTransition probability (psi) estimates (mean):\n")
    print(x$psi$mean, digits = digits)
    cat("+/- SE:\n")
    print(x$psi$se, digits = digits)
    cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile):\n", sep = "")
    print(array(paste(format(x$psi$simpleCI[1,,],digits = digits, trim = TRUE),
                      format(x$psi$simpleCI[2,,],digits = digits, trim = TRUE),
                      sep = ' - '), dim = dim(x$psi$mean),
                dimnames = list(x$input$originNames, x$input$targetNames)),
          quote = FALSE)
  }
  if (inherits(x, "estMC")) {
    cat("\nMC estimate (mean):", format(x$MC$mean, digits = digits), "+/- (SE)",
        format(x$MC$se, digits = digits), '\n')
    cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile): ",
        paste(format(x$MC$simpleCI, digits = digits, trim = TRUE),
              collapse = ' - '), '\n', sep = "")
  }
  if (inherits(x, "estMantel")) {
    if (is.null(x$corr)) {
      x$corr <- list(mean = x$meanCorr, se = x$seCorr,
                     simpleCI = x$simpleCICorr)
      x$input <- list(alpha = x$alpha)
    }
    cat("\nrM estimate (mean):", format(x$corr$mean, digits = digits), "+/- (SE)",
        format(x$corr$se, digits = digits), '\n')
    cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
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
  if (inherits(x, "estGamma")) {
    cat("\nReverse transition probability (gamma) estimates (mean):\n")
    print(x$gamma$mean, digits = digits)
    cat("+/- SE:\n")
    print(x$gamma$se, digits = digits)
    cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile):\n", sep = "")
    print(array(paste(format(x$gamma$simpleCI[1,,],digits = digits, trim = TRUE),
                      format(x$gamma$simpleCI[2,,],digits = digits, trim = TRUE),
                      sep = ' - '), dim = dim(x$gamma$mean),
                dimnames = list(x$input$targetNames, x$input$originNames)),
          quote = FALSE)
  }
  if (inherits(x, "estTargetRelAbund")) {
    cat("\nTarget site relative abundance estimates (mean):\n")
    print(x$targetRelAbund$mean, digits = digits)
    cat("+/- SE:\n")
    print(x$targetRelAbund$se, digits = digits)
    cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile):\n", sep = "")
    print(array(paste(format(x$targetRelAbund$simpleCI[1,],digits = digits, trim = TRUE),
                      format(x$targetRelAbund$simpleCI[2,],digits = digits, trim = TRUE),
                      sep = ' - '), dim = dim(x$gamma$mean)[1],
                dimnames = list(x$input$targetNames)),
          quote = FALSE)
  }
  if (inherits(x, "estPi")) {
    cat("\nOrigin/target site combination probability (pi) estimates (mean):\n")
    print(x$pi$mean, digits = digits)
    cat("+/- SE:\n")
    print(x$pi$se, digits = digits)
    cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
        "% confidence interval (simple quantile):\n", sep = "")
    print(array(paste(format(x$pi$simpleCI[1,,],digits = digits, trim = TRUE),
                      format(x$pi$simpleCI[2,,],digits = digits, trim = TRUE),
                      sep = ' - '), dim = dim(x$pi$mean),
                dimnames = list(x$input$originNames, x$input$targetNames)),
          quote = FALSE)
  }
  cat("\nThis is a subset of what's available inside this estMigConnectivity output.\n")
  if (inherits(x, "estPi"))
    cat("For more info, try ?reverseTransition or str(obj_name, max.level = 2).\n")
  else if (inherits(x, "estPsi"))
    cat("For more info, try ?estTransition or str(obj_name, max.level = 2).\n")
  else if (inherits(x, "estMC"))
    cat("For more info, try ?estStrength or ?estMC or str(obj_name, max.level = 2).\n")
  else if (inherits(x, "estGamma"))
    cat("For more info, try ?reverseTransition or str(obj_name, max.level = 2).\n")
  else if (inherits(x, "estMantel"))
    cat("For more info, try ?estMantel or str(obj_name, max.level = 2).\n")
  else # In case we left something out...
    cat("For more info, try str(obj_name, max.level = 2).\n")
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
#' if (interactive()) {
#'   OVENdist <- terra::vect("data-raw/Spatial_Layers/OVENdist.shp")
#'   OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
#'   terra::crs(OVENdist) <-
#'     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#'
#'   OVENvals <- read.csv("data-raw/deltaDvalues.csv")
#'
#'   b <- isoAssign(isovalues = OVENvals[,2],
#'                  isoSTD = 12,
#'                  intercept = -10,
#'                  slope = 0.8,
#'                  odds = NULL,
#'                  restrict2Likely = TRUE,
#'                  nSamples = 1000,
#'                  sppShapefile = OVENdist,
#'                  assignExtent = c(-179,-60,15,89),
#'                  element = "Hydrogen",
#'                  surface = FALSE,
#'                  period = "Annual")
#'
#'   plot(b, map = "population")
#' }
#'
#' @export
plot.intrinsicAssign <- function(x,map,...){
  if(inherits(x,"isoAssign")){
    if(!(map %in% c("probability","population","odds"))){
      stop("map must be either probability, population, or odds")}
    op <- graphics::par(no.readonly = TRUE)
    if(map == "population"){
      terra::plot(x$popassign,horiz = TRUE,...)
      }
    if(map == "probability"){
      for(i in 1:terra::nlyr(x$probassign)){
        terra::plot(x$probassign[[i]],horiz = TRUE,...)
        graphics::par(ask = TRUE)
      }
    graphics::par(op)
    }
    if(map == "odds"){
      for(i in 1:terra::nlyr(x$probassign)){
        terra::plot(x$oddsassign[[i]],horiz = TRUE,...)
        graphics::par(ask = TRUE)
      }
    graphics::par(op)
  }
  on.exit(graphics::par(op))
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

#plot <- function(x,...) UseMethod("plot")
#' Basic plot function for estMigConnectivity objects
#'
#' @param x an estMigConnectivity object (output of estTransition, estStrength,
#'   estMC, or estMantel)
#' @param plot.which which parameter (psi, MC, rM, or r) to graph. Defaults to
#'   psi for estMC objects, to rM (Mantel correlation) otherwise
#' @param point points on graph can represent mean, median, or point estimates
#'   (not considering error). Defaults to mean, the standard estimate from
#'   resampling
#' @param range lines / error bars drawn around points can represent simple
#'   quantile-based confidence intervals (simpleCI), bias-corrected quantile-
#'   based confidence intervals (bcCI), or +- standard error (se). Defaults to
#'   simpleCI
#' @param xlab label for the x-axis. Defaults to "Origin" for psi, otherwise ""
#' @param ylab label for the y-axis. Defaults to the parameter being plotted
#' @param originNames names of the origin sites (for plotting psi). If left
#'   NULL, the function attempts to get these from the estimate
#' @param targetNames names of the target sites (for plotting psi or r). If left
#'   NULL, the function attempts to get these from the estimate
#' @param ageNames names of the age classes (for plotting r with more than one
#'   age). If left NULL, the function uses 1:[number of ages]
#' @param col colors to use for labeling transition probabilities for
#'   different target sites. If left NULL, defaults to 1:[number of target sites]
#' @param pch symbols to use for labeling transition probabilities for
#'   different target sites. If left NULL, defaults to 21:25, then
#'   0:([number of target sites]-5)
#' @param las style of axis labels (0-3). We set the default at 1 (always
#'   horizontal) here, but if you prefer your labels parallel to the axis, set
#'   at 0
#' @param gap space left between the center of the error bar and the lines
#'   marking the error bar in units of the height (width) of the letter "O".
#'   Defaults to 0
#' @param sfrac width of "crossbar" at the end of error bar as a fraction of the
#'   x plotting region. Defaults to 0, unless range is set to "se", in which
#'   case it defaults to 0.01
#' @param legend leave as FALSE to not print a legend. Otherwise the position
#'   of the legend (for psi or r (multi-age) only; one of "bottomright",
#'   "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", or
#'   "center")
#' @param map placeholder for eventually allowing users to plot psi estimates
#'   on a map
#' @param ... Additional parameters passed to \code{\link{plotCI}}
#'
#' @seealso \code{\link{estMC}}, \code{\link{estMantel}}
#'
#' @export
plot.estMigConnectivity <- function(x,
                                    plot.which = ifelse(inherits(x, "estPsi"),
                                                        "psi",
                                                        ifelse(inherits(x, "estMC"),
                                                               "MC",
                                                               ifelse(inherits(x, "estGamma"),
                                                                      "gamma",
                                                                      "rM"))),
                                    point = c("mean", "median", "point"),
                                    range = c("simpleCI", "bcCI", "se"),
                                    xlab = NULL, ylab = plot.which,
                                    originNames = NULL, targetNames = NULL,
                                    ageNames = NULL,
                                    col = NULL, pch = NULL, las = 1,
                                    gap = 0,
                                    sfrac = ifelse(range=="se", 0.01, 0),
                                    legend = FALSE, map = FALSE, ...) {
  if ((plot.which %in% c("corr", "Mantel")))
    plot.which <- "rM"
  if ((plot.which %in% c("transition", "Transition")))
    plot.which <- "psi"
  if ((plot.which %in% c("strength", "Strength")))
    plot.which <- "MC"
  if ((plot.which %in% c("abund", "abundance", "relAbund")))
    plot.which <- "targetRelAbund"
  if (!(plot.which %in% c("psi", "MC", "rM", "r", "gamma", "pi",
                          "targetRelAbund")))
    stop("Set plot.which to psi, MC, rM, gamma, pi, targetRelAbund, or r")
  point <- match.arg(point)
  range <- match.arg(range)
  if (inherits(x, "estMC")) {
    if (is.null(x$psi) && !is.null(x$samplePsi)) {
      bcCIPsi <- array(NA, dim = c(2, dim(x$samplePsi)[2], dim(x$samplePsi)[3]))
      for (i in 1:dim(x$samplePsi)[2]) {
        for (j in 1:dim(x$samplePsi)[3]) {
          psi.z0 <- qnorm(sum(x$samplePsi[, i, j] < mean(x$samplePsi[, i, j],
                                                         na.rm = TRUE)) /
                            length(which(!is.na(x$samplePsi[, i, j]))))
          bcCIPsi[ , i, j] <- quantile(x$samplePsi[, i, j],
                                       pnorm(2 * psi.z0 +
                                               qnorm(c(x$alpha/2,1-x$alpha/2))),
                                       na.rm=TRUE, names = FALSE)
        }
      }
      x$psi <- list(sample = x$samplePsi,
                    mean = apply(x$samplePsi, 2:3, mean),
                    se = apply(x$samplePsi, 2:3, sd),
                    simpleCI = apply(x$samplePsi, 2:3, quantile,
                                     probs = c(x$alpha/2, 1-x$alpha/2),
                                     na.rm=TRUE, names = FALSE),
                    bcCI = bcCIPsi,
                    median = apply(x$samplePsi, 2:3, median),
                    point = x$pointPsi)
      x$MC <- list(mean = x$meanMC, se = x$seMC, simpleCI = x$simpleCI,
                   bcCI = x$bcCI, median = x$medianMC, point = x$pointMC)
      x$input <- list(originNames = dimnames(x$samplePsi)[[2]],
                      targetNames = dimnames(x$samplePsi)[[3]])
    }
  }
  else if (plot.which %in% c("psi", "MC") && !inherits(x, "estPsi"))
    stop("This estimate does not include psi or MC - try setting plot.which to rM")
  if (plot.which == "gamma" && !inherits(x, "estGamma"))
    stop("This estimate does not include gamma - try setting plot.which to something else")
  if (plot.which == "pi" && !inherits(x, "estPi"))
    stop("This estimate does not include pi - try setting plot.which to something else")
  if (plot.which == "targetRelAbund" && !inherits(x, "estTargetRelAbund"))
    stop("This estimate does not include targetRelAbund - try setting plot.which to something else")
  if (inherits(x, "estMantel")) {
    if (is.null(x$corr)) {
      x$corr <- list(sample = x$sampleCorr,
                     mean = x$meanCorr,
                     se = x$seCorr,
                     simpleCI = x$simpleCICorr,
                     bcCI = x$bcCICorr,
                     median = x$medianCorr,
                     point = x$pointCorr)
      x$input <- list(alpha = x$alpha)
    }
  }
  else if (plot.which == "rM")
    stop("This estimate does not include rM - try setting plot.which to MC or psi")
  if (plot.which=="MC") {
    y <- ifelse(point == "mean", x$MC$mean,
                ifelse(point == "median", x$MC$median, x$MC$point))
    ests.df <- data.frame(y = y,
                          lower = ifelse(range == "bcCI", x$MC$bcCI[1],
                                         ifelse(range == "simpleCI",
                                                x$MC$simpleCI[1],
                                                y - x$MC$se)),
                          upper = ifelse(range == "bcCI", x$MC$bcCI[2],
                                         ifelse(range == "simpleCI",
                                                x$MC$simpleCI[2],
                                                y + x$MC$se)))
  }
  else if (plot.which == "rM") {
    y <- ifelse(point == "mean", x$corr$mean,
                ifelse(point == "median", x$corr$median, x$corr$point))
    ests.df <- data.frame(y = y,
                          lower = ifelse(range == "bcCI", x$corr$bcCI[1],
                                         ifelse(range == "simpleCI",
                                                x$corr$simpleCI[1],
                                                y - x$corr$se)),
                          upper = ifelse(range == "bcCI", x$corr$bcCI[2],
                                         ifelse(range == "simpleCI",
                                                x$corr$simpleCI[2],
                                                y + x$corr$se)))
  }
  else if (plot.which %in% c("psi", "gamma", "pi")) {
    if (is.null(originNames)) {
      if (is.null(x$input$originNames)) {
        if (is.null(dimnames(x$psi$sample)[2])) {
          originNames <- LETTERS[1:dim(x$psi$sample)[2]]
        }
        else {
          originNames <- dimnames(x$psi$sample)[[2]]
        }
      }
      else
        originNames <- x$input$originNames
    }
    if (is.null(targetNames)) {
      if (is.null(x$input$targetNames)) {
        if (is.null(dimnames(x$psi$sample)[3])) {
          targetNames <- 1:dim(x$psi$sample)[3]
        }
        else {
          targetNames <- dimnames(x$psi$sample)[[3]]
        }
      }
      else
        targetNames <- x$input$targetNames
    }
    nTargetSites <- ifelse(is.null(x$psi), dim(x$gamma$sample)[2], dim(x$psi$sample)[3])
    nOriginSites <- ifelse(is.null(x$psi), dim(x$gamma$sample)[3], dim(x$psi$sample)[2])
    if (length(originNames)==1) {
      originNames <- rep(originNames, nOriginSites)
      originNamesHidden <- 1:nOriginSites
    }
    else {
      originNamesHidden <- originNames
    }
    if (length(targetNames)==1) {
      targetNames <- rep(targetNames, nTargetSites)
      targetNamesHidden <- 1:nTargetSites
    }
    else {
      targetNamesHidden <- targetNames
    }
    if (plot.which == "psi") {
      y <- switch(point,
                  mean = x$psi$mean,
                  median = x$psi$median,
                  point = x$psi$point)
      yrange <- switch(range,
                       bcCI = x$psi$bcCI,
                       simpleCI = x$psi$simpleCI,
                       se = aperm(array(c(y - x$psi$se, y + x$psi$se),
                                        c(dim(y), 2)), c(3, 1, 2)))
      ests.df <- data.frame(y = c(y),
                            To = factor(rep(targetNamesHidden, each = nOriginSites),
                                        levels = targetNamesHidden),
                            FromTo = rep(1:nOriginSites, nTargetSites) +
                              rep(1:nTargetSites, each = nOriginSites) /
                              (nTargetSites * 2) - 0.3,
                            lower = c(yrange[1,,]),
                            upper = c(yrange[2,,]))
    }
    else if (plot.which == "pi") {
      y <- switch(point,
                  mean = x$pi$mean,
                  median = x$pi$median,
                  point = x$pi$point)
      yrange <- switch(range,
                       bcCI = x$pi$bcCI,
                       simpleCI = x$pi$simpleCI,
                       se = aperm(array(c(y - x$pi$se, y + x$pi$se),
                                        c(dim(y), 2)), c(3, 1, 2)))
      ests.df <- data.frame(y = c(y),
                            To = factor(rep(targetNamesHidden, each = nOriginSites),
                                        levels = targetNamesHidden),
                            FromTo = rep(1:nOriginSites, nTargetSites) +
                              rep(1:nTargetSites, each = nOriginSites) /
                              (nTargetSites * 2) - 0.3,
                            lower = c(yrange[1,,]),
                            upper = c(yrange[2,,]))
    }
    else {
      y <- switch(point,
                  mean = x$gamma$mean,
                  median = x$gamma$median,
                  point = x$gamma$point)
      yrange <- switch(range,
                       bcCI = x$gamma$bcCI,
                       simpleCI = x$gamma$simpleCI,
                       se = aperm(array(c(y - x$gamma$se, y + x$gamma$se),
                                        c(dim(y), 2)), c(3, 1, 2)))
      ests.df <- data.frame(y = c(y),
                            To = factor(rep(originNamesHidden, each = nTargetSites),
                                        levels = originNamesHidden),
                            FromTo = rep(1:nTargetSites, nOriginSites) +
                              rep(1:nOriginSites, each = nTargetSites) /
                              (nOriginSites * 2) - 0.3,
                            lower = c(yrange[1,,]),
                            upper = c(yrange[2,,]))

    }
  }
  else if (plot.which == "r" || plot.which == "targetRelAbund") {
    if (is.null(targetNames)) {
      if (is.null(x$input$targetNames)) {
        if (is.null(dimnames(x$psi$sample)[3])) {
          targetNames <- 1:dim(x$psi$sample)[3]
        }
        else {
          targetNames <- dimnames(x$psi$sample)[[3]]
        }
      }
      else
        targetNames <- x$input$targetNames
    }
    nTargetSites <- length(targetNames)
    if (plot.which == "r") {
      y <- switch(point,
                  mean = x$r$mean,
                  median = x$r$median)
      if (length(dim(y)) < 2){
        nAges <- 1
        yrange <- switch(range,
                         bcCI = x$r$bcCI,
                         simpleCI = x$r$simpleCI,
                         se = aperm(array(c(y - x$r$se, y + x$r$se),
                                          c(length(y), 2)), c(2, 1)))
        ests.df <- data.frame(y = c(y),
                              To = 1:nTargetSites,
                              lower = c(yrange[1,]),
                              upper = c(yrange[2,]))
      }
      else{
        nAges <- dim(y)[1]
        if (is.null(ageNames))
          ageNames <- 1:nAges
        yrange <- switch(range,
                         bcCI = x$r$bcCI,
                         simpleCI = x$r$simpleCI,
                         se = aperm(array(c(y - x$r$se, y + x$r$se),
                                          c(dim(y), 2)), c(3, 1, 2)))
        ests.df <- data.frame(y = c(y),
                              Age = rep(1:nAges, nTargetSites),
                              To = rep(1:nTargetSites, each = nAges),
                              AgeTo = rep(1:nAges, nTargetSites) +
                                rep(1:nTargetSites, each = nAges) /
                                (nTargetSites * 2) - 0.3,
                              lower = c(yrange[1,,]),
                              upper = c(yrange[2,,]))
      }
    }
    else { # plotting targetRelAbund
      y <- switch(point,
                  mean = x$targetRelAbund$mean,
                  median = x$targetRelAbund$median)
      yrange <- switch(range,
                       bcCI = x$targetRelAbund$bcCI,
                       simpleCI = x$targetRelAbund$simpleCI,
                       se = aperm(array(c(y - x$targetRelAbund$se,
                                          y + x$targetRelAbund$se),
                                        c(length(y), 2)), c(2, 1)))
      ests.df <- data.frame(y = c(y),
                            To = 1:nTargetSites,
                            lower = c(yrange[1,]),
                            upper = c(yrange[2,]))
      nAges <- 1

    }
  }
  if (map) {
    warning("Map plotting not yet available")
  }
  if (plot.which %in% c("psi", "pi")) {
    if (is.null(col)) {
      col <- 1:nTargetSites
    }
    else if (length(col) < nTargetSites)
      col <- rep_len(col, nTargetSites)
    if (is.null(pch)) {
      if (nTargetSites > 5)
        pch <- c(21:25, 0:(nTargetSites - 6))
      else
        pch <- 20 + 1:nTargetSites
    }
    else if (length(pch) < nTargetSites)
      pch <- rep_len(pch, nTargetSites)
    if (is.null(xlab))
      xlab <- "Origin"
    gplots::plotCI(ests.df$FromTo, ests.df$y, li = ests.df$lower,
                   ui = ests.df$upper,
                   pch = pch[as.integer(ests.df$To)],
                   col = col[as.integer(ests.df$To)],
                   pt.bg = col[as.integer(ests.df$To)],
                   ylab = ylab, xaxt = "n", xlab = xlab, gap = gap,
                   sfrac = sfrac, las = las, ...)
    graphics::axis(1, at = seq(from = 1, by = 1, length.out = nOriginSites),
         labels = originNames, ...)
    if (!isFALSE(legend))
      legend(legend, legend = targetNames, col = col, pch = pch,
             pt.bg = col)
  }
  else if (plot.which=="gamma") {
    if (is.null(col)) {
      col <- 1:nOriginSites
    }
    else if (length(col) < nOriginSites)
      col <- rep_len(col, nOriginSites)
    if (is.null(pch)) {
      if (nOriginSites > 5)
        pch <- c(21:25, 0:(nOriginSites - 6))
      else
        pch <- 20 + 1:nOriginSites
    }
    else if (length(pch) < nOriginSites)
      pch <- rep_len(pch, nOriginSites)
    if (is.null(xlab))
      xlab <- "Target"
    gplots::plotCI(ests.df$FromTo, ests.df$y, li = ests.df$lower,
                   ui = ests.df$upper,
                   pch = pch[as.integer(ests.df$To)],
                   col = col[as.integer(ests.df$To)],
                   pt.bg = col[as.integer(ests.df$To)],
                   ylab = ylab, xaxt = "n", xlab = xlab, gap = gap,
                   sfrac = sfrac, las = las, ...)
    graphics::axis(1, at = seq(from = 1, by = 1, length.out = nTargetSites),
                   labels = targetNames, ...)
    if (!isFALSE(legend))
      legend(legend, legend = originNames, col = col, pch = pch,
             pt.bg = col)
  }
  else if (plot.which=="r" || plot.which=="targetRelAbund") {
    if (nAges == 1) {
      if (is.null(col)) {
        col <- "black"
      }
      if (is.null(pch)) {
        pch <- 20
      }
      if (is.null(xlab))
        xlab <- "Target"
      gplots::plotCI(ests.df$To, ests.df$y, li = ests.df$lower,
                     ui = ests.df$upper,
                     pch = pch,
                     col = col,
                     pt.bg = col,
                     ylab = ylab, xlab = xlab, gap = gap, xaxt = "n",
                     sfrac = sfrac, las = las, ...)
      graphics::axis(1, at = seq(from = 1, by = 1, length.out = nTargetSites),
           labels = targetNames, ...)
    }
    else {
      if (is.null(col)) {
        col <- 1:nTargetSites
      }
      else if (length(col) < nTargetSites)
        col <- rep_len(col, nTargetSites)
      if (is.null(pch)) {
        if (nTargetSites > 5)
          pch <- c(21:25, 0:(nTargetSites - 6))
        else
          pch <- 20 + 1:nTargetSites
      }
      else if (length(pch) < nTargetSites)
        pch <- rep_len(pch, nTargetSites)
      if (is.null(xlab))
        xlab <- "Age"
      gplots::plotCI(ests.df$AgeTo, ests.df$y, li = ests.df$lower,
                     ui = ests.df$upper,
                     pch = pch[as.integer(ests.df$To)],
                     col = col[as.integer(ests.df$To)],
                     pt.bg = col[as.integer(ests.df$To)],
                     ylab = ylab, xaxt = "n", xlab = xlab, gap = gap,
                     sfrac = sfrac, las = las, ...)
      graphics::axis(1, at = 1:nAges, labels = ageNames, ...)
      if (!isFALSE(legend))
        legend(legend, legend = targetNames, col = col, pch = pch,
               pt.bg = col)
    }
  }
  else {
    if (is.null(col)) {
      col <- "black"
    }
    if (is.null(pch)) {
      pch <- 20
    }
    if (is.null(xlab))
      xlab <- ""
    gplots::plotCI(0, ests.df$y, li = ests.df$lower,
                   ui = ests.df$upper,
                   pch = pch[1],
                   col = col[1], gap = gap, sfrac = sfrac,
                   ylab = ylab, xaxt = "n", xlab = xlab, las = las, ...)
  }
}
map.estPsi <- function(x, originSites, targetSites, xOffset = NULL,
                       yOffset = NULL, col = NULL, maxWidth = 100000,
                       alpha.range = 0.2, alpha.point = 0,
                       subsetOrigin = NULL, subsetTarget = NULL,
                       doubled = FALSE) {
# #
# # data(OVENdata) # Ovenbird
# #
# M<-estMC(isGL=OVENdata$isGL, # Logical vector: light-level geolocator(T)/GPS(F)
#          geoBias = OVENdata$geo.bias, # Light-level geolocator location bias
#          geoVCov = OVENdata$geo.vcov, #Light-level geolocator covariance matrix
#          targetDist = OVENdata$targetDist, # Target location distance matrix
#          originDist = OVENdata$originDist, # Origin location distance matrix
#          targetSites = OVENdata$targetSites, # Non-breeding / target sites
#          originSites = OVENdata$originSites, # Breeding / origin sites
#          originPoints = OVENdata$originPoints, # Capture Locations
#          targetPoints = OVENdata$targetPoints, # Target locations from devices
#          originRelAbund = OVENdata$originRelAbund, # Origin relative abundances
#          resampleProjection = terra::crs(OVENdata$targetPoints),
#          verbose = 0,   # output options - see help ??estMC
#          nSamples = 10000) # This is set low for example
#
  nTargetSites <- ncol(x$psi$mean)
  nOriginSites <- nrow(x$psi$mean)
  if (is.null(subsetOrigin))
    subsetOrigin <- 1:nOriginSites
  if (is.null(subsetTarget))
    subsetTarget <- 1:nTargetSites
  originNames <- x$input$originNames
  targetNames <- x$input$targetNames
  if (is.null(xOffset))
    xOffset <- matrix(0, nOriginSites, nTargetSites)
  if (is.null(yOffset))
    yOffset <- matrix(0, nOriginSites, nTargetSites)
  allSites <- rbind(originSites[, "geometry"], targetSites[, "geometry"])
  if (is.null(col)) {
    col <- 1:nTargetSites
  }
  # png('psi_plot1.png', width = 6, height = 6, units = 'in', res = 1200)
  op <- graphics::par(mar=c(0,0,0,0))
  extent <- sf::st_bbox(allSites)
  plot(allSites, xlim=c(extent[1],extent[3]),
       ylim=c(extent[2], extent[4]), lwd = 1.5)
  # plot(OVENdata$originSites[1],add=TRUE,lwd=1.75)
  # plot(OVENdata$originSites[2],add=TRUE,lwd=1.75)
  # plot(OVENdata$targetSites,add=TRUE,lwd=1.5,col=c("gray70","gray35","gray10"))

  # legend("topleft",legend=paste("MC =",round(M$meanMC,2), "\u00b1", round(M$seMC,2)),bty="n",cex=1.8,bg="white",xjust=0)
  for (i in subsetOrigin) {
    xO <- sf::st_coordinates(sf::st_centroid(originSites[i,]))[,1]
    yO <- sf::st_coordinates(sf::st_centroid(originSites[i,]))[,2]
    for (j in subsetTarget) {
      if (x$psi$simpleCI[2,i,j] > 0) {
        xT <- sf::st_coordinates(sf::st_centroid(targetSites[j,]))[,1] + xOffset[i, j]
        yT <- sf::st_coordinates(sf::st_centroid(targetSites[j,]))[,2] + yOffset[i, j]
        angle <- atan((yT - yO)/(xT - xO))
        if (is.nan(angle))
          angle <- 0
        if (xT < xO)
          angle <- angle + pi
        cosa <- cos(angle)
        sina <- sin(angle)
        if (doubled) {
          graphics::polygon(c(xO + x$psi$simpleCI[2,i,j] * sina * maxWidth,
                    xT + x$psi$simpleCI[2,i,j] * sina * maxWidth,
                    xT + x$psi$simpleCI[1,i,j] * sina * maxWidth,
                    xO + x$psi$simpleCI[1,i,j] * sina * maxWidth),
                  c(yO - x$psi$simpleCI[2,i,j] * cosa * maxWidth,
                    yT - x$psi$simpleCI[2,i,j] * cosa * maxWidth,
                    yT - x$psi$simpleCI[1,i,j] * cosa * maxWidth,
                    yO - x$psi$simpleCI[1,i,j] * cosa * maxWidth),
                  col = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                            1 - (i + j) / (nOriginSites + nTargetSites),
                            alpha=alpha.range), border = NA)
          graphics::polygon(c(xO - x$psi$simpleCI[2,i,j] * sina * maxWidth,
                    xT - x$psi$simpleCI[2,i,j] * sina * maxWidth,
                    xT - x$psi$simpleCI[1,i,j] * sina * maxWidth,
                    xO - x$psi$simpleCI[1,i,j] * sina * maxWidth),
                  c(yO + x$psi$simpleCI[2,i,j] * cosa * maxWidth,
                    yT + x$psi$simpleCI[2,i,j] * cosa * maxWidth,
                    yT + x$psi$simpleCI[1,i,j] * cosa * maxWidth,
                    yO + x$psi$simpleCI[1,i,j] * cosa * maxWidth),
                  col = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                            1 - (i + j) / (nOriginSites + nTargetSites),
                            alpha=alpha.range), border = NA)
          graphics::polygon(c(xO - x$psi$mean[i,j] * sina * maxWidth,
                    xT - x$psi$mean[i,j] * sina * maxWidth,
                    xT + x$psi$mean[i,j] * sina * maxWidth,
                    xO + x$psi$mean[i,j] * sina * maxWidth),
                  c(yO + x$psi$mean[i,j] * cosa * maxWidth,
                    yT + x$psi$mean[i,j] * cosa * maxWidth,
                    yT - x$psi$mean[i,j] * cosa * maxWidth,
                    yO - x$psi$mean[i,j] * cosa * maxWidth),
                  col = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                            1 - (i + j) / (nOriginSites + nTargetSites),
                            alpha=alpha.point),
                  border = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                               1 - (i + j) / (nOriginSites + nTargetSites),
                               alpha=1))
        }
        else {
          graphics::polygon(c(xO - x$psi$mean[i,j] * sina * maxWidth / 2,
                    xT - x$psi$mean[i,j] * sina * maxWidth / 2,
                    xT + x$psi$mean[i,j] * sina * maxWidth / 2,
                    xO + x$psi$mean[i,j] * sina * maxWidth / 2),
                  c(yO + x$psi$mean[i,j] * cosa * maxWidth / 2,
                    yT + x$psi$mean[i,j] * cosa * maxWidth / 2,
                    yT - x$psi$mean[i,j] * cosa * maxWidth / 2,
                    yO - x$psi$mean[i,j] * cosa * maxWidth / 2),
                  col = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                            1 - (i + j) / (nOriginSites + nTargetSites),
                            alpha=alpha.point),
                  border = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                               1 - (i + j) / (nOriginSites + nTargetSites),
                               alpha=1))
          graphics::polygon(c(xO - sina * maxWidth * (x$psi$simpleCI[2,i,j] -
                                              x$psi$mean[i,j] / 2),
                    xT - (x$psi$simpleCI[2,i,j] - x$psi$mean[i,j] / 2) *
                      sina * maxWidth,
                    xT - (x$psi$simpleCI[1,i,j] - x$psi$mean[i,j] / 2) *
                      sina * maxWidth,
                    xO - (x$psi$simpleCI[1,i,j] - x$psi$mean[i,j] / 2) *
                      sina * maxWidth),
                  c(yO+(x$psi$simpleCI[2,i,j]-x$psi$mean[i,j]/2)*cosa*maxWidth,
                    yT+(x$psi$simpleCI[2,i,j]-x$psi$mean[i,j]/2)*cosa*maxWidth,
                    yT+(x$psi$simpleCI[1,i,j]-x$psi$mean[i,j]/2)*cosa*maxWidth,
                    yO+(x$psi$simpleCI[1,i,j]-x$psi$mean[i,j]/2)*cosa*maxWidth),
                  col = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                            1 - (i + j) / (nOriginSites + nTargetSites),
                            alpha=alpha.range), border = NA)

        }
        shape::Arrowhead(xT, yT, angle / pi * 180, arr.width = x$psi$mean[i,j], arr.length = 1/8,
                         arr.type = 'curved', npoint = 15,
                         lcol = grDevices::rgb(i/nOriginSites, j/nTargetSites,
                                    1 - (i + j) / (nOriginSites + nTargetSites),
                                    alpha=1), arr.adj = 0)
      }
    }
  }
# dev.off()
# shape::Arrows(gCentroid(OVENdata$originSites[1])@coords[,1],
#               gCentroid(OVENdata$originSites[1])@coords[,2],
#               gCentroid(OVENdata$targetSites[2])@coords[,1]+80000,
#               extent(OVENdata$targetSites[2])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,1,],2,mean)[2]*10),
#               lty=1)
#
# shape::Arrows(gCentroid(OVENdata$originSites[1])@coords[,1],
#               gCentroid(OVENdata$originSites[1])@coords[,2],
#               gCentroid(OVENdata$targetSites[3])@coords[,1],
#               extent(OVENdata$targetSites[3])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,1,],2,mean)[3]*10),
#               lty=1)
#
# shape::Arrows(gCentroid(OVENdata$originSites[2])@coords[,1],
#               gCentroid(OVENdata$originSites[2])@coords[,2],
#               gCentroid(OVENdata$targetSites[1])@coords[,1],
#               extent(OVENdata$targetSites[1])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,2,],2,mean)[1]*10),
#               lty=1)
#
# shape::Arrows(gCentroid(OVENdata$originSites[2])@coords[,1],
#               gCentroid(OVENdata$originSites[2])@coords[,2],
#               (gCentroid(OVENdata$targetSites[2])@coords[,1]-80000),
#               extent(OVENdata$targetSites[2])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,2,],2,mean)[2]*10))
#
# box(which="plot")
# #
  graphics::par(op)
}
