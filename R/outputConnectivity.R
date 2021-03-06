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
    # if (is.null(x$psi) && !is.null(x$samplePsi)) {
    #   x$psi <- list(mean = apply(x$samplePsi, 2:3, mean),
    #                 se = apply(x$samplePsi, 2:3, sd),
    #                 simpleCI = apply(x$samplePsi, 2:3, quantile,
    #                                  probs = c(alpha/2, 1-alpha/2),
    #                                  na.rm=TRUE, type = 8, names = F))
    #   x$MC <- list(mean = x$meanMC, se = x$seMC, simpleCI = x$simpleCI)
    #   x$input <- list(alpha = x$alpha)
    #   if (is.null(dimnames(x$samplePsi)[2])){
    #     x$input$originNames <- LETTERS[1:dim(x$samplePsi)[2]]
    #     originNamesFilled <- TRUE
    #   }
    #   else {
    #     x$input$originNames <- dimnames(x$samplePsi)[[2]]
    #     originNamesFilled <- FALSE
    #   }
    #   if (is.null(dimnames(x$samplePsi)[3])){
    #     x$input$targetNames <- 1:dim(x$samplePsi)[3]
    #     targetNamesFilled <- TRUE
    #   }
    #   else {
    #     x$input$targetNames <- dimnames(x$samplePsi)[[3]]
    #     targetNamesFilled <- FALSE
    #   }
    # }
    # else {
    #   originNamesFilled <- targetNamesFilled <- FALSE
    # }
    # if (!is.null(x$psi)) {
    #   dimnames(x$psi$mean) <- dimnames(x$psi$se) <- list(x$input$originNames,
    #                                                      x$input$targetNames)
    #   cat("\nTransition probability (psi) estimates (mean):",
    #       ifelse(originNamesFilled, "(Arbitrary origin site labels used)", ""),
    #       ifelse(targetNamesFilled, "(Arbitrary target site labels used)", ""),
    #       "\n")
    #   print(x$psi$mean, digits = digits)
    #   cat("+/- SE:\n")
    #   print(x$psi$se, digits = digits)
    #   cat(ifelse(is.null(x$input$alpha), "", 100 * (1 - x$input$alpha)),
    #       "% confidence interval (simple quantile):\n", sep = "")
    #   print(array(paste(format(x$psi$simpleCI[1,,],digits = digits, trim = TRUE),
    #                     format(x$psi$simpleCI[2,,],digits = digits, trim = TRUE),
    #                     sep = ' - '), dim = dim(x$psi$mean),
    #               dimnames = list(x$input$originNames, x$input$targetNames)),
    #         quote = FALSE)
    # }
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
  cat("\nThis is a subset of what's available inside this estMigConnectivity output.\n")
  if (inherits(x, "estPsi"))
    cat("For more info, try ?estTransition or ?estMC or str(obj_name, max.levels = 2).\n")
  else if (inherits(x, "estMC"))
    cat("For more info, try ?estStrength or ?estMC or str(obj_name, max.levels = 2).\n")
  else
    cat("For more info, try ?estMantel or str(obj_name, max.levels = 2).\n")
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
                                                               "rM")),
                                    point = c("mean", "median", "point"),
                                    range = c("simpleCI", "bcCI", "se"),
                                    xlab = NULL, ylab = plot.which,
                                    originNames = NULL, targetNames = NULL,
                                    ageNames = NULL,
                                    col = NULL, pch = NULL,
                                    gap = 0,
                                    sfrac = ifelse(range=="se", 0.01, 0),
                                    legend = FALSE, map = FALSE, ...) {
  if (map) {
    warning("Map plotting not yet available")
  }
  if ((plot.which %in% c("corr", "Mantel")))
    plot.which <- "rM"
  if ((plot.which %in% c("transition", "Transition")))
    plot.which <- "psi"
  if ((plot.which %in% c("strength", "Strength")))
    plot.which <- "MC"
  if (!(plot.which %in% c("psi", "MC", "rM", "r")))
    stop("Set plot.which to psi, MC, rM, or r")
  point <- match.arg(point)
  range <- match.arg(range)
  if (inherits(x, "estMC")) {
    if (is.null(x$psi)) {
      bcCIPsi <- array(NA, dim = c(2, dim(x$samplePsi)[2], dim(x$samplePsi)[3]))
      for (i in 1:dim(x$samplePsi)[2]) {
        for (j in 1:dim(x$samplePsi)[3]) {
          psi.z0 <- qnorm(sum(x$samplePsi[, i, j] < mean(x$samplePsi[, i, j],
                                                         na.rm = T)) /
                            length(which(!is.na(x$samplePsi[, i, j]))))
          bcCIPsi[ , i, j] <- quantile(x$samplePsi[, i, j],
                                       pnorm(2 * psi.z0 +
                                               qnorm(c(x$alpha/2,1-x$alpha/2))),
                                       na.rm=TRUE, type = 8, names = F)
        }
      }
      x$psi <- list(sample = x$samplePsi,
                    mean = apply(x$samplePsi, 2:3, mean),
                    se = apply(x$samplePsi, 2:3, sd),
                    simpleCI = apply(x$samplePsi, 2:3, quantile,
                                     probs = c(x$alpha/2, 1-x$alpha/2),
                                     na.rm=TRUE, type = 8, names = F),
                    bcCI = bcCIPsi,
                    median = apply(x$samplePsi, 2:3, median),
                    point = x$pointPsi)
      x$MC <- list(mean = x$meanMC, se = x$seMC, simpleCI = x$simpleCI,
                   bcCI = x$bcCI, median = x$medianMC, point = x$pointMC)
      x$input <- list(originNames = dimnames(x$samplePsi)[[2]],
                      targetNames = dimnames(x$samplePsi)[[3]])
    }
  }
  else if (plot.which != "rM" && !inherits(x, "estPsi"))
    stop("This estimate does not include psi or MC - try setting plot.which to rM")
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
  else if (plot.which == "psi") {
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
    nTargetSites <- length(targetNames)
    nOriginSites <- length(originNames)
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
                          From = factor(rep(originNames, nTargetSites),
                                        levels = originNames),
                          To = factor(rep(targetNames, each = nOriginSites),
                                      levels = targetNames),
                          FromTo = rep(1:nOriginSites, nTargetSites) +
                            rep(1:nTargetSites, each = nOriginSites) /
                            (nTargetSites * 2) - 0.3,
                          lower = c(yrange[1,,]),
                          upper = c(yrange[2,,]))
  }
  else if (plot.which == "r") {
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
  if (plot.which=="psi") {
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
                   sfrac = sfrac, ...)
    graphics::axis(1, at = seq(from = 1, by = 1, length.out = nOriginSites),
         labels = originNames)
    if (!isFALSE(legend))
      legend(legend, legend = targetNames, col = col, pch = pch,
             pt.bg = col)
  }
  else if (plot.which=="r") {
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
                     sfrac = sfrac, ...)
      graphics::axis(1, at = seq(from = 1, by = 1, length.out = nTargetSites),
           labels = targetNames)
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
                     sfrac = sfrac, ...)
      graphics::axis(1, at = 1:nAges, labels = ageNames)
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
                   ylab = ylab, xaxt = "n", xlab = xlab, ...)
  }
}

