#' @export
is.isoAssign <- function(x) inherits(x, "isoAssign")
is.estMC <- function(x) inherits(x, "estMC")

#' @export
print <- function(x,...) UseMethod("print")
print.estMigConnectivity <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Migratory Connectivity Strength Estimate(s)\n")
  if (inherits(x, "estMC")) {
    cat("MC estimate (mean):", format(x$meanMC, digits = digits), "+/- (SE)",
        format(x$seMC, digits = digits), '\n')
    cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
        "% confidence interval (simple quantile): ",
        paste(format(x$simpleCI, digits = digits, trim = TRUE),
              collapse = ' - '), '\n', sep = "")
    cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
        "% confidence interval (bias-corrected): ",
        paste(format(x$bcCI, digits = digits, trim = TRUE), collapse = ' - '),
        '\n', sep = "")
    cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
        "% credible interval (highest posterior density): ",
        paste(format(x$hpdCI, digits = digits, trim = TRUE), collapse = ' - '),
        '\n', sep = "")
    cat("   median:", format(x$medianMC, digits = digits), '\n')
    if (!is.na(x$pointMC))
      cat("   point calculation (not considering error):",
          format(x$pointMC, digits = digits), '\n')
  }
  if (inherits(x, "estMantel")) {
    cat("rM estimate (mean):", format(x$meanCorr, digits = digits), "+/- (SE)",
        format(x$seCorr, digits = digits), '\n')
    cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
        "% confidence interval (simple quantile): ",
        paste(format(x$simpleCICorr, digits = digits, trim = TRUE),
              collapse = ' - '), '\n', sep = "")
    cat("   ", ifelse(is.null(x$alpha), "", 100 * (1 - x$alpha)),
        "% confidence interval (bias-corrected): ",
        paste(format(x$bcCICorr, digits = digits, trim = TRUE), collapse = ' - '),
        '\n', sep = "")
    cat("   median:", format(x$medianCorr, digits = digits), '\n')
    if (!is.na(x$pointCorr))
      cat("   point calculation (not considering error):",
          format(x$pointCorr, digits = digits), '\n')
  }
}

#' @export
summary <- function(x,...) UseMethod("summary")

summary.estMigConnectivity <- function(x, ...){print.estMigConnectivity(x, ...)}

summary.isoAssign<-function(x, ...){
  cat("Individual Probability Surfaces \n")
  print(x$probassign,...)
  cat("\n Individual likely/unlikely Surfaces \n")
  print(x$oddsassign,...)
  cat("\n Population-level assignment Surface \n")
  print(x$popassign,...)
  cat("\n Individual Probability data frame* \n")
  print(x$probDF[1:6,1:5])
  cat("\n Individual likely/unlikely data frame* \n")
  print(x$oddsDF[1:6,1:5])
  cat("\n Individual single cell assignment \n")
  str(x$SingleCell)
  cat("\n Random number seed set to: \n")
  print(x$RandomSeed)
  cat("\n * only first few columns are printed")

}

#' basic plot function for the different isoAssign outputs
#' @export
plot <- function(x,...) UseMethod("plot")
plot.isoAssign <- function(x,map,...){
  if(!(map %in% c("probability","population","odds"))){
    stop("map must be either probability, population, or odds")}
  op <- par(no.readonly = TRUE)
  if(map == "population"){
    raster::plot(x$popassign,horiz = TRUE,...)
  }
  if(map == "probability"){
    for(i in 1:raster::nlayers(x$probassign)){
      raster::plot(x$probassign[[i]],horiz = TRUE,...)
      par(ask = TRUE)
    }
    par(op)
  }
  if(map == "odds"){
    for(i in 1:raster::nlayers(x$probassign)){
      raster::plot(x$oddsassign[[i]],horiz = TRUE,...)
      par(ask = TRUE)
    }
    par(op)
  }
  on.exit(par(op))
}
