#' @export
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
summary.estMigConnectivity <- function(x, ...)
{
  print.estMigConnectivity(x, ...)
}

