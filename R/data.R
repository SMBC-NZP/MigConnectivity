#' Example transition probabilities (psis) between origin and target sites
#'
#' A dataset containing example psi matrices used in Cohen et al. (in prep.).
#'
#' @format A named list with 8 transition probability matrices in it. The
#' direction is from origin site (rows) to target sites (columns), so each
#' row of each matrix sums to 1.  The psi matrices are:
#' \itemize{
#'   \item Full Mix: Full mixing from all origin sites to all target sites
#'   \item Avoid One Site: All origin sites have the same transition
#'      probabilities, mostly avoiding target site 4
#'   \item Full Connectivity: Each origin site transitions to only one target
#'      site
#'   \item Half Mix: Origin sites A and B mix fully between target sites 1 and
#'      2, but don't move to target sites 3 or 4, while origin sites C and D
#'      mix fully between target sites 3 and 4, but don't move to target sites
#'      1 or 2
#'   \item Low: Simulation scenario labelled "Moderate Connectivity" in Cohen
#'      et al. (2014)
#'   \item Medium: Simulation scenario labelled "Strong Connectivity" in Cohen
#'      et al. (2014)
#'   \item One Site Preference: Three origin sites have full mixing, but origin
#'      site D only goes to target site 4
#'   \item Negative: Artificial transition probability scenario developed to
#'      produce a negative MC value under some circumstances
#' }
"samplePsis"

#' Example origin and target site positions and distances on a 2-D plane
#'
#' \code{sampleOriginPos} is a dataset containing example origin site positions
#' from 12 scenarios used in Cohen et al. (in prep.).  For the same 12
#' scenarios, \code{sampleOriginDist} contains the origin site distances,
#' \code{sampleTargetPos} contains the target site positions, and
#' \code{sampleTargetDist} contains the target site distances.
#'
#' @format Each dataset is a named list with 12 matrices in it, representing 12
#' scenarios. The position matrices each have 2 columns (x and y position) and
#' 4 rows (for each origin or target site).  The distance matrices are
#' symmetrical and 4 x 4. The 12 scenarios are:
#' \itemize{
#'   \item Linear: Both origin and target sites arranged in horizontal linear
#'      fashion, with equal distances between each adjacent site
#'   \item B Dist BC*2: Linear, but the central origin sites are twice as far
#'      from each other as the edge sites are from the adjacent origin sites
#'   \item B Dist BC/2: Linear, but the central origin sites are half as far
#'      from each other as the edge sites are from the adjacent origin sites
#'   \item B Dist CD*2: Linear, but the last two origin sites are twice as far
#'      from each other as the other adjacent origin sites
#'   \item B Dist CD/2: Linear, but the last two origin sites are half as far
#'      from each other as the other adjacent origin sites
#'   \item B Grid: Origin sites arranged on a grid, target sites arranged
#'      linearly, both with all adjacent sites (excluding diagonals)
#'      equidistant
#'   \item NB Dist 23*2: Linear, but the central target sites are twice as far
#'      from each other as the edge sites are from the adjacent target sites
#'   \item NB Dist 23/2: Linear, but the central target sites are half as far
#'      from each other as the edge sites are from the adjacent target sites
#'   \item NB Dist 34*2: Linear, but the last two target sites are twice as far
#'      from each other as the other adjacent target sites
#'   \item NB Dist 34/2: Linear, but the last two target sites are half as far
#'      from each other as the other adjacent target sites
#'   \item NB Grid: Target sites arranged on a grid, origin sites arranged
#'      linearly, both with all adjacent sites (excluding diagonals)
#'      equidistant
#'   \item B/NB Grid: Origin and target sites each arranged on a grid, both
#'      with all adjacent sites (excluding diagonals) equidistant
#' }
"sampleOriginPos"

#' @rdname sampleOriginPos
"sampleOriginDist"

#' @rdname sampleOriginPos
"sampleTargetPos"

#' @rdname sampleOriginPos
"sampleTargetDist"

#' Example origin site abundances and relative abundances
#'
#' \code{sampleOriginN} is a dataset containing example origin site abundances
#' from 5 scenarios used in Cohen et al. (in prep.).  For the same 5
#' scenarios, \code{sampleOriginRelN} contains the relative abundances.
#'
#' @format Each dataset is a named list with 5 vectors in it. Each vector has
#' 4 elements (for the 4 origin sites).  The relative abundance vectors each
#' sum to 1. The 5 scenarios are:
#' \itemize{
#'   \item Base: Equal abundance at each origin site
#'   \item B Doub: The second origin site has twice the abundance of the other
#'      three sites
#'   \item B Half: The second origin site has half the abundance of the other
#'      three sites
#'   \item D Doub: The last origin site has twice the abundance of the other
#'      three sites
#'   \item D Half: The last origin site has half the abundance of the other
#'      three sites
#' }
"sampleOriginN"

#' @rdname sampleOriginN
"sampleOriginRelN"
