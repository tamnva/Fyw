#' Find alpha value using bisection approach
#'
#' @param phiP Phase shift of isotope concentration in precipitation (in radian).
#' phiP MUST be within the range of [0, 2*pi]
#' @param phiS Phase shift of isotope concentration in streamflow (in radian). \cr
#' phiS MUST be within the range of [0, 2*pi] \cr
#' phiS MUST be smaller than phiP
#' @param Fyw Youngwater fraction, equals to the amplitude ratio AS/AP where AS
#' and AP are the amplitudes of the fitted sine wave to isotope concentrations in
#' streamflow and in precipitation, respectively.
#' @param eps The minimum search distance
#' @param maxIter Maximum number of iteration when using bisection approach
#
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' findAlpha <- function(phiS = 2*pi, phiP = pi, Fyw = 0.3, eps = 1e-6,
#' maxIter = 1000)
#
findAlpha <- function(phiS = NULL, phiP = NULL, Fyw = NULL, eps = 1e-6,
                      maxIter = 1000){

  # alpha function (Eq. 11 in Kirchner 2016, HESS)
  f <- function(alpha){
    phiS - phiP - alpha * atan(sqrt(Fyw^(-2/alpha) - 1))
  }

  # Using bisection algorithm, search range 1.0e-6 and 20
  alphaMin <- eps
  alphaMax <- 20.0
  alpha <- (alphaMax + alphaMin)/2

  minVal <- f(alphaMin)
  maxVal <- f(alphaMax)
  midVal <- f(alpha)

  iter <- 0

  while (((alphaMax - alphaMin) > eps) & (iter < maxIter)){
    # check if the function does not change sign within the search range
    if ((minVal * midVal) > 0 & (midVal * maxVal) > 0){
      stop("Cannot find root of this equation within the given alpha range")
      alpha <- NA
    }

    # check if found exact solution
    if (midVal == 0){
      stop("Find exact solution")
    }

    # update min max aplpha value
    if ((minVal * midVal) < 0){
      alphaMax <- alpha
      maxVal <- midVal
    } else {
      alphaMin <- alpha
      minVal <- midVal
    }

    alpha <- (alphaMax + alphaMin)/2
    midVal <- f(alpha)
    iter <- iter + 1
  }

  return(alpha)
}
