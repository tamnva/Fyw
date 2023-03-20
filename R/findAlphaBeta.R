#' Find findAlphaBeta using bisection approach
#' @description Find alpha (shape) and beta (scale) factor of the gamma
#' distriubtion using the bisection approach. Solve equations (10) and (11) in
#' Kirchner (2016).
#' @param phiP Phase shift of isotope concentration in precipitation (in radian).
#' phiP MUST be within the range of [0, 2*pi]
#' @param phiS Phase shift of isotope concentration in streamflow (in radian). \cr
#' phiS MUST be within the range of [0, 2*pi] \cr
#' phiS MUST be smaller than phiP (will be automatically checked inside this function)
#' @param Fyw Youngwater fraction, equals to the amplitude ratio AS/AP where AS
#' and AP are the amplitudes of the fitted sine wave to isotope concentrations in
#' streamflow and in precipitation, respectively.
#' @param eps The minimum search distance
#' @param maxIter Maximum number of iteration when using bisection approach
#' @return A list object containing the shape (alpha), scale (beta) factors,
#' and also mean transit time
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' result <- findAlphaBeta(phiS = 2*pi, phiP = pi, Fyw = 0.3, eps = 1e-6, maxIter = 1000)
#' # alpha (shape) factor
#' result$alpha
#' # beta (scale) factor
#' result$beta
#' # Mean transit time (years)
#' result$meanTT
#' @export

findAlphaBeta <- function(phiS = NULL, phiP = NULL, Fyw = NULL, eps = 1e-6){

  # NOTE: Shape factor is alpha, scale factor is beta

  # alpha function (Eq. 11 in Kirchner 2016, HESS)
  f <- function(alpha){
    phiS - phiP - alpha * atan(sqrt(Fyw^(-2/alpha) - 1))
  }

  # Check if phase shift in isotope signal
  if(phiS < phiP){
    alpha <- NA
    stop("phiS MUST BE > phiP")
  }

  # Check if Fyw > 1
  if(Fyw > 1){
    alpha <- NA
    stop("Fyw must be smaller or equalt to 1")
  }

  # Using bisection algorithm, search range 1.0e-6 and 20
  alphaMin <- eps
  alphaMax <- 100.0
  alpha <- (alphaMax + alphaMin)/2

  minVal <- f(alphaMin)
  maxVal <- f(alphaMax)
  midVal <- f(alpha)

  while (alphaMax - alphaMin > eps){
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
  }

  output <- list()
  output$alpha <- alpha

  # Also find beta using equation 10 (Kirchner, 2016), f = 1 for seasonal cycle
  output$beta <- (1/(2*pi*1))*sqrt(Fyw^(-2/alpha) - 1)
  output$meanTT <- output$alpha * output$beta

  return(output)
}
