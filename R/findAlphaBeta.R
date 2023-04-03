#' Find findAlphaBeta using bisection approach
#' @description Find alpha (shape) and beta (scale) factors of the gamma
#' distriubtion using the bisection approach. Solve equations (10) and (11) in
#' Kirchner (2016).
#' @param phiS_phiP Phase shift of isotope concentration in streamflow (in radian). \cr
#' phiS_phiP MUST be > 0
#' @param AS The amplitude of the sine-wave function of isotope concentration in
#' precipitation. AS MUST be > 0
#' @param AP The amplitude of the sine-wave function of isotope concentration in
#' streamflow. AP MUST be > 0
#' @param eps The minimum search distance
#' @return A list object containing the shape (alpha), scale (beta) factors,
#' and also mean transit time
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' result <- findAlphaBeta(phiS_phiP = pi,
#'                         AS = 1.0,
#'                         AP = 2.0,
#'                         eps = 1e-6)
#' # alpha (shape) factor
#' result$alpha
#' # beta (scale) factor
#' result$beta
#' # Mean transit time (years)
#' result$meanTT
#' @export

findAlphaBeta <- function(phiS_phiP = NULL,
                          AS = NULL,
                          AP = NULL,
                          eps = 1e-6){

  # NOTE: Shape factor is alpha, scale factor is beta

  # alpha function (Eq. 11 in Kirchner 2016, HESS)
  f <- function(alpha){
    phiS_phiP - alpha * atan(sqrt((AS/AP)^(-2/alpha) - 1))
  }

  # Check if phase shift in isotope signal
  if(phiS_phiP <= 0){
    stop("phiS_phiP must be > 0")
  }

  # Check if Fyw > 1
  if(AS/AP > 1){
    alpha <- NA
    stop("Fyw must be < 1")
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
  output$beta <- (1/(2*pi*1))*sqrt((AS/AP)^(-2/alpha) - 1)
  output$meanTT <- output$alpha * output$beta

  return(output)
}
