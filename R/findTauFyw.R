#' Find findtauyw using bisection approach
#'
#' @param alpha A scalar or vector of the shape factor of the gamma distribution
#' @param beta A scalar or vector of the scale factor of the gamma distribution
#' phiS MUST be smaller than phiP
#' @param Fyw A scalar or vector of the o the youngwater fraction
#' @param eps The minimum search distance, default value is 1e-6. In other words,
#' if the next tauyw values changes within this range compared to the previous value
#' the function will stop and return the next tauyw value
#' @return Return a data frame object of the age threshold of youngwater fraction
#' (tauyw) in years. This data frame also contains input.
#' @details In this function, we solve equation 13 (Kirchner et al., 2013) with
#' given alpha (shape), beta (scale), and youngwater fraction (Fyw). The result is
#' the age threshold of youngwater fraction (tauyw) in years
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' findtauyw(alpha = 1.0,
#'           beta = 2.0,
#'           Fyw = 0.25)
#' @export

findtauyw <- function(alpha = NULL,
                      beta = NULL,
                      Fyw = NULL,
                      eps = 1e-6){

  # check alpha factor is null
  if (is.null(alpha)){
    stop("alpha factor cannot be NULL")
  }

  # check if beta factor is null
  if (is.null(beta)){
    stop("beta factor cannot be NULL")
  }

  # check if Fyw is null
  if (is.null(Fyw)){
    stop("Fyw cannot be NULL")
  }

  # check if all Fyw <= 1
  if (!all(Fyw <= 1)){
    stop("Fyw must be <= 1")
  }

  # check if all Fyw >= 0
  if (!all(Fyw >= 0)){
    stop("Fyw must be >= 0")
  }

  # check if length of alpha, beta, and Fyw are equal
  if (length(alpha) != length(beta) | length(alpha) != length(Fyw)){
    stop("Input vectors alpha, beta, and Fyw must have the same length")
  }

  # define function need to be found the root
  f <- function(tauyw, i){
    pgamma(tauyw, shape = alpha[i], scale = beta[i], lower = TRUE) - Fyw[i]
  }

  # Initial output (tauyw)
  output <- c()

  # Loop over the length of alpha
  for (i in 1:length(alpha)){
    tauywMin <- 0
    tauywMax <- 100.0
    tauyw <- (tauywMax + tauywMin)/2

    minVal <- f(tauywMin, i)
    maxVal <- f(tauywMax, i)
    midVal <- f(tauyw, i)

    while (tauywMax - tauywMin > eps){

      if (maxVal * minVal > 0){
        stop("Cannot find tauyw within the range of [0,100] years")
      }

      # check if found exact solution
      if (midVal == 0){
        stop("Find exact solution")
      }

      # update min max aplpha value
      if ((minVal * midVal) < 0){
        tauywMax <- tauyw
        maxVal <- midVal
      } else {
        tauywMin <- tauyw
        minVal <- midVal
      }

      tauyw <- (tauywMax + tauywMin)/2
      midVal <- f(tauyw, i)
    }

    output <- c(output, tauyw)
  }

  # output as data frame, including input
  output <- data.frame(alpha = alpha, beta = beta, Fyw = Fyw, tauyw = output)

  return(output)
}
