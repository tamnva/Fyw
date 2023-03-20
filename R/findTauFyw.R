#' Find findtauyw using bisection approach
#'
#' @param shape Shape factor (or alpha) of the gamma distribution
#' @param scale Scale factor (or beta) of the gamma distribution
#' phiS MUST be smaller than phiP
#' @param Fyw Youngwater fraction, equals to the amplitude ratio AS/AP where AS
#' and AP are the amplitudes of the fitted sine wave to isotope concentrations in
#' streamflow and in precipitation, respectively.
#' @return Return the age threshold of youngwater fraction tauyw in years.
#' @details In this function, we solve equation 13 (Kirchner et al., 2013) with
#' given shape (alpha), scale (beta), and youngwater fraction (Fyw). The result is
#' the age threshold of youngwater fraction (tauyw) in years
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' findtauyw(shape = 1.0, scale = 2.0, Fyw = 0.25)
#' @export

findtauyw <- function(shape = NULL, scale = NULL, Fyw = NULL){

  # check shape factor is null
  if (is.null(shape)){
    stop("shape factor cannot be NULL")
  }

  # check if scale factor is null
  if (is.null(scale)){
    stop("scale factor cannot be NULL")
  }

  # check if Fyw is null
  if (is.null(Fyw)){
    stop("Fyw cannot be NULL")
  }

  # check if Fyw is out of the [0,1] range
  if ((Fyw >= 1) | (Fyw <= 0)){
    stop("Fyw must be within the rang of 0 and 1")
  }

  # define function need to be found the root
  f <- function(tauyw){
    pgamma(tauyw, shape = shape, scale = scale, lower = TRUE) - Fyw
  }

  tauywMin <- 0
  tauywMax <- 100.0
  tauyw <- (tauywMax + tauywMin)/2

  minVal <- f(tauywMin)
  maxVal <- f(tauywMax)
  midVal <- f(tauyw)

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
    midVal <- f(tauyw)
  }

  return(tauyw)
}
