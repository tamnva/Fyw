#' Find phiP within the range of 0 and 2pi
#'
#' @param aP A scalar or vector of the shape factor of the gamma distribution
#' @param bP A scalar or vector of the scale factor of the gamma distribution
#' phiS MUST be smaller than phiP
#' @return A single phiP value within the raange of [0, 2*pi]
#' @details This function find root of the following system of equations for phiP
#' that within the range of [0, 2*pi] given values of aP and bP
#' aP <- - sqrt(aP^2 + bP^2) * sin(phiP) = -Ap*sin(phiP)
#' bP <- sqrt(aP^2 + bP^2) * cos(phiP) = Ap*cos(phiP)
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' Ap <- 2.5
#' phiP <- 0.1
#' findPhi(aP = -Ap*sin(phiP), bP = Ap*cos(phiP))
#' # This function should return phiP = 0.1
#' phiP <- 2*pi - 0.1
#' findPhi(aP = -Ap*sin(phiP), bP = Ap*cos(phiP))
#' # This function should return phiP = 2*pi - 0.1
#' @export

findPhi <- function(aP = NULL, bP = NULL){
  phiP <- atan(-aP/bP)

  # if phiP is negative then get positive value
  if (phiP < 0){
    if (phiP < - pi){
      phiP <- phiP + 2*pi
    } else {
      phiP <- phiP + pi
    }
  }

  # Then get two solutions of pi within [0, 2*pi]
  if (sin(phiP)*aP > 0){
    if(phiP + pi > 2*pi){
      phiP <- phiP - pi
    } else {
      phiP <- phiP + pi
    }
  }

  # Return result
  return(phiP)

}
