#' Find \eqn{\varphi} within the range of \eqn{[0, 2\pi]}
#'
#' @param aP A first scalar (see Details)
#' @param bP A second scalar (see Details)
#'
#' @return A single phiP value within the range of \eqn{[0, 2\pi]}.
#'
#' @details
#' Given the sine function (e.g., Eq. 4 in Kirchner, (2016)): \cr
#' \deqn{c_P = A_P \cdot sin(2 \pi t - \varphi_P) + k_P}{}
#'
#' which can be rewritten as follows (e.g., Eq. 5 in Kirchner, (2016)): \cr
#'  \deqn{c_P = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}{}
#'
#' where: \cr
#' \deqn{a_P = -A_p \cdot sin(\varphi_P)}{}
#' \deqn{b_P = A_p \cdot cos(\varphi_P)}{}
#'
#' Assuming that \eqn{a_P, b_P} are known, find \eqn{\varphi_P} within the range
#' of \eqn{[0, 2\pi]}{}
#'
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems –
#' Part 1: Seasonal tracer cycles quantify young water fractions, but not mean
#' transit times, in spatially heterogeneous catchments. \emph{Hydrol. Earth Syst. Sci.},
#' 20, 279–297. \url{https://doi.org/10.5194/hess-20-279-2016}.
#'
#' @examples
#' Ap <- 2.5
#' phiP <- 0.1
#'
#' findPhi(aP = -Ap*sin(phiP), bP = Ap*cos(phiP))
#'
#' # This function should return phiP = 0.1
#'
#' phiP <- 2*pi - 0.1
#'
#' findPhi(aP = -Ap*sin(phiP), bP = Ap*cos(phiP))
#'
#' # This function should return phiP = 2*pi - 0.1
#'
#' @export

findPhi <- function(aP = NULL, bP = NULL){

  # Check if input is null
  if (is.null(aP) | is.null(bP)) stop("aP and bP cannot be NULL")

  # Find first value of phi
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
