#' Find phase and phase shift
#'
#' @description
#' Find phase of the sine-wave describing isotope in precipitation (\eqn{\varphi_P})
#' and streamflow (\eqn{\varphi_S}) and the phase shift \eqn{\varphi_S - \varphi_P}
#'
#'
#' @param aP The slope coefficient of the linear sine wave: \cr
#' \eqn{c_P = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}.
#'
#' @param bP The slope coefficient of the linear sine wave: \cr
#' \eqn{c_P = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}.
#'
#' @return A list object containing \eqn{\varphi_P} (rad), \eqn{\varphi_S} (rad),
#' and \eqn{\varphi_S - \varphi_P} (rad)
#'
#' @details
#' Given the sine functions describing isotope in precipitation \eqn{c_P(t)} and
#' streamflow \eqn{c_S(t)} (e.g., Eq. 4 in Kirchner, (2016)): \cr
#'
#' \deqn{c_P(t) = A_P \cdot sin(2 \pi t - \varphi_P) + k_P}{}
#' \deqn{c_S(t) = A_S \cdot sin(2 \pi t - \varphi_S) + k_S}{}
#'
#' which can be rewritten as follows (e.g., Eq. 5 in Kirchner, (2016)): \cr
#'  \deqn{c_P(t) = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}{}
#'  \deqn{c_S(t) = a_S \cdot cos(2 \pi t) + b_S \cdot sin(2 \pi t) + k_S}{}
#'
#' where: \cr
#' \deqn{a_P = -A_P \cdot sin(\varphi_P)}{}
#' \deqn{b_P = A_P \cdot cos(\varphi_P)}{}
#'
#' and
#' \deqn{a_S = -A_S \cdot sin(\varphi_S)}{}
#' \deqn{b_S = A_S \cdot cos(\varphi_S)}{}
#'
#' Assuming that \eqn{a_P, b_P, a_S, b_S} are known, find \eqn{\varphi_P} within
#' the range of \eqn{[0, 2\pi]}, \eqn{\varphi_S} within the range of \eqn{[0, 2\pi]}
#' if \eqn{\varphi_S > \varphi_P}, otherwise, find  \eqn{\varphi_S} within the
#' range of \eqn{[2\pi, 4\pi]}
#'
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems –
#' Part 1: Seasonal tracer cycles quantify young water fractions, but not mean
#' transit times, in spatially heterogeneous catchments. \emph{Hydrol. Earth Syst. Sci.},
#' 20, 279–297.
#'
#' @examples
#' AP <- 2.5
#' AS <- 2.0
#' phiP <- 0.5*pi
#' phiS <- 1.5*pi
#' findPhiPhaseShift(aP = -AP*sin(phiP),
#'                   bP = AP*cos(phiP),
#'                   aS = -AS*sin(phiS),
#'                   bS = AS*cos(phiS))
#'
#' AP <- 2.5
#' AS <- 2.0
#' phiP <- pi
#' phiS <- 2.5*pi
#' findPhiPhaseShift(aP = -AP*sin(phiP),
#'                   bP = AP*cos(phiP),
#'                   aS = -AS*sin(phiS),
#'                   bS = AS*cos(phiS))
#'
#' @export

findPhiPhaseShift <- function(aP = NULL, bP = NULL, aS = NULL, bS = NULL){

  # output as list object
  output <- list()

  # Check if input is null
  if (is.null(aP) | is.null(bP) | is.null(aS) | is.null(bS)) {
    stop("aP or bP or aS or bS cannot be NULL")
  }

  # Function to find phi between [0, 2pi]
  findPhi <- function(a = NULL, b = NULL){

    # Find first value of phi
    phi <- atan(-a/b)

    # if phi is negative then get positive value
    if (phi < 0){
      if (phi < - pi){
        phi <- phi + 2*pi
      } else {
        phi <- phi + pi
      }
    }

    # Then get two solutions of pi within [0, 2*pi]
    if (sin(phi)*a > 0){
      if(phi + pi > 2*pi){
        phi <- phi - pi
      } else {
        phi <- phi + pi
      }
    }

    # Return result
    return(phi)

  }

  # Find phiP within [0, 2pi]
  output$phiP <- findPhi(a = aP, b = bP)

  # Find phiS within [0, 2pi]
  output$phiS <- findPhi(a = aS, b = bS)

  # Update phiS if phiS < phiP (the output signal appear in the year later)
  # Please see the vignettes for a detail explaination
  if (output$phiS < output$phiP) output$phiS <- output$phiS + 2*pi

  # Phase shift
  output$phiS_phiP <- output$phiS - output$phiP

  # Return output
  return(output)
}
