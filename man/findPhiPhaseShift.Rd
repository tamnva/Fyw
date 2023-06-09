% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findPhiPhaseShift.R
\name{findPhiPhaseShift}
\alias{findPhiPhaseShift}
\title{Find phase and phase shift}
\usage{
findPhiPhaseShift(aP = NULL, bP = NULL, aS = NULL, bS = NULL)
}
\arguments{
\item{aP}{The slope coefficient of the linear sine wave: \cr
\eqn{c_P = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}.}

\item{bP}{The slope coefficient of the linear sine wave: \cr
\eqn{c_P = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}.}
}
\value{
A list object containing \eqn{\varphi_P} (rad), \eqn{\varphi_S} (rad),
and \eqn{\varphi_S - \varphi_P} (rad)
}
\description{
Find phase of the sine-wave describing isotope in precipitation (\eqn{\varphi_P})
and streamflow (\eqn{\varphi_S}) and the phase shift \eqn{\varphi_S - \varphi_P}
}
\details{
Given the sine functions describing isotope in precipitation \eqn{c_P(t)} and
streamflow \eqn{c_S(t)} (e.g., Eq. 4 in Kirchner, (2016)): \cr

\deqn{c_P(t) = A_P \cdot sin(2 \pi t - \varphi_P) + k_P}{}
\deqn{c_S(t) = A_S \cdot sin(2 \pi t - \varphi_S) + k_S}{}

which can be rewritten as follows (e.g., Eq. 5 in Kirchner, (2016)): \cr
 \deqn{c_P(t) = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P}{}
 \deqn{c_S(t) = a_S \cdot cos(2 \pi t) + b_S \cdot sin(2 \pi t) + k_S}{}

where: \cr
\deqn{a_P = -A_P \cdot sin(\varphi_P)}{}
\deqn{b_P = A_P \cdot cos(\varphi_P)}{}

and
\deqn{a_S = -A_S \cdot sin(\varphi_S)}{}
\deqn{b_S = A_S \cdot cos(\varphi_S)}{}

Assuming that \eqn{a_P, b_P, a_S, b_S} are known, find \eqn{\varphi_P} within
the range of \eqn{[0, 2\pi]}, \eqn{\varphi_S} within the range of \eqn{[0, 2\pi]}
if \eqn{\varphi_S > \varphi_P}, otherwise, find  \eqn{\varphi_S} within the
range of \eqn{[2\pi, 4\pi]}
}
\examples{
AP <- 2.5
AS <- 2.0
phiP <- 0.5*pi
phiS <- 1.5*pi
findPhiPhaseShift(aP = -AP*sin(phiP),
                  bP = AP*cos(phiP),
                  aS = -AS*sin(phiS),
                  bS = AS*cos(phiS))

AP <- 2.5
AS <- 2.0
phiP <- pi
phiS <- 2.5*pi
findPhiPhaseShift(aP = -AP*sin(phiP),
                  bP = AP*cos(phiP),
                  aS = -AS*sin(phiS),
                  bS = AS*cos(phiS))

}
\references{
Kirchner, J. W. (2016). Aggregation in environmental systems –
Part 1: Seasonal tracer cycles quantify young water fractions, but not mean
transit times, in spatially heterogeneous catchments. \emph{Hydrol. Earth Syst. Sci.},
20, 279–297.
}
