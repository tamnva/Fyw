% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitSineNL.R
\name{fitSineNL}
\alias{fitSineNL}
\title{Non-linear sine wave fitting}
\usage{
fitSineNL(
  observed = NULL,
  a = c(0, 10),
  phi = c(0, 2 * pi),
  k = c(-20, 0),
  t = NULL,
  nIter = 50000,
  nBestIter = 30,
  weight = rep(1, length(observed)),
  setSeed = NULL
)
}
\arguments{
\item{observed}{A numeric vector containing values of the dependent variables, this
could be observed isotope in precipitation (or streamflow).}

\item{a}{A numeric vector with positive numbers, representing the maximum and
minimum search range for the amplitude A of isotope cycle in precipitation
(or streamflow).}

\item{phi}{A numeric vector, representing the maximum and minimum search range
for the phase shift of isotope  in precipitation (or streamflow). The range
must be within \eqn{[0, 2\pi]}.}

\item{k}{A numeric vector, representing the maximum and minimum search range
for the phase shift of the constant value in the fitted function.}

\item{t}{A date vector in format format yyyy-mm-dd (e.g., 2020-12-31) representing
the date of observed isotope.}

\item{nIter}{Maximum number of iterations (number of parametersets), higher
number higher chance to have a better fit with observed data, however, longer
computation time.}

\item{nBestIter}{Number of best parameter sets. For example, when running with
50000 parameter sets, the function will rank which parameter set give a better
fit with observed based on square error, then save all outputs from nBestIter.}

\item{weight}{A weight vector (same length with observed) when users want to weight
with streamflow or precipitation (to have the weighted Fyw). Otherwise, if no
input is given, all of the weight is 1 (meaning equal weights/unweighted - this
case for calculating the unweighted Fyw).}

\item{setSeed}{An integer number used for set.seed() to have a reproducible
results.}
}
\value{
A list of object containing four tables: (1) a table of simulated results (e.g.,
isotope concentration in precipitation and streamflow), (2) a table of observed data
used for fitting, (3) a table of the goodness-of-fit (here is the weighted or unweighted
mean square error between simulated and observed data), and (4) a table of the parameters
of the fitted sine-wave.
}
\description{
Fitting a non-linear sine-wave function to observed data.
}
\details{
The non-linear sine wave function used for fitting has the following
form (Eq. 5; Kirchner, 2016): \cr

\deqn{c = a \cdot sin(2 \pi t - \varphi_P) + k_P}{}
}
\examples{
# Get isotope data in streamflow of the Alp catchment (from the example dataset)
isotopeS_Alp <- subset(isotopeData, catchment == "Alp" & variable == "streamflow")

# Fit a non-linear sine wave function to observed isotope data in streamflow
fitSineS <- fitSineNL(observed = isotopeS_Alp$delta_18O,
                      a = c(0,10),
                      phi = c(0, 2*pi),
                      k = c(-20,0),
                      t = isotopeS_Alp$date,
                      nIter = 5000,
                      nBestIter = 2,
                      weight = isotopeS_Alp$water_flux_mm)
}
\references{
Kirchner, J. W. (2016). Aggregation in environmental systems –
Part 1: Seasonal tracer cycles quantify young water fractions, but not mean
transit times, in spatially heterogeneous catchments. \emph{Hydrol. Earth Syst. Sci.},
20, 279–297.
}
