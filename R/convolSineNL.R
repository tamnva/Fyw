#' Convolve sine wave isotope concentration in precipitation using gamma function
#'
#' @param AP Amplitude of the sine wave (cP = AP*sin(2*pi*t - phiP) + kP) of
#' isotope concentration in precipitation.
#' @param phiP Phase shift of the sine wave (cP = AP*sin(2*pi*t - phiP) + kP) of
#' isotope concentration in precipitation.
#' @param kP The constant factor of sine wave (cP = AP*sin(2*pi*t - phiP) + kP)
#' of isotope concentration in precipitation.
#' @param estAlpha The 'user' estimated alpha (shape) of the gamma distribution
#' function used for the convolution approach
#' @param estBeta The estimated beta (scale) of the gamma distribution function
#' use for the convolution approach
#' @param simulatedDate The date vector (must be in format "YYYY-MM-DD") and as
#' DATE of the time isotope concentration in streamflow should be estimated
#' @param nWarmupYears Number of warm-up years before the starting date of the
#' simulatedDate. This should be a few years (e.g.,5 years) to remove the effect
#' of initial condition on the estimated results.
#' @param printAll A logical variable, if TRUE then the results (simulated isotope
#' concentration in streamflow) is printed for the entire simulation period,
#' including the warm up period. If FALSE then the results is only for the
#' simulatedDate
#' @return A data frame of the date and simulated instream isotope concentration
#' @details This function convolute the input sine wave (cP = AP*sin(2*pi*t - phiP)
#' + kP) with the given gamma distribution (user given the shape and scale factor)
#' and return the simulated isotope concentrations in streamflow
#' @references
#' Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal
#' tracer cycles quantify young water fractions, but not mean transit times,
#' in spatially heterogeneous catchments, Hydrol. Earth Syst. Sci., 20, 279–297,
#' https://doi.org/10.5194/hess-20-279-2016.
#' @examples
#' simC <- convolSineNL(AP = 1.909, phiP = 2.016, kP = -11.05995, estAlpha = 0.802,
#'                      estBeta = 0.127, simulatedDate = isotopeS$date,
#'                      nWarmupYears = 5, printAll = FALSE)
#' @export

convolSineNL <- function(AP = NULL, phiP = NULL, kP = NULL, estAlpha = NULL,
                         estBeta = NULL, simulatedDate = NULL,
                         nWarmupYears = NULL, printAll = FALSE){

  # Sine function in form c = A sin (2*pi*t - phi) + k
  sineNL <- function(time, AP, phiP, kP){
    AP * sin(2*pi*time - phiP) + kP
  }

  # check input data is null
  if (is.null(simulatedDate) | is.null(AP) | is.null(phiP) | is.null(kP) |
      is.null(estAlpha)| is.null(estBeta)){
    stop("One of the input variable is NULL")
  }

  # check if nWarumUpYears is too short
  if (nWarmupYears < 5){
    stop("nWarmupYears should be larger than 5 to avoid effect of initial conditions")
  }

  # simulation date
  simDate <- seq(from = simulatedDate[1] - nWarmupYears*365,
                 to = simulatedDate[length(simulatedDate)],
                 by = "days")

  # Number of simulation dates
  nDates <- length(simDate)

  # Convert SIMulation DATE to DECimal e.g., 2001.00, 2001.99
  simDateDec <- lubridate :: decimal_date(simDate)

  # Simulation date in decimal but starting date is 0
  simDateDec0 <- simDateDec- simDateDec[1]

  # Convert OBServed DATE to DECimal
  obsDateDec <- lubridate :: decimal_date(simulatedDate)

  # Find Location of COMmon DATE between simDateDec and obsDateDec
  comDateLoc <- which(simDateDec %in% obsDateDec)

  # Convert simulation date to range of [0,1]
  simDate01 <- simDateDec - trunc(simDateDec)

  # Find the observed isotope in precipitation (estimated from sineNL functoin )
  obsIsoPre <- sineNL(simDate01, AP = AP, phiP = phiP, kP = kP)
  #plot(x = simDateDec, y = obsIsoPre, type = "l")


  # Initialized the CONVOLuted ISOtope CONCentration
  convolIsoConc <- rep(0,nDates)

  for (i in 1:nDates){
    simIsoConc <- obsIsoPre[i] * pgamma(simDateDec0, shape = estAlpha,
                                        scale = estBeta, lower = TRUE)
    simIsoConc <- diff(c(0,simIsoConc))
    convolIsoConc[i:nDates] <- convolIsoConc[i:nDates] + simIsoConc[1:(nDates - i + 1)]
  }

  if(printAll){
    output <- data.frame(date = simDate, simIsoConc = convolIsoConc)
  } else {
    output <- data.frame(date = simulatedDate,simIsoConc = convolIsoConc[comDateLoc])
  }

  return(output)

}
