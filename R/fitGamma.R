#' Fit gamma distribution to observed instream tracer concentratoin
#'
#' @param AP Amplitude of the sine wave (cP = AP*sin(2*pi*t - phiP) + kP) of
#' isotope concentration in precipitation.
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
#'
#' #Get isotope data in streamflow of the Alp catchment (from the example dataset)
#' isotopeS_Alp <- subset(isotopeData, catchment == "Alp" & variable == "streamflow")
#'
#' simC <- fitGamma(AP = 1.909, phiP = 2.016, kP = -11.05995, alphaRange = c(0.01,5),
#'                      betaRange = c(0.01,10), simulatedDate = isotopeS_Alp$date,
#'                      fittedData = isotopeS_Alp$delta_18O,
#'                      weight = isotopeS_Alp$water_flux_mm, nIter = 8,
#'                      nBestIter = 3, nCores = 1,
#'                      nWarmupYears = 5)
#' @export

fitGamma <- function(AP = NULL, phiP = NULL, kP = NULL, alphaRange = NULL,
                     betaRange = NULL, simulatedDate = NULL, fittedData = NULL,
                     weight = NULL, nIter = 5000, nBestIter = 30, nCores = 1,
                     nWarmupYears = NULL){

  # Generate nIter number of parametersets within the range of [0,1]
  parameterSet <- lhs :: randomLHS(n = nIter, k = 2)

  # Convert the generated parameter set to user-defined range
  parameterSet[,1] <- alphaRange[1] + parameterSet[,1] * (alphaRange[2] - alphaRange[1])
  parameterSet[,2] <- betaRange[1] + parameterSet[,2] * (betaRange[2] - betaRange[1])

  # Make cluster
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)

  # Get instream tracer concentrations by running on nCores
  simulatedC <- foreach(i = 1:nIter, .combine = cbind, .export=c("convolSineNL")) %dopar% {
    simC <- convolSineNL(AP = AP, phiP = phiP, kP = kP, estAlpha = parameterSet[i,1],
                         estBeta = parameterSet[i,2], simulatedDate = simulatedDate,
                         nWarmupYears = nWarmupYears, printAll = FALSE)
    simC$simIsoConc
  }

  # Close cluster
  stopCluster(cl)

  # Calculate square error
  squareErr <- rep(NA, nIter)
  for (i in 1:nIter){
    if(is.null(weight)){
      squareErr[i] <- sum((simulatedC[,i] - fittedData)^2)
    } else {
      squareErr[i] <- sum(weight * (simulatedC[,i] - fittedData)^2)
    }
  }

  # Ranking results by decreasing square error
  rankingDecreasing <- order(squareErr, decreasing = FALSE)

  # Sort simulated results by decreasing square error
  simulatedC <- as.data.frame(simulatedC[,rankingDecreasing])
  colnames(simulatedC) <- paste0("simulation_", c(1:ncol(simulatedC)))

  date <- simulatedDate
  parameterSet <- as.data.frame(parameterSet[rankingDecreasing,])
  colnames(parameterSet) <- c("alpha", "beta")

  # Create output object
  output <- list()
  output$simulatedC <- cbind(date, simulatedC[,1:nBestIter])
  output$parameterSet <- parameterSet[1:nBestIter,]

  return(output)
}
