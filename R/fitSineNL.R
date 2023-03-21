#' Fit a nonlinear sine wave function
#'
#' @param obsC A numeric vectors containing values of the dependent variables, this
#' could be measure isotope concentration in precipitation (or streamflow).
#' @param A A numeric vector with positive numbers, representing the maximum and
#' minimum search range for the amplitude A of cycle of isotope concentrations
#' in precipitation (or streamflow).
#' @param phi A numeric vector, representing the maximum and minimum search range
#' for the phase shift of isotope concentrations in precipitation (or
#' streamflow). The range must be within [0, 2*pi].
#' @param k A numeric vector, representing the maximum and minimum search range
#' for the phase shift of the constant value in the fitted function.
#' @param t A date (or character) vector in format format yyyy-mm-dd
#' (e.g., 2020-12-31)representing the date of observed isotope concentrations
#' @details Fitting non linear sine function of the following form: \cr
#' y = a * sin(2 * pi * t - phi) + k
#' @return A list of object
#' @seealso
#' toDecimal, randomLHS
#' @examples
#' fitSineS <- fitSineNL(obsC = isotopeS$O18,
#'                       a = c(0,10),
#'                       phi = c(0, 2*pi),
#'                       k = c(-20,0),
#'                       t = isotopeS$date,
#'                       nIter = 5000,
#'                       nBestIter = 2,
#'                       weights = isotopeS$streamflow_mm)
#' @export

 fitSineNL <- function(obsC = NULL,
                       a = c(0,10),
                       phi = c(0, 2*pi),
                       k = c(-20,0),
                       t = NULL,
                       nIter = 50000,
                       nBestIter = 30,
                       weights = rep(1, length(obsC))){

   # Output as list object
   output <- list()

   # Check if date input is supplied as date (yyyy-mm-dd)
   if (is.null(t)){
     stop("t must be supplied")
   } else {
     if(!all(!is.na(as.Date(t, format = "%Y-%m-%d")))){
       stop("t must be in format: yyyy-mm-dd")
     }
   }

   # Generate random parameter set using uniform latin hypercube sampling
   parameterSet <- lhs :: randomLHS(nIter,3)
   parameterSet <- as.data.frame(parameterSet)
   colnames(parameterSet) <- c("a", "phi", "k")

   # Convert to user-defined range
   parameterSet$a <- a[1] +  parameterSet$a * (a[2] - a[1])
   parameterSet$phi <- phi[1] +  parameterSet$phi * (phi[2] - phi[1])
   parameterSet$k <- k[1] +  parameterSet$k * (k[2] - k[1])

   # Conver t to decimal number range
   t <- as.Date(t, format = "%Y-%m-%d")
   tDecimal <- lubridate :: decimal_date(t)
   tDecimal <- tDecimal - trunc(tDecimal)

   # Initial the weighted least square value
   weightedLeastSquare <- c()

   # Loop over the parameter set values
   for (i in 1:nIter){

     # Calculate the predict isotope concentration
     cPredicted <- parameterSet$a[i] * sin(2*pi*tDecimal - parameterSet$phi[i]) + parameterSet$k[i]

     # Calculate sum or square error (with weights)
     weightedLeastSquare <- c(weightedLeastSquare,
                              sum(weights * (cPredicted - obsC)^2))
   }

   # Select best saveBestIter models
   sortDecreasing <- order(weightedLeastSquare, decreasing = FALSE)
   parameterSet <- parameterSet[sortDecreasing, ]

   # Return the best nBestIter
   predictedC <- c()
   simulation <- c()
   r2 <- c()

   for (i in 1:nBestIter){
     predictedVariable <- parameterSet$a[i] * sin(2*pi*tDecimal - parameterSet$phi[i]) + parameterSet$k[i]
     predictedC <- c(predictedC, predictedVariable)
     simulation <- c(simulation, rep(paste0("simulation_",i), length(obsC)))
     r2 <- c(r2, cor(x = obsC, y =  predictedVariable)^2)
   }

   output$predictedC <- data.frame(date = rep(t,nBestIter),predictedC = predictedC,
                                   simulation = simulation)
   output$performance <- data.frame(simulation = unique(simulation), weightedLeastSquare =
                                      weightedLeastSquare[sortDecreasing[1:nBestIter]],
                                    r2 = r2)
   output$parameter <- data.frame(cbind(data.frame(simulation = c(1:nBestIter)),
                                        parameterSet[1:nBestIter, ]))
   output$observed <- data.frame(date = t, obsC = obsC)

   # Return output
   return(output)

 }
