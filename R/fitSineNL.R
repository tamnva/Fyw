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
#' y = a \cdot sin(2 \cdot pi \cdot t - phi) + k
#' @return A list of object
#' @seealso
#' toDecimal, randomLHS
#' @examples
#' fitSineNonLinear(y = runif(10), t = runif(10))

 fitSineNL <- function(obsC = NULL, a = c(0,5), phi = c(0, 2*pi), k = c(-20,-10),
                       t = NULL, nIter = 10000, nBestIter = 100,
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
   t <- lubridate :: decimal_date(as.Date(t, format = "%Y-%m-%d"))
   t <- t - trunc(t)

   # Initial the weighted least square value
   weightedLeastSquare <- c()

   # Loop over the parameter set values
   for (i in 1:nIter){

     # Calculate the predict isotope concentration
     cPredicted <- parameterSet$a[i] * sin(2*pi*t - parameterSet$phi[i])
     + parameterSet$k[i]

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
   performance <- c()

   for (i in 1:nBestIter){
     predictedVariable <- parameterSet$a[i] * sin(2*pi*t - parameterSet$phi[i]) + parameterSet$k[i]
     predictedC <- c(predictedC, predictedVariable)
     simulation <- c(simulation, rep(paste0("simulation_",i)), length(c))
     weightedLeastSquare <-  weightedLeastSquare[sortDecreasing[i]]
     r2 <- cor(x = obsC, y =  predictedVariable)^2

     if (i == 1) {
       performance <- data.frame(simulation = paste0("simulation_",i),#
                                weightedLeastSquare = weightedLeastSquare,
                                r2 = r2)
     } else {
       performance <- rbind(performance, c(paste0("simulation_",i),weightedLeastSquare, r2))
     }

   }

   output$predictedC <- tibble::tibble(date = rep(t,nBestIter),
                                       predictedC = predictedC,
                                       simulation = simulation)
   output$performance <- tibble::as_tibble(performance)
   output$parameter <- tibble::as.tibble(cbind(data.frame(simulation = c(1:nBestIter)),
                                               parameterSet[1:nBestIter, ]))

   # Return output
   return(output)

 }
