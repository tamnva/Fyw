## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
# Load require packages
library(Fyw)

# Uncomment if users want to see the plots
# library(ggplot2)   

## ---- message=FALSE-----------------------------------------------------------
# Get isotope in precipitation (P) and streamflow (S) from the the Alp catchment
isotopeP <- subset(isotopeData, catchment == "Alp" & variable == "precipitation")
isotopeS <- subset(isotopeData, catchment == "Alp" &  variable == "streamflow")

## ---- message=FALSE, fig.height=2.5, fig.width=7------------------------------
#------------------------------------------------------------------------------#
#             Fit sine wave to observed isotope in precipitation               #
#------------------------------------------------------------------------------#
fitSineP <- fitSineNL(observed = isotopeP$delta_18O,
                      a = c(1,3),
                      phi = c(0, 2*pi),
                      k = c(-12,-8), 
                      t = isotopeP$date,
                      nIter = 10000,
                      nBestIter = 10,
                      weight = isotopeP$water_flux_mm,
                      setSeed = 2023)

# Plot the results (fitted sine wave and observed)
# ggplot(fitSineP$simulated)+
#  geom_line(aes(x = date, y = simulated, color = simulation)) +         # 10 best fitted sine waves
#  geom_point(data = fitSineP$observed, aes(x = date, y = observed)) +   # observed data
#  labs(x = "", y = expression(paste(delta^{18}, "O in precipitation (‰)"))) +
#  theme(legend.position = "none")

## ---- message=FALSE, message=FALSE, fig.height=2.5, fig.width=7, eval = FALSE----
#  # NOTE: Due to high computational demand, here I set nIter = 4 and nWarmupYears = 5
#  #       but to have a better fit nIter should be large (e.g, 2000),
#  #       and nWarmupYears should be large (e.g., 10) to reduce the impact of initial conditions
#  #       and the alphaRange, betaRange should be larger (e.g., from 0.01 to 5)
#  
#  # Here we randomly combined 10 best fit sine wave functions (Step 1) with 3 randomly
#  #       generated gamma parameters and extract results of the 2 best fitted gamma functions
#  
#  # Parameters of the fitted sine wave to isotope in precipitation
#  fitGamma <- fitGamma(AP = fitSineP$parameter$a,
#                       phiP = fitSineP$parameter$phi,
#                       kP = fitSineP$parameter$k,
#                       alphaRange = c(0.75,0.9),
#                       betaRange = c(0.2,0.3),
#                       simulatedDate = isotopeS$date,
#                       fittedData = isotopeS$delta_18O,
#                       weight = isotopeS$water_flux_mm,
#                       nIter = 3,
#                       nBestIter = 2,
#                       nCores = 1,
#                       nWarmupYears = 5,
#                       setSeed = 2023)
#  
#  # Plot to see the quality of fit
#  # ggplot()+
#  #   geom_point(data = isotopeS, aes(x = date, y = delta_18O), alpha = 0.35) +                  # Observed
#  #  geom_line(data = fitGamma$simulated, aes(x = date, y = simulated, color = simulation)) +   # Simulated
#  #  labs(x = "", y = expression(paste(delta^{18}, "O in streamflow (‰)"))) +
#  #  theme(legend.position = "none")

## ---- message=FALSE, message=FALSE, fig.height=2.5, fig.width=7, eval = FALSE----
#  
#  # Let's say we want to find the Fyw fraction with age threshold of 2 months
#  tauyw <- 2/12
#  
#  # Find Fyw
#  Fywater <- c()
#  for (i in 1:nrow(fitGamma$parameterSet)){
#    Fywater <- c(Fywater, pgamma(q = tauyw,
#                                 shape = fitGamma$parameterSet$alpha[i],
#                                 scale = fitGamma$parameterSet$beta[i],
#                                 lower.tail = TRUE))
#  }
#  
#  # show Fyw values
#  Fywater
#  

