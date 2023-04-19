## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
# Load require packages
library(Fyw)
library(lubridate)
#library(ggplot2)  # uncomment to see plot

## ---- message=FALSE-----------------------------------------------------------
# Get isotope in precipitation (P) and streamflow (S) from the the Alp catchment
isotopeP <- subset(isotopeData, catchment == "Alp" & variable == "precipitation")
isotopeS <- subset(isotopeData, catchment == "Alp" &  variable == "streamflow")

## ---- message=FALSE-----------------------------------------------------------
#------------------------------------------------------------------------------#
#             Fit sine wave to observed isotope in precipitation               #
#------------------------------------------------------------------------------#
# First, convert date to decimal, ranging from [0, 1]
t <- decimal_date(isotopeP$date) - trunc(decimal_date(isotopeP$date))
  
# Fit to sine wave to observed isotope in P using IRLS method
sineP <- IRLS(Y = isotopeP$delta_18O,
              X = data.frame(cos = cos(2*pi*t), sin = sin(2*pi*t)),
              pweights = isotopeP$water_flux_mm)

# Uncomment to see the plot
# ggplot()+
#   geom_point(data = isotopeP, aes(x = date, y = delta_18O, color = "cP (observed)"), alpha = 0.5) + 
#   geom_line(aes(x = isotopeP$date, y = as.numeric(sineP$fitted.values), color = "cP (fitted sine-wave)"), linewidth = 0.75)

#------------------------------------------------------------------------------#
#             Fit sine wave to observed isotope in streamflow                  #
#------------------------------------------------------------------------------#
# First, convert date to decimal, ranging from [0, 1]
t <- decimal_date(isotopeS$date) - trunc(decimal_date(isotopeS$date))
  
# Fit to sine wave to observed isotope in P using IRLS method
sineS <- IRLS(Y = isotopeS$delta_18O,
              X = data.frame(cos = cos(2*pi*t), sin = sin(2*pi*t)),
              pweights = isotopeS$water_flux_mm)

# Uncomment to see the plot
# ggplot()+
#   geom_point(data = isotopeS, aes(x = date, y = delta_18O, color = "cS (observed)"), alpha = 0.5) + 
#   geom_line(aes(x = isotopeS$date, y = as.numeric(sineS$fitted.values), color = "cS (fitted sine-wave)"), linewidth = 0.75)

## ---- message=FALSE-----------------------------------------------------------
# Parameters of the fitted sine wave to isotope in precipitation
kP <- as.numeric(sineP$coefficients[1])
aP <- as.numeric(sineP$coefficients[2])
bP <- as.numeric(sineP$coefficients[3])

# Parameters of the fitted sine wave to isotope in precipitation
kS <- as.numeric(sineS$coefficients[1])
aS <- as.numeric(sineS$coefficients[2])
bS <- as.numeric(sineS$coefficients[3])

# Amplitudes of the fitted sine wave to precipitation and streamflow (please see Eq. 6; Kirchner, (2016))
AP <- sqrt(aP^2 + bP^2)
AS <- sqrt(aS^2 + bS^2)

# Phase 
phiP <- findPhi(aP = aP, bP = bP)
phiS <- findPhi(aP = aS, bP = bS)

# Phase shift (please see Eq. 11; Kirchner, (2016)))
phiS_phiP <- phiS - phiP

# Parameters of the gamma (transit time distribution) (please see Eqs. 10, 11; Kirchner, (2016))
alphaBetaMTT <- findAlphaBetaMTT(phiS_phiP = phiS_phiP,
                                 AS = AS,
                                 AP = AP)

# Age threshold of the young water fraction (please see Eq. 14; Kirchner, (2016)))
tauyw <- findtauyw(alpha = alphaBetaMTT$alpha, method = "approximation")$tauyw

# Young water fraction (equals to the amplitude ratio)
Fyw_1 <- AS/AP

# Or the young water fraction can be calculated using the lower incomplete gamma function (please see Eq. 13; Kirchner, (2016))
Fyw_2 <- pgamma(q = tauyw,
                shape = alphaBetaMTT$alpha,
                scale = alphaBetaMTT$beta,
                lower.tail = TRUE)

# Show the results in data frame
t(data.frame(AP = AP, phiP = phiP, kP = kP, AS = AS, phiS = phiS, kS = kS, 
             Fyw_1 = Fyw_1, Fyw_2 = Fyw_2, tauyw = tauyw))

