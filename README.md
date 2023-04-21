# Fyw <a href="https://github.com/tamnva/Fyw/blob/master/vignettes/icon.svg"><img src="vignettes/icon.svg" align="right" height="120" /></a>

[![R-CMD-check](https://github.com/tamnva/Fyw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tamnva/Fyw/actions/workflows/R-CMD-check.yaml)[![DOI](https://zenodo.org/badge/615738927.svg)](https://zenodo.org/badge/latestdoi/615738927)


## 1. Overview

- Fyw package provides several functions for finding the youngwater fraction and its related parameters (e.g., alpha, beta of the gamma distribution, age threshold of the young water fraction).

## 2. Installation

``` r

# Install devtools package if needed
if (!require(devtools)) install.packages("devtools")
library(devtools)

# If the package is in use, detach it before installing
# detach("package:Fyw", unload = TRUE)

# Install Fyw package for github
install_github("tamnva/Fyw", build_vignettes = TRUE)
```

### 3. Quick start!

This section describes how the young water fraction (Fyw) is calculated based on the traditional procedure (Kirchner, 2016) and revised procedure (this package). 

First we need to load require package (or install if you have not installed these packages)

```r
library(Fyw)         # This package
library(ggplot2)     # For plotting the results
library(tibble)      # For displaying the results in tibble format
library(lubridate)   # Convert date to decimal
```

Then, we need data for this demonstration. Let's take the isotope data from the Alp catchment.

```r
# Get isotope in precipitation (P) and streamflow (S) from the the Alp catchment
isotopeP <- subset(isotopeData, catchment == "Alp" & variable == "precipitation")
isotopeS <- subset(isotopeData, catchment == "Alp" &  variable == "streamflow")
```

#### 3.1. Fyw estimation using the traditional procedure
In the traditional procedure (Kirchner, 2016), Fyw is estimated using the following steps:

-   Step 1: Fit the sine wave to observed isotope in precipitation and streamflow using the Iteratively Least Squares (IRLS) regression approach. The fitted sine wave has the following form (please see Eq. 5; Kirchner, (2016)):

$$ c_P = a_P \cdot cos(2 \pi t) + b_P \cdot sin(2 \pi t) + k_P $$

```r
#------------------------------------------------------------------------------#
#             Fit sine wave to observed isotope in precipitation               #
#------------------------------------------------------------------------------#
# First, convert date to decimal, ranging from [0, 1]
t <- decimal_date(isotopeP$date) - trunc(decimal_date(isotopeP$date))
  
# Fit to sine wave to observed isotope in P using IRLS method
sineP <- IRLS(Y = isotopeP$delta_18O,
              X = data.frame(cos = cos(2*pi*t), sin = sin(2*pi*t)),
              pweights = isotopeP$water_flux_mm)

# Visualize results (fitted sine wave to isotope in precipitation)
ggplot()+
  geom_point(data = isotopeP, aes(x = date, y = delta_18O, color = "cP (observed)"), alpha = 0.5) + 
  geom_line(aes(x = isotopeP$date, y = as.numeric(sineP$fitted.values), color = "cP (fitted sine-wave)"), linewidth = 0.75)

#------------------------------------------------------------------------------#
#             Fit sine wave to observed isotope in streamflow                  #
#------------------------------------------------------------------------------#
# First, convert date to decimal, ranging from [0, 1]
t <- decimal_date(isotopeS$date) - trunc(decimal_date(isotopeS$date))
  
# Fit to sine wave to observed isotope in P using IRLS method
sineS <- IRLS(Y = isotopeS$delta_18O,
              X = data.frame(cos = cos(2*pi*t), sin = sin(2*pi*t)),
              pweights = isotopeS$water_flux_mm)
              
# Visualize results (fitted sine wave to isotope in streamflow)
ggplot()+
 geom_point(data = isotopeS, aes(x = date, y = delta_18O, color = "cS (observed)"), alpha = 0.5) + 
 geom_line(aes(x = isotopeS$date, y = as.numeric(sineS$fitted.values), color = "cS (fitted sine-wave)"), linewidth = 0.75)
```

-   Step 2: Find Fyw from the two fitted sine waves

```r
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

phi <- findPhiPhaseShift(aP = aP, bP = bP, aS = aS, bS = bS)
phiP <- phi$phiP
phiS <- phi$phiS
phiS_phiP <- phi$phiS_phiP

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
```

#### 3.2. Fyw estimation using the revised procedure
In the revised procedure, Fyw is estimated using the following steps:

-   **Step 1**: Fit the sine wave to observed isotope in precipitation. The fitted sine 
wave has the following form (please see Eq. 4; Kirchner, (2016)):


$$ c_P = A_P \cdot sin(2 \pi t - \varphi_P)  + k_P $$

```r
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
ggplot(fitSineP$simulated)+
  geom_line(aes(x = date, y = simulated, color = simulation)) +         # 10 best fitted sine waves
  geom_point(data = fitSineP$observed, aes(x = date, y = observed)) +   # observed data
  labs(x = "", y = expression(paste(delta^{18}, "O in precipitation (‰)"))) +
  theme(legend.position = "none")
```

-   Step 2: Find parameters of the gamma distribution (alpha, beta) by convolution of the fitted sine wave to isotope in precipitation (step 1) with the gamma functions (please see Eq. 4; Kirchner, (2016)) and evaluate the goodness-of-fit with observed isotope in streamflow.

```r
# NOTE: This require very high computational demand

# Parameters of the fitted sine wave to isotope in precipitation
fitGamma <- fitGamma(AP = fitSineP$parameter$a,
                     phiP = fitSineP$parameter$phi,
                     kP = fitSineP$parameter$k,
                     alphaRange = c(0.75,0.9),
                     betaRange = c(0.2,0.3),
                     simulatedDate = isotopeS$date,
                     fittedData = isotopeS$delta_18O,
                     weight = isotopeS$water_flux_mm,
                     nIter = 30,
                     nBestIter = 5,
                     nCores = 4,
                     nWarmupYears = 5,
                     setSeed = 2023)

# Plot to see the quality of fit
ggplot()+   
  geom_point(data = isotopeS, aes(x = date, y = delta_18O), alpha = 0.35) +                  # Observed
  geom_line(data = fitGamma$simulated, aes(x = date, y = simulated, color = simulation)) +   # Simulated
  labs(x = "", y = expression(paste(delta^{18}, "O in streamflow (‰)"))) +
  theme(legend.position = "none")
```

-   Step 3: Now calculate the young water fraction with user defined age threshold (please see Eq. 13; Kirchner, (2016))

```r

# Let's say we want to find the Fyw fraction with age threshold of 2 months
tauyw <- 2/12 

# Find Fyw 
Fywater <- c()
for (i in 1:nrow(fitGamma$parameterSet)){
  Fywater <- c(Fywater, pgamma(q = tauyw,
                               shape = fitGamma$parameterSet$alpha[i],
                               scale = fitGamma$parameterSet$beta[i],
                               lower.tail = TRUE))
}

# show Fyw values
Fywater
```

### References

1. Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal tracer cycles quantify young water fractions, but not mean transit times, in spatially heterogeneous catchments. *Hydrol. Earth Syst. Sci.*, 20, 279–297. https://doi.org/10.5194/hess-20-279-2016.

2. von Freyberg, J., Rücker, A., Zappa, M., Schlumpf, A., Studer, B., Kirchner, J. W.  (2022). Four years of daily stable water isotope data in stream water and precipitation from three Swiss catchments. *Sci. Data*, 9(46). https://doi.org/10.1038/s41597-022-01148-1.

3. von Freyberg, J., Allen, S. T., Seeger, S., Weiler, M., and Kirchner, J. W. (2018). Sensitivity of young water fractions to hydro-climatic forcing and landscape properties across 22 Swiss catchments. *Hydrol. Earth Syst. Sci.*, 22, 3841–3861. https://doi.org/10.5194/hess-22-3841-2018-supplement.

