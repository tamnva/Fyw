# Fyw

## Overview

- Fyw package provides several functions for finding the youngwater fraction and its related parameters (e.g., alpha, beta of the gamma distribution, age threshold of the youngwater fraction)

## Installation

``` r

# Install devtools package
install.packages("devtools")
library(devtools)

# If the package is in use, detach it before installing
detach("package:Fyw", unload = TRUE)

# Install Fyw package for github
install_github("tamnva/Fyw")
library(Fyw)
```

## Usage

1. Example data set

``` r
?isotopeP
?isotopeQ
```
2. Fit observed O18 in precipitation to sine wave function in form (weighted with precipatation volume): C = A * sin(2*pi*t - phi) + k

``` r
fitSineP <- fitSineNL(obsC = isotopeP$O18, a = c(0,10), phi = c(0, 2*pi),
                      k = c(-20,0), t = isotopeP$date, nIter = 50000,
                      nBestIter = 10, weights = isotopeP$precippitation_mm)

# Plot observed isotope in precipitation and the fitted sine wave
library(ggplot2)
ggplot(fitSineP$predictedC)+
  geom_line(aes(x = date, y = predictedC, color = simulation))+
  geom_point(data = fitSineP$observed, aes(x = date, y = obsC))+
  scale_y_continuous(limits = c(-30,10))+
  theme(legend.position = "none")
  
# Fitted parameters (top nBestIter)
fitSineP$parameter
```

3. Fit observed O18 in streamflow to sine wave function in form (flow weighted): C = A * sin(2*pi*t - phi) + k

``` r
fitSineQ <- fitSineNL(obsC = isotopeQ$O18, a = c(0,10), phi = c(0, 2*pi),
                      k = c(-20,0), t = isotopeQ$date, nIter = 50000,
                      nBestIter = 10, weights = isotopeQ$precippitation_mm)

ggplot(fitSineQ$predictedC)+
  geom_line(aes(x = date, y = predictedC, color = simulation))+
  geom_point(data = fitSineQ$observed, aes(x = date, y = obsC))+
  scale_y_continuous(limits = c(-30,10))+
  theme(legend.position = "none")
  
# Fitted parameters (top nBestIter)
fitSineQ$parameter
```

4. Youngwater fraction
``` r
fitSineQ$parameter$a/fitSineP$parameter$a
```
