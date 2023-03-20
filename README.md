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

``` r
# Observed isotope (O18) concentrations in rainfall - remove missing values

isotopeP <- na.omit(isotopeP)


# Observed isotope (O18) concentrations in streamflow - remove missing values

isotopeQ <- na.omit(isotopeQ)

# Fit observed data to sine wave function in form 
# C = A * sin(2*pi*t - phi) + k

fitSineP <- fitSineNL(obsC = isotopeP$O18, a = c(0,10), phi = c(0, 2*pi), k = c(-20,0), t = isotopeP$date, nIter = 50000, nBestIter = 10, weights = isotopeP$precippitation_mm)

```
