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
?isotopeS
```
2. Fit observed O18 in precipitation to sine wave function in form (weighted with precipatation volume): $c = a \cdot sin(2 \cdot \pi \cdot t - phi) + k$

``` r
fitSineP <- fitSineNL(obsC = isotopeP$O18, a = c(0,10), phi = c(0, 2*pi),
                      k = c(-20,0), t = isotopeP$date, nIter = 50000,
                      nBestIter = 10, weights = isotopeP$precippitation_mm)
                      
# remove 'weights = isotopeP$precippitation_mm' for unweighted Fyw

# Plot observed isotope in precipitation and the fitted sine wave
library(ggplot2)
ggplot(fitSineP$predictedC)+
  geom_line(aes(x = date, y = predictedC, color = simulation))+
  geom_point(data = fitSineP$observed, aes(x = date, y = obsC), size = 0.75)+
  scale_y_continuous(limits = c(-30,10))+ylab("precipitation O18 concentration")+
  theme(legend.position = "none")
  
# Fitted parameters (top nBestIter)
fitSineP$parameter
```

3. Fit observed O18 in streamflow to sine wave function in form (flow weighted): $c = a \cdot sin(2 \cdot \pi \cdot t - phi) + k$

``` r
fitSineS <- fitSineNL(obsC = isotopeS$O18, a = c(0,10), phi = c(0, 2*pi),
                      k = c(-20,0), t = isotopeS$date, nIter = 50000,
                      nBestIter = 10, weights = isotopeS$streamflow_mm)

# remove 'weights = isotopeS$streamflow_mm' for unweighted Fyw

ggplot(fitSineS$predictedC)+
  geom_line(aes(x = date, y = predictedC, color = simulation))+
  geom_point(data = fitSineS$observed, aes(x = date, y = obsC), size = 0.75)+
  scale_y_continuous(limits = c(-30,10))+ ylab("instream O18 concentration")+
  theme(legend.position = "none")
  
# Fitted parameters (top nBestIter)
fitSineS$parameter
```

4. Youngwater fraction (weighted). To get the Fyw (unweighted) just remove the ```weights``` variables while fitting fitSineP and fitSineS as mentioned in the code above. Now we have different fitted parameters while fitting the sine wave function to isotope concentrations in  precipitation (P) and streamflow (S). We can have different combinations of these fitted parameters to have different Fyw values. This can be done with the ```Fyw``` function.
``` r
Fyw_weighted <- Fyw(AP = fitSineP$parameter$a, phiP = fitSineP$parameter$phi,
                    AS = fitSineS$parameter$a, phiS = fitSineS$parameter$phi)

```
Note: The values ```Fyw_weighted$P[i]``` and ```Fyw_weighted$P[i]``` indicate that ```Fyw_weighted$Fyw[i]``` and their associates values (e.g, beta, alpha,...) are estimated from the input combination ```AP[Fyw_weighted$P[1]]``` and ```phiP[Fyw_weighted$P[1]]``` with ```AS[Fyw_weighted$P[1]]``` and ```phiS[Fyw_weighted$P[1]]```, which are also given in the ```Fyw_weighted``` data frame. 
