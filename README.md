# Fyw
[![R-CMD-check](https://github.com/tamnva/Fyw/workflows/R-CMD-check/badge.svg)](https://github.com/tamnva/Fyw/actions)

## 1. Overview

- Fyw package provides several functions for finding the youngwater fraction and its related parameters (e.g., alpha, beta of the gamma distribution, age threshold of the youngwater fraction)

## 2. Installation

``` r

# Install devtools package
install.packages("devtools")
library(devtools)

# If the package is in use, detach it before installing
detach("package:Fyw", unload = TRUE)

# Install Fyw package for github
install_github("tamnva/Fyw")
```

### 3. Example dataset (included in Fyw packaage)
These data are isotope (O18) concentrations in precipitation and streamflow

``` r
# Load Fyw package
library(Fyw)

# Example dataset
?isotopeP
?isotopeS

# Other packages to run the example Rscript 
library(ggplot2)
library(lubridate)
```

## 4. Approach 1: Fyw estimated using non-linear sine wave fitting
The fitted sine wave is: $c = a \cdot sin(2 \cdot \pi \cdot t - phi) + k$ where $a$, $phi$, and $k$ are parameters need to be estimated. In this approach, user need to be checked if there is outlier in the observed data and remove beforehand.


### 4.1. Fit observed O18 in precipitation to sine wave function

``` r
fitSineP <- fitSineNL(obsC = isotopeP$O18, a = c(0,10), phi = c(0, 2*pi),
                      k = c(-20,0), t = isotopeP$date, nIter = 50000,
                      nBestIter = 10, weights = isotopeP$precippitation_mm)
                      
# remove 'weights = isotopeP$precippitation_mm' for unweighted Fyw

# Plot observed isotope in precipitation and the fitted sine wave
ggplot(fitSineP$predictedC)+
  geom_line(aes(x = date, y = predictedC, color = simulation))+
  geom_point(data = fitSineP$observed, aes(x = date, y = obsC), size = 0.75)+
  scale_y_continuous(limits = c(-30,10))+ylab("precipitation O18 concentration")+
  theme(legend.position = "none")
  
# Fitted parameters (top nBestIter)
fitSineP$parameter
```

### 4.2. Fit observed O18 in streamflow to sine wave function

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

### 4.3. Youngwater fraction (weighted) associated values
Note: To get the Fyw (unweighted) just remove the ```weights``` variables while fitting fitSineP and fitSineS as mentioned in the code above. Now we have different fitted parameters while fitting the sine wave function to isotope concentrations in  precipitation (P) and streamflow (S). We can have different combinations of these fitted parameters to have different Fyw values. This can be done with the ```Fyw``` function.

``` r
# Get youngwater fraction
Fyw <- findFyw(AP = fitSineP$parameter$a, phiP = fitSineP$parameter$phi,
               AS = fitSineS$parameter$a, phiS = fitSineS$parameter$phi)

# Example of finding the age threshold
tauyw <- findtauyw(alpha = Fyw$alpha, beta = Fyw$beta, Fyw = Fyw$Fyw)

# Plot results
ggplot(stack(tauyw), aes(y = values, fill=ind)) +
  geom_boxplot(alpha = 0.5)+
  facet_wrap(~ind, scales="free_y", ncol = 4)+
  xlab("") + ylab("") + 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())
```
Note: The values ```Fyw$P[i]``` and ```Fyw$P[i]``` indicate that ```Fyw$Fyw[i]``` and their associates values (e.g, beta, alpha,...) are estimated from the input combination ```AP[Fyw$P[i]]``` and ```phiP[Fyw$P[i]]``` with ```AS[Fyw$P[i]]``` and ```phiS[Fyw$P[i]]```, which are also given in the ```Fyw``` data frame. 

## 5. Approach 2: Fyw estimated using linear sine wave fitting 
The fitted sine wave is: $c = a \cdot cos(2 \cdot \pi \cdot t) + b \cdot sin(2 \cdot \pi \cdot t) + k$ where $a$, $b$, and $k$ are parameters need to be estimated. In this approach, we use the Iteratively reweighted least squares (IRLS) approach, which can help in limiting the influence of outiler (Kirchner, 2016).

### 5.1. Fit observed O18 in precipitation to sine wave function

``` r
# Convert data to decimal date
tP <-  decimal_date(isotopeP$date)
tP <- tP - trunc(tP)
tS <-  decimal_date(isotopeS$date)
tS <- tS - trunc(tS)

# Fit to sine wave function
fitSinePre <- IRLS(Y = isotopeP$O18,
                   X = data.frame(cos = cos(2*pi*tP),  sin = sin(2*pi*tP)),
                   pweights = isotopeP$precippitation_mm)

# NOTE: remove 'pweights = isotopeP$precippitation_mm' for unweighted Fyw

# Plot fitted sine wave
ggplot()+
  geom_line(aes(x = isotopeP$date, y = as.numeric(fitSinePre$fitted.values), color = "Fitted line"))+
  geom_point(data = isotopeP, aes(x = date, y = O18, color = "observed"))+
  scale_color_manual(values = c("Fitted line" = "red","observed" = "black"))+
  labs(x = "", y = "O18 concentrations in precipitation", color = "")+
  theme(legend.position = "top")
  
# Get amplitude and phase shift
aP <- fitSinePre$coefficients[2]
bP <- fitSinePre$coefficients[3]
AP <- sqrt(sum(aP^2 + bP^2))
phiP <- atan(bP/aP)
```

### 5.2. Fit observed O18 in precipitation to sine wave function

``` r
# Fit to sine wave function
fitSineStr <- IRLS(Y = isotopeS$O18,
                   X = data.frame(cos = cos(2*pi*tS),  sin = sin(2*pi*tS)),
                   pweights = isotopeS$streamflow_mm)

# NOTE: remove 'pweights = isotopeS$streamflow_mm' for unweighted Fyw

# Plot fitted sine wave
ggplot()+
  geom_line(aes(x = isotopeS$date, y = as.numeric(fitSineStr$fitted.values), color = "Fitted line"))+
  geom_point(data = isotopeS, aes(x = date, y = O18, color = "observed"))+
  scale_color_manual(values = c("Fitted line" = "red","observed" = "black"))+
  labs(x = "", y = "O18 concentrations in streamflow", color = "")+
  theme(legend.position = "top")
  
  
# Get amplitude and phase shift
aS <- fitSineStr$coefficients[2]
bS <- fitSineStr$coefficients[3]
AS <- sqrt(sum(aS^2 + bS^2))
phiS <- atan(bS/aS)
```
### 5.3. Find the youngwater fraction and associated values

``` r
# Youngwater fraction (weighted)
Fyw_new <- findFyw(AP = AP, phiP = phiP,
                   AS = AS, phiS = phiS)

# Age theshold of Fyw_new
tauyw_new <- findtauyw(alpha = Fyw_new$alpha, beta = Fyw_new$beta,
                       Fyw = Fyw_new$Fyw)

# Uncertainty estimation
# TODO?


```



