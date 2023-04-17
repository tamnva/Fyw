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

### 3. How to use this package

Please see the package vignettes

### References

Kirchner, J. W. 2016. “Aggregation in Environmental Systems – Part 1: Seasonal 
Tracer Cycles Quantify Young Water Fractions, but Not Mean Transit Times, in 
Spatially Heterogeneous Catchments.” Hydrology and Earth System Sciences 20 (1): 279–97. https://doi.org/10.5194/hess-20-279-2016.
