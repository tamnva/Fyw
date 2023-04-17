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

1. Kirchner, J. W. (2016). Aggregation in environmental systems – Part 1: Seasonal tracer cycles quantify young water fractions, but not mean transit times, in spatially heterogeneous catchments. *Hydrol. Earth Syst. Sci.*, 20, 279–297. https://doi.org/10.5194/hess-20-279-2016.

2. von Freyberg, J., Rücker, A., Zappa, M., Schlumpf, A., Studer, B., Kirchner, J. W.  (2022). Four years of daily stable water isotope data in stream water and precipitation from three Swiss catchments. *Sci. Data*, 9(46). https://doi.org/10.1038/s41597-022-01148-1.

3. von Freyberg, J., Allen, S. T., Seeger, S., Weiler, M., and Kirchner, J. W. (2018). Sensitivity of young water fractions to hydro-climatic forcing and landscape properties across 22 Swiss catchments. *Hydrol. Earth Syst. Sci.*, 22,
   3841–3861. https://doi.org/10.5194/hess-22-3841-2018-supplement.

