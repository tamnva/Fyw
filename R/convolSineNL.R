

convolSineNL <- function(sDate = NULL, eDate = NULL, AP = NULL,
                         phiP = NULL, kP = NULL, shape = NULL,
                         scale = NULL, isotopeS = NULL, weights = NULL){

  # Sine function in form c = A sin (2*pi*t - phi) + k
  SineNL <- function(t, A, phi, k){
    A * sin(2*pi*t - phi) + k
  }

  tRange <- seq(from = as.Date("2010-01-01", format = "%Y-%m-%d"),
                to = as.Date("2020-12-31", format = "%Y-%m-%d"),
                by = "days")
  tRangeAbs <- decimal_date(tRange)
  tRangeDec <- tRangeAbs - trunc(tRangeAbs)
  tSRangeAbs <- decimal_date(isotopeS$date)
  tSLoc <- which(tRangeAbs %in% tSRangeAbs)
  tRangeFromOrigin <- tRangeAbs- tRangeAbs[1]
  ntimeStep <- length(tRangeAbs)

  # Get cummulative concentration using cummulative gamma
  isoC <- SineNL(tRangeDec, A = 1.909, phi = 2.016, k = -11.05995)
  plot(x = tRangeAbs, y = isoC, type = "l")
  convoluteC <- rep(0,ntimeStep)

  for (i in 1:ntimeStep){
    simC <- isoC[i] * pgamma(tRangeFromOrigin, shape = 0.802,
                             scale = 0.127, lower = TRUE)
    simC <- diff(c(0,simC))
    convoluteC[i:ntimeStep] <- convoluteC[i:ntimeStep] + simC[1:(ntimeStep - i + 1)]
  }

  #Calculate mean square
  plot(x = tRangeAbs, y = convoluteC, type = "l", ylim = c(-20,-5))
  lines(x = tRangeAbs[tSLoc], y = isotopeS$O18)


  # Young water fraction 0.2 years
  pgamma(0.2, shape = 0.802, scale = 0.127, lower = TRUE)

}
