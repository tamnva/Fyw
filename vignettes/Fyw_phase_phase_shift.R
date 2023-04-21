## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, eval = TRUE----------------------------------------
# Load require packages
library(ggplot2)

## ---- message=FALSE, eval = TRUE, fig.width=7, fig.height=3-------------------
# Generate time vector from 0 to 1
# t = [0,1] corresponding to year [0,1] and phase [0,2pi]
t <- seq(0,2,0.05)

# Assume that the phase of precipitation and streamflow is pi/2 and pi, repsectively
phiP <- pi
phiS <- 2*pi + pi/2

# The phase shift in this case is
phiS - phiP

# Assume that the amplitude of precipitation is 2 and of streamflow is 1.5
AP <- 2
AS <- 1.5

# Assume that kP = kS = -10
kP = -10
kS = -10

# Now get the sine wave of precipitation and streamflow
cP <- AP * sin (2*pi * t - phiP) + kP
cS <- AS * sin (2*pi * t - phiS) + kS

# Assemble the data for plot
df <- data.frame(t = t, cP = cP, cS = cS)
ps <- data.frame(t = c(0.5, 1.25), phaseShift = c(-10,-10))

# Plot to see the phase shift
ggplot(df) +
  geom_line(aes(x = t, y = cP, color = "cP"), linewidth = 0.75) +
  geom_line(aes(x = t, y = cS, color = "cS"), linewidth = 0.75) +
  geom_line(data = ps, aes(x = t, y = phaseShift, color = "phase shift"), linewidth = 2) +
  scale_x_continuous("time (year)", 
                     sec.axis = sec_axis(~. *2*pi, expression(paste(varphi, " (rad)" )), 
                                         breaks = seq(0,2,0.5)*2*pi,
                                         labels = c(0, 
                                                    expression(paste(pi)), 
                                                    expression(paste("2",pi)), 
                                                    expression(paste("3",pi)), 
                                                    expression(paste("4",pi)))
                     )
  ) + 
  scale_color_manual(labels = c(expression(paste("   ", c[P], " (",varphi[P], " = ", pi ,")")), 
                                expression(paste("   ", c[S], " (",varphi[S], " = 2.5", pi ,")")), 
                                expression(paste("   ", varphi[S], " - ", varphi[P], " = 1.5", pi))),
                     values = c("red", "blue", "green")) +
  labs(y = expression(paste(c[P], "  and  ", c[S])), color = "") +
  theme_bw()

