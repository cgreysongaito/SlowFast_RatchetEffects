rm(list=ls())

library(QPot)
library(zeallot)
library(tidyverse)
library(viridis)

setwd("/home/chrisgg/Documents/Guelph/PhD/TimeDelays/")

## Set up isoclines
resiso_data_05 <- read_csv("data/resiso_data_05.csv")
resiso_data_07 <- read_csv("data/resiso_data_07.csv")

#Model Set Up - NOTE g and f parameters added to allow for differential stochastic noise
res_model <- "( r * x * (1 - ( ( x * g ) / k ) ) ) - ( a * x * y * f ) / ( 1 +  ( a * h * x * g ) ) "
con_model <- "p * ( ( e * a * x * g * y ) / (1 + ( a * h * x * g ) ) - ( m * y ) )"

#General Functions
param_model <- function(param) {
  c(eff, ep, res_g, con_f) %<-% param
  model_parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = eff, m = 0.4, p = ep, g = res_g, f = con_f)
  res_model_param <- Model2String(res_model, parms = model_parms, supress.print = TRUE)
  con_model_param <- Model2String(con_model, parms = model_parms, supress.print = TRUE)
  parms_eqn <- list(res_model_param, con_model_param)
  return(parms_eqn)
}

r_convert <- 1/sqrt(17) #1:4 different noise
c_convert <- 4/sqrt(17) 
# Set bounds and step numbers for QPotential calculation
bounds.x <- c(0, 4/r_convert)
bounds.y <- c(0, 4/c_convert)
step.number.x <- 1000
step.number.y <- 1000

## Section for manuscript
#Quasi-potential with different noise (1:4) for epsilon of 1/12.7 & efficiency of 0.5
param_symnoise_ep127_eff05 <- list(eff = 0.5,
                                   ep = 1/12.7,
                                   res_g = r_convert,
                                   con_f = c_convert)

u0_symnoise_ep127_eff05 <- c(x = 2.02/r_convert, y = 1.65/c_convert)

qp_symnoise_ep127_eff05 <- QPotential(x.rhs = param_model(param_symnoise_ep127_eff05)[[1]],
                                      x.start = u0_symnoise_ep127_eff05[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep127_eff05)[[2]],
                                      y.start = u0_symnoise_ep127_eff05[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)
par(mfrow=c(1,1), mgp = c(2.6,1,0))
png("figs/qp_symnoise_ep127_eff05.png")
filled.contour(x = seq(0,4, length.out = nrow(qp_symnoise_ep127_eff05)),
               y = seq(0,4, length.out = nrow(qp_symnoise_ep127_eff05)),
               qp_symnoise_ep127_eff05,
               xlim = c(0,3),
               ylim = c(0,3),
               cex.lab = 2,
               color.palette = viridis,
               nlevels = 8,
               xlab = "Resource",
               ylab = "Consumer",
               plot.axes = { axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(resiso_data_05$xrange, resiso_data_05$resiso, col = "blue", lwd = 3);
                 axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(c(2.02,2.02), c(0.0, 3.0), col = "orange", lwd = 3);
                 })
dev.off()

#Quasi-potential with different noise (1:4) for epsilon of 1/12.7 & efficiency of 0.7
param_symnoise_ep127_eff07 <- list(eff = 0.7,
                                   ep = 1/12.7,
                                   res_g = r_convert,
                                   con_f = c_convert)

u0_symnoise_ep127_eff07 <- c(x = 0.96/r_convert, y = 2.28/c_convert)

qp_symnoise_ep127_eff07 <- QPotential(x.rhs = param_model(param_symnoise_ep127_eff07)[[1]],
                                      x.start = u0_symnoise_ep127_eff07[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep127_eff07)[[2]],
                                      y.start = u0_symnoise_ep127_eff07[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)


filled.contour(x = seq(0,4, length.out = nrow(qp_symnoise_ep127_eff07)),
               y = seq(0,4, length.out = nrow(qp_symnoise_ep127_eff07)),
               qp_symnoise_ep127_eff07,
               xlim = c(0,3),
               ylim = c(0,3),
               cex.lab = 2,
               color.palette = viridis,
               nlevels = 8,
               xlab = "Resource",
               ylab = "Consumer",
               plot.axes = { axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(resiso_data_07$xrange, resiso_data_07$resiso, col = "blue", lwd = 3);
                 axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(c(0.96, 0.96), c(0.0, 3.0), col = "orange", lwd = 3)})


#Quasi-potential with different noise (1:4) for epsilon of 1/250 & efficiency of 0.5
param_symnoise_ep0004_eff05 <- list(eff = 0.5,
                                   ep = 1/250,
                                   res_g = r_convert,
                                   con_f = c_convert)

u0_symnoise_ep0004_eff05 <- c(x = 2.02/r_convert, y = 1.65/c_convert)

qp_symnoise_ep0004_eff05 <- QPotential(x.rhs = param_model(param_symnoise_ep0004_eff05)[[1]],
                                      x.start = u0_symnoise_ep0004_eff05[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep0004_eff05)[[2]],
                                      y.start = u0_symnoise_ep0004_eff05[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

filled.contour(x = seq(0,4, length.out = nrow(qp_symnoise_ep0004_eff05)),
               y = seq(0,4, length.out = nrow(qp_symnoise_ep0004_eff05)),
               qp_symnoise_ep0004_eff05,
               xlim = c(0,3),
               ylim = c(0,3),
               cex.lab = 2,
               color.palette = viridis,
               nlevels = 8,
               xlab = "Resource",
               ylab = "Consumer",
               plot.axes = { axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(resiso_data_05$xrange, resiso_data_05$resiso, col = "blue", lwd = 3);
                 axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(c(2.02, 2.02), c(0.0, 3.0), col = "orange", lwd = 3)})

#Quasi-potential with different noise (1:4) for epsilon of 1/250 & efficiency of 0.7
param_symnoise_ep0004_eff07 <- list(eff = 0.7,
                                    ep = 1/250,
                                    res_g = r_convert,
                                    con_f = c_convert)

u0_symnoise_ep0004_eff07 <- c(x = 0.96/r_convert, y = 2.28/c_convert)

qp_symnoise_ep0004_eff07 <- QPotential(x.rhs = param_model(param_symnoise_ep0004_eff07)[[1]],
                                      x.start = u0_symnoise_ep0004_eff07[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep0004_eff07)[[2]],
                                      y.start = u0_symnoise_ep0004_eff07[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

filled.contour(x = seq(0,4, length.out = nrow(qp_symnoise_ep0004_eff07)),
               y = seq(0,4, length.out = nrow(qp_symnoise_ep0004_eff07)),
               qp_symnoise_ep0004_eff07,
               xlim = c(0,3),
               ylim = c(0,3),
               cex.lab = 2,
               color.palette = viridis,
               nlevels = 8,
               xlab = "Resource",
               ylab = "Consumer",
               plot.axes = { axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(resiso_data_05$xrange, resiso_data_05$resiso, col = "blue", lwd = 3);
                 axis(1, cex.axis=1.5); axis(2, cex.axis=1.5); lines(c(0.96, 0.96), c(0.0, 3.0), col = "orange", lwd = 3)})
