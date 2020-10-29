rm(list=ls())

library(QPot)
library(zeallot)
library(plot3D)
library(viridis)

setwd("/home/chrisgg/julia/TimeDelays/")

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

sto_realization <- function(u0, param, sigma = 0.01, time = 1000, deltat = 1) {
  model.state <- u0
  model.sigma <- sigma
  model.time <- time
  model.deltat <- deltat
  ts.ex <- TSTraj(y0 = model.state,
                   time = model.time,
                   deltat = model.deltat,
                   x.rhs = param_model(param)[[1]],
                   y.rhs = param_model(param)[[2]],
                   sigma = model.sigma,
                   lower.bound = 0)
  return(ts.ex)
}

vector_decomp_plot <- function(qpotential, bounds.x, bounds.y, param) {
  VDAll <- VecDecomAll(surface = qpotential, x.rhs = param_model(param)[[1]], y.rhs = param_model(param)[[2]],x.bound = bounds.x, y.bound = bounds.y)
  ## Plot the deterministic skeleton vector field.
  ## VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
  par(mfrow=c(1,2))
  ## Plot the gradient vector field.
  VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
  ## Plot the remainder vector field.
  VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)
}

unit_conv <- 1/sqrt(2)

# Set bounds and step numbers for QPotential calculation
bounds.x <- c(0, 4/unit_conv)
bounds.y <- c(0, 4/unit_conv)
step.number.x <- 1000
step.number.y <- 1000



### Section for Gabriel to look over ###

#Quasi-potential with same noise for epsilon of 1 (with cycling) - using the R-M model parameterized above.
param_symnoise_ep1 <- list(eff = 0.85,
                           ep = 1.0,
                           res_g = 1.0,
                           con_f = 1.0)

# u0_symnoise_ep1 <- c(x = 1.1435249, y = 3.061860) # - starting values "on" the limit cycle
u0_symnoise_ep1 <- c(x = 0.686, y = 2.249) # starting values near the fixed point

ts_symnoise_ep1 <- sto_realization(u0 = u0_symnoise_ep1, param = param_symnoise_ep1)

TSPlot(ts_symnoise_ep1, deltat = 1)
TSPlot(ts_symnoise_ep1, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep1, dim = 1)
TSDensity(ts_symnoise_ep1, dim = 2)

qp_symnoise_ep1 <- QPotential(x.rhs = param_model(param_symnoise_ep1)[[1]],
                              x.start = u0_symnoise_ep1[1],
                              x.bound = bounds.x,
                              x.num.steps = step.number.x,
                              y.rhs = param_model(param_symnoise_ep1)[[2]],
                              y.start = u0_symnoise_ep1[2],
                              y.bound = bounds.y,
                              y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep1, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5)
# Gabe, this creates the 2d visualisation of the QPotential (where the consumer isocline would cross the resource isocline to the left of the hopf - so cycling). But as can be seen there is no bump in the middle around the (unstable) fixed point.

png(filename="figs/decomp_symnoise_ep1.png")
vector_decomp_plot(qp_symnoise_ep1, bounds.x, bounds.y, param_symnoise_ep1)
dev.off()

png(filename="figs/3dQPot_symnoise_ep1.png")
persp3D(z = qp_symnoise_ep1, x = seq(0,4,length.out = 1000), y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,3), col = viridis(n = 100, option = "A"), contour=TRUE,   xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 70)
dev.off()

## Now we compare the above QPotential with the Qpotential created by Abbott et al for a Consumer-Resource model (paramterized slightly differently). In the Nolting, Abbott article, this QPotential shows the bump around the unstable fixed point and the cycling.
## Note, I could not find the R code that created the qpotential in the article (just the mathematica code to make the figure). So this was my attempt.
## One issue is I'm not sure where they place the starting point for the ordered upwind method. Below I tried both on the cycle and next to the fixed point.
cycle.eqn.x <- "a * x * ( 1 - ( x / b ) ) - ( d * x * y ) / ( h + x )"
cycle.eqn.y <- "( e * x * y ) / ( h + x ) - m * y"
model.state <- c(x = 3, y = 3) # for the stochastic stochastic realization
model.parms <- c(a = 1.5, b = 45, d = 5, h = 18, m = 4, e = 10)
model.sigma <- 0.1
model.time <- 1000 # we used 2500 in the figures
model.deltat <- 1
ts.ex2 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat,x.rhs = cycle.eqn.x, y.rhs = cycle.eqn.y, parms = model.parms, sigma = model.sigma)
TSPlot(ts.ex2, deltat = model.deltat)                                  # Figure 8
TSPlot(ts.ex2, deltat = model.deltat, dim = 2, line.alpha = 25)        # Figure 9a
TSDensity(ts.ex2, dim = 1)                                             # Histogram
TSDensity(ts.ex2, dim = 2)

u0_abbottCR_cycle <- c(x = 14.10788, y = 15.24985) #starting values on the cycle
u0_abbottCR_fixedpoint <- c(x = 12, y = 6.6) #starting balues on the fixed point

eqn.x <- Model2String(cycle.eqn.x, parms = model.parms)
eqn.y <- Model2String(cycle.eqn.y, parms = model.parms)
eq1.qp <- QPotential(x.rhs = eqn.x,
                      x.start = u0_abbottCR_fixedpoint[1],
                      x.bound = c(-0.5, 30),
                      x.num.steps = 4000,
                      y.rhs = eqn.y,
                      y.start = u0_abbottCR_fixedpoint[2],
                      y.bound = c(-0.5, 20),
                      y.num.steps = 4000)

QPContour(eq1.qp, dens = c(1000, 1000), x.bound = c(-0.5, 30),y.bound = c(-0.5, 20), c.parm = 10)

png(filename="figs/3dqpot_cycle.png")
persp3D(z = eq1.qp,x = seq(-0.5,30,length.out = 4000), y = seq(-0.5,20,length.out = 4000), xlim=c(0,30), ylim=c(0,20), col = viridis(n = 100, option = "A"), contour=TRUE,   xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 70)
dev.off()

# From the 2d and 3d quasipotential functions, again looks like not getting the bump around the fixed point.

# Calculate all three vector fields.
VDAll <- VecDecomAll(surface = eq1.qp, x.rhs = eqn.x, y.rhs = eqn.y,x.bound = c(-0.5, 30), y.bound = c(-0.5, 20))
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = c(-0.5, 30), y.bound = c(-0.5, 20, xlim = c(0, 7.5), ylim = c(0, 7.5),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
png(filename="figs/decomp_cycle.png")
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = c(-0.5, 30), y.bound = c(-0.5, 20), arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
dev.off()

## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = c(-0.5, 30), y.bound = c(-0.5, 20), arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)

#Now if we compare to Example 3 in QPot: An R Package for Stochastic Differential Equation Quasi-Potential Analysis by Moore et al. ( a limit cycle ). We can see the bump around the fixed point.
cycle.eqn.x <- "- (y - beta) + mu * (x - alpha) * (1 - (x - alpha)^2 - (y - beta)^2)"
cycle.eqn.y <- "(x - alpha) + mu * (y - beta) * (1 - (x - alpha)^2 - (y - beta)^2)"
model.state <- c(x = 3, y = 3)
model.parms <- c(alpha = 4, beta = 5, mu = 0.2)
model.sigma <- 0.1
model.time <- 1000 # we used 2500 in the figures
model.deltat <- 0.005
ts.ex2 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat,x.rhs = cycle.eqn.x, y.rhs = cycle.eqn.y, parms = model.parms, sigma = model.sigma)
TSPlot(ts.ex2, deltat = model.deltat)                                  # Figure 8
TSPlot(ts.ex2, deltat = model.deltat, dim = 2, line.alpha = 25)        # Figure 9a
TSDensity(ts.ex2, dim = 1)                                             # Histogram
TSDensity(ts.ex2, dim = 2)

eqn.x <- Model2String(cycle.eqn.x, parms = model.parms)
eqn.y <- Model2String(cycle.eqn.y, parms = model.parms)
eq1.qp <- QPotential(x.rhs = eqn.x, x.start = 4.15611, x.bound = c(-0.5, 7.5),x.num.steps = 4000, y.rhs = eqn.y, y.start = 5.98774, y.bound = c(-0.5, 7.5),y.num.steps = 4000)

QPContour(eq1.qp, dens = c(1000, 1000), x.bound = c(-0.5, 7.5),y.bound = c(-0.5, 7.5), c.parm = 10)

png(filename="figs/3dqpot_limitcycle.png")
persp3D(z =eq1.qp, x = seq(-0.5,7.5,length.out = 4000), y = seq(-0.5,7.5,length.out = 4000),zlim = c(0,0.8), xlim=c(0,7), ylim=c(0,7), col = viridis(n = 100, option = "A"), contour=TRUE, xlab="X", ylab="Y", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 80)
dev.off()

# Calculate all three vector fields.
VDAll <- VecDecomAll(surface = eq1.qp, x.rhs = eqn.x, y.rhs = eqn.y,x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5))
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), xlim = c(0, 7.5), ylim = c(0, 7.5),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)




## Section for manuscript
#Quasi-potential with same noise for epsilon of 0.01 (with efficiency of 0.8 - canard)
param_symnoise_ep001_eff08 <- list(eff = 0.8,
                                   ep = 0.01,
                                   res_g = 1.0,
                                   con_f = 1.0)

u0_symnoise_ep001_eff08 <- c(x = 0.9131818, y = 2.28127)

ts_symnoise_ep001_eff08 <- sto_realization(u0 = u0_symnoise_ep001_eff08, param = param_symnoise_ep001_eff08)

TSPlot(ts_symnoise_ep001_eff08, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep001_eff08, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep001_eff08, dim = 1)
TSDensity(ts_symnoise_ep001_eff08, dim = 2)



qp_symnoise_ep001_eff08 <- QPotential(x.rhs = param_model(param_symnoise_ep001_eff08)[[1]],
                                      x.start = u0_symnoise_ep001_eff08[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep001_eff08)[[2]],
                                      y.start = u0_symnoise_ep001_eff08[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

#not working - need to change where qp began?
QPContour(surface = qp_symnoise_ep001_eff08, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)

#persp3D(z = qp_symnoise_ep001_eff08, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)

vector_decomp_plot(qp_symnoise_ep001_eff08, bounds.x, bounds.y, param_symnoise_ep001_eff08)

#TODO
#1. non canard limit cycle -> centre just a flat surface - shouldnt be that
#4. coordinate transform to get different noise for C and R

#Quasi-potential with same noise for epsilon of 0.1 (with efficiency of 0.5)
param_symnoise_ep01_eff05 <- list(eff = 0.5,
                                   ep = 0.1,
                                   res_g = unit_conv,
                                   con_f = unit_conv)

u0_symnoise_ep01_eff05 <- c(x = 2.02/unit_conv, y = 1.65/unit_conv)

ts_symnoise_ep01_eff05 <- sto_realization(u0 = u0_symnoise_ep01_eff05, param = param_symnoise_ep01_eff05)

TSPlot(ts_symnoise_ep01_eff05, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep01_eff05, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep01_eff05, dim = 1)
TSDensity(ts_symnoise_ep01_eff05, dim = 2)

qp_symnoise_ep01_eff05 <- QPotential(x.rhs = param_model(param_symnoise_ep01_eff05)[[1]],
                                      x.start = u0_symnoise_ep01_eff05[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep01_eff05)[[2]],
                                      y.start = u0_symnoise_ep01_eff05[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep01_eff05, dens = c(1000, 1000), x.bound = c(0,4),y.bound = c(0,4), c.parm = 5)

png(filename="figs/3dqpot_symnoise_ep01_eff05.png")
persp3D(z = qp_symnoise_ep01_eff05, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)
dev.off()

vector_decomp_plot(qp_symnoise_ep001_eff08, bounds.x, bounds.y, param_symnoise_ep01_eff05)

#Quasi-potential with same noise for epsilon of 0.5 (with efficiency of 0.5)
param_symnoise_ep05_eff05 <- list(eff = 0.5,
                                   ep = 0.5,
                                   res_g = unit_conv,
                                   con_f = unit_conv)

u0_symnoise_ep05_eff05 <- c(x = 2.02/unit_conv, y = 1.65/unit_conv)

ts_symnoise_ep05_eff05 <- sto_realization(u0 = u0_symnoise_ep05_eff05, param = param_symnoise_ep05_eff05)

TSPlot(ts_symnoise_ep05_eff05, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep05_eff05, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep05_eff05, dim = 1)
TSDensity(ts_symnoise_ep05_eff05, dim = 2)

qp_symnoise_ep05_eff05 <- QPotential(x.rhs = param_model(param_symnoise_ep05_eff05)[[1]],
                                      x.start = u0_symnoise_ep05_eff05[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep05_eff05)[[2]],
                                      y.start = u0_symnoise_ep05_eff05[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep05_eff05, dens = c(1000, 1000), x.bound = c(0,4),y.bound = c(0,4), c.parm = 5)

png(filename="figs/3dqpot_symnoise_ep05_eff05.png")
persp3D(z = qp_symnoise_ep05_eff05, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)
dev.off()

vector_decomp_plot(qp_symnoise_ep05_eff05, bounds.x, bounds.y, param_symnoise_ep05_eff05)

#Quasi-potential with same noise for epsilon of 0.9 (with efficiency of 0.5)
param_symnoise_ep09_eff05 <- list(eff = 0.5,
                                   ep = 0.9,
                                   res_g = unit_conv,
                                   con_f = unit_conv)

u0_symnoise_ep09_eff05 <- c(x = 2.02/unit_conv, y = 1.65/unit_conv)

ts_symnoise_ep09_eff05 <- sto_realization(u0 = u0_symnoise_ep09_eff05, param = param_symnoise_ep09_eff05)

TSPlot(ts_symnoise_ep09_eff05, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep09_eff05, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep09_eff05, dim = 1)
TSDensity(ts_symnoise_ep09_eff05, dim = 2)

qp_symnoise_ep09_eff05 <- QPotential(x.rhs = param_model(param_symnoise_ep09_eff05)[[1]],
                                      x.start = u0_symnoise_ep09_eff05[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep09_eff05)[[2]],
                                      y.start = u0_symnoise_ep09_eff05[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep09_eff05, dens = c(1000, 1000), x.bound = c(0,4),y.bound = c(0,4), c.parm = 5)

png(filename="figs/3dqpot_symnoise_ep09_eff05.png")
persp3D(z = qp_symnoise_ep09_eff05, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)
dev.off()

vector_decomp_plot(qp_symnoise_ep09_eff05, bounds.x, bounds.y, param_symnoise_ep09_eff05)

#Quasi-potential with same noise for epsilon of 0.04 (with efficiency of 0.6)
param_symnoise_ep004_eff06 <- list(eff = 0.6,
                                   ep = 0.04,
                                   res_g = unit_conv,
                                   con_f = unit_conv)

u0_symnoise_ep004_eff06 <- c(x = 1.30/unit_conv, y = 2.21/unit_conv)

ts_symnoise_ep004_eff06 <- sto_realization(u0 = u0_symnoise_ep004_eff06, param = param_symnoise_ep004_eff06)

TSPlot(ts_symnoise_ep004_eff06, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep004_eff06, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep004_eff06, dim = 2)

qp_symnoise_ep004_eff06 <- QPotential(x.rhs = param_model(param_symnoise_ep004_eff06)[[1]],
                                      x.start = u0_symnoise_ep004_eff06[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep004_eff06)[[2]],
                                      y.start = u0_symnoise_ep004_eff06[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep004_eff06, dens = c(1000, 1000), x.bound = c(0,4),y.bound = c(0,4), c.parm = 5)

png(filename="figs/3dqpot_symnoise_ep004_eff06.png")
persp3D(z = qp_symnoise_ep004_eff06, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)
dev.off()

vector_decomp_plot(qp_symnoise_ep004_eff06, bounds.x, bounds.y, param_symnoise_ep004_eff06)

#Quasi-potential with same noise for epsilon of 0.15 (with efficiency of 0.6)
param_symnoise_ep015_eff06 <- list(eff = 0.6,
                                   ep = 0.15,
                                   res_g = unit_conv,
                                   con_f = unit_conv)

u0_symnoise_ep015_eff06 <- c(x = 1.30/unit_conv, y = 2.21/unit_conv)

ts_symnoise_ep015_eff06 <- sto_realization(u0 = u0_symnoise_ep015_eff06, param = param_symnoise_ep015_eff06)

TSPlot(ts_symnoise_ep015_eff06, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep015_eff06, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep015_eff06, dim = 1)
TSDensity(ts_symnoise_ep015_eff06, dim = 2)

qp_symnoise_ep015_eff06 <- QPotential(x.rhs = param_model(param_symnoise_ep015_eff06)[[1]],
                                      x.start = u0_symnoise_ep015_eff06[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep015_eff06)[[2]],
                                      y.start = u0_symnoise_ep015_eff06[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep015_eff06, dens = c(1000, 1000), x.bound = c(0,4),y.bound = c(0,4), c.parm = 5)

png(filename="figs/3dqpot_symnoise_ep015_eff06.png")
persp3D(z = qp_symnoise_ep015_eff06, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)
dev.off()

vector_decomp_plot(qp_symnoise_ep015_eff06, bounds.x, bounds.y, param_symnoise_ep015_eff06)

#Quasi-potential with same noise for epsilon of 0.9 (with efficiency of 0.6)
param_symnoise_ep09_eff06 <- list(eff = 0.6,
                                   ep = 0.9,
                                   res_g = unit_conv,
                                   con_f = unit_conv)

u0_symnoise_ep09_eff06 <- c(x = 1.30/unit_conv, y = 2.21/unit_conv)

ts_symnoise_ep09_eff06 <- sto_realization(u0 = u0_symnoise_ep09_eff06, param = param_symnoise_ep09_eff06)

TSPlot(ts_symnoise_ep09_eff06, deltat = 1, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts_symnoise_ep09_eff06, deltat = 1, dim = 2)
TSDensity(ts_symnoise_ep09_eff06, dim = 1)
TSDensity(ts_symnoise_ep09_eff06, dim = 2)

qp_symnoise_ep09_eff06 <- QPotential(x.rhs = param_model(param_symnoise_ep09_eff06)[[1]],
                                      x.start = u0_symnoise_ep09_eff06[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_symnoise_ep09_eff06)[[2]],
                                      y.start = u0_symnoise_ep09_eff06[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_symnoise_ep09_eff06, dens = c(1000, 1000), x.bound = c(0,4),y.bound = c(0,4), c.parm = 5)

png(filename="figs/3dqpot_symnoise_ep09_eff06.png")
persp3D(z = qp_symnoise_ep09_eff06, x = seq(0,4,length.out = 1000),y = seq(0,4,length.out = 1000), xlim = c(0,3), ylim = c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab = "Resource", ylab = "Consumer", zlab = "Quasipotential", ticktype = "detailed", theta = 20, phi = 20)
dev.off()

vector_decomp_plot(qp_symnoise_ep09_eff06, bounds.x, bounds.y, param_symnoise_ep09_eff06)

#Quasi-potential with different noise (1:2) for epsilon of 0.04 (with efficiency of 0.6)
param_asymnoise_ep004_eff06 <- list(eff = 0.6,
                                   ep = 0.04,
                                   res_g = 0.447,
                                   con_f = 0.894)

u0_asymnoise_ep004_eff06 <- c(x = 1.30/0.447, y = 2.21/0.894)

ts_asymnoise_ep004_eff06 <- sto_realization(u0 = u0_asymnoise_ep004_eff06, param = param_asymnoise_ep004_eff06)

TSPlot(ts_asymnoise_ep004_eff06, deltat = 1, ylim = c(0,8), xlim = c(0,8))
TSPlot(ts_asymnoise_ep004_eff06, deltat = 1, dim = 2)
TSDensity(ts_asymnoise_ep004_eff06, dim = 2)

bounds.x <- c(0, 8)
bounds.y <- c(0, 8)

qp_asymnoise_ep004_eff06 <- QPotential(x.rhs = param_model(param_asymnoise_ep004_eff06)[[1]],
                                      x.start = u0_asymnoise_ep004_eff06[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_asymnoise_ep004_eff06)[[2]],
                                      y.start = u0_asymnoise_ep004_eff06[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_asymnoise_ep004_eff06, dens = c(1000, 1000), x.bound = c(0,8*0.447),y.bound = c(0,8*0.8894), c.parm = 5)

png(filename="figs/3dqpot_asymnoise_ep004_eff06.png")
persp3D(z =qp_asymnoise_ep004_eff06,x=seq(0,8*0.447,length.out = 1000),y=seq(0,8*0.894,length.out = 1000), xlim=c(0,3), ylim=c(0,3), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 40, phi = 20)
dev.off()

vector_decomp_plot(qp_asymnoise_ep004_eff06, bounds.x, bounds.y, param_symnoise_ep09_eff06)

#Quasi-potential with different noise (1:4) for epsilon of 0.04 (with efficiency of 0.6)
param_asymnoise14_ep09_eff06 <- list(eff = 0.6,
                                   ep = 0.04,
                                   res_g = 1/sqrt(17),
                                   con_f = 4/sqrt(17))

u0_asymnoise14_ep09_eff06 <- c(x = 1.30/(1/sqrt(17)), y = 2.21/(4/sqrt(17)))

ts_asymnoise14_ep09_eff06 <- sto_realization(u0 = u0_asymnoise_ep09_eff06, param = param_asymnoise_ep09_eff06)

TSPlot(ts_asymnoise14_ep09_eff06, deltat = 1, ylim = c(0,8), xlim = c(0,8))
TSPlot(ts_asymnoise14_ep09_eff06, deltat = 1, dim = 2)
TSDensity(ts_asymnoise14_ep09_eff06, dim = 2)

bounds.x <- c(0, 4/(1/sqrt(17)))
bounds.y <- c(0, 4/(4/sqrt(17)))

qp_asymnoise14_ep09_eff06 <- QPotential(x.rhs = param_model(param_asymnoise14_ep09_eff06)[[1]],
                                      x.start = u0_asymnoise14_ep09_eff06[1],
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_asymnoise14_ep09_eff06)[[2]],
                                      y.start = u0_asymnoise14_ep09_eff06[2],
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

QPContour(surface = qp_asymnoise14_ep09_eff06, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)

png(filename="figs/3dqpot_asymnoise14_ep004_eff06.png")
persp3D(z =qp_asymnoise14_ep09_eff06,x=seq(0,4,length.out = 1000),y=seq(0,4,length.out = 1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 40, phi = 20)
dev.off()

vector_decomp_plot(qp_asymnoise14_ep09_eff06, bounds.x, bounds.y, param_asymnoise14_ep09_eff06)


# Quasi-potential testing

var.eqn.x <- "(alpha * x) * (1 - (x / beta)) - ((delta * (x^2) * y) / (kappa + (x^2)))"
var.eqn.y <- "((gamma * (x^2) * y) / (kappa + (x^2))) - mu * (y^2)"
model.parms <- c(alpha = 1.54, beta = 10.14, delta = 1, gamma = 0.476,kappa = 1, mu = 0.112509)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms)
## Do not print to screen.
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)
model.state <- c(x = 1, y = 2)
model.sigma <- 0.05
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 0.025
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat,x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma)


TSPlot(ts.ex1, deltat = model.deltat)                     # Figure 2
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)            # Figure 3a
TSDensity(ts.ex1, dim = 1)                                # like Figure 2 histogram
TSDensity(ts.ex1, dim = 2)                                # Figure 3b

eq1.x <- 1.40491
eq1.y <- 2.80808
eq2.x <- 4.9040
eq2.y <- 4.06187
bounds.x <- c(-0.5, 20.0)
bounds.y <- c(-0.5, 20.0)
step.number.x <- 1000
step.number.y <- 1000

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = eq1.x, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = eq1.y,y.bound = bounds.y, y.num.steps = step.number.y)

eq2.local <- QPotential(x.rhs = parms.eqn.x, x.start = eq2.x, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = eq2.y,y.bound = bounds.y, y.num.steps = step.number.y)

ex1.global <- QPGlobal(local.surfaces = list(eq1.local, eq2.local),unstable.eq.x = c(0, 4.2008), unstable.eq.y = c(0, 4.0039),x.bound = bounds.x, y.bound = bounds.y)

QPContour(surface = ex1.global, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)

persp3D(z =ex1.global, xlim=c(0,0.5), ylim=c(0,0.5), col = viridis(n = 100, option = "A"), contour=TRUE)

# Calculate all three vector fields.
VDAll <- VecDecomAll(surface = ex1.global, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 11), ylim = c(0, 6),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)

## Deterministic Skeleton Analysis
library(deSolve)
library(phaseR)
library(rSymPy)


sympy("sympify(( r * x * (1 - ( x / k ) ) ) -  ( ( a * x * y ) / ( 1 +  ( a * h * x ) ) ))")
sympy("sympify(( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ))")
r <- Var("r")
x <- Var("x")
k <- Var("k")
a <- Var("a")
h <- Var("h")
m <- Var("m")
y <- Var("y")
e <- Var("e")


sympy("solve( ( r * x * (1 - ( x / k ) ) )  -  ( ( a * x * y ) / ( 1 +  ( a * h * x ) ) ), y )")
sympy("solve(( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) , x)")


# Ensure model is working correctly
RM <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- ( r * X * (1 - ( X / k ) ) ) -  ( ( a * X * Y ) / ( 1 +  ( a * h * X ) ) )
    dY <- ( ( e * a * X * Y ) / (1 + ( a * h * X ) ) ) - ( m * Y )

    list(c(dX, dY))
  })
}

RMep <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- ( r * X * (1 - ( X / k ) ) ) -  ( ( a * X * Y ) / ( 1 +  ( a * h * X ) ) )
    dY <- p * ( ( ( e * a * X * Y ) / (1 + ( a * h * X ) ) ) - ( m * Y ) )

    list(c(dX, dY))
  })
}

RM_abbott <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X * ( 1 - ( X / b ) ) - ( d * X * Y ) / ( h + X )
    dY <- ( e * X * Y ) / ( h + X ) - m * Y

    list(c(dX, dY))
  })
}

times <- seq(0, 600, by = 1)
initstate = c(X = 0.75, Y = 2.26)
parmslc <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.85, m = 0.4)

initstate_abbott = c(X = 3, Y = 3)
parms_abbott <- c(a = 1.5, b = 45, d = 5, h = 18, m = 4, e = 10)

parmslcep <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.8, m = 0.4, p = 0.01)
solved <- ode(y = initstate, times = times, func = RM, parms = parmslc)
plot(solved[,2],solved[,3])



# Questions:
#When does ordered upwind method stop - when does Quasipotential end?
