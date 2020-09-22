library(QPot)
library(deSolve)
library(phaseR)
library(rSymPy)
library(plot3D)
library(viridis)

# Questions:
#When does ordered upwind method stop - when does Quasipotential end?


#Quasi-potential with same noise for epsilon of 1 (with cycling)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.72, e = 0.8, m = 0.4)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 1.1922515, y= 2.367524)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 0.01
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma)

TSPlot(ts.ex1, deltat = model.deltat)
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)   
bounds.x <- c(-0.5, 5.0)
bounds.y <- c(-0.5, 5.0)
step.number.x <- 1000
step.number.y <- 1000

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 1.1922515, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 2.367524, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
# From gradient vector field does not look like upwind method correctly tested area within limit cycle - why?
persp3D(z =eq1.local, x=seq(0,4,length.out = 1000), y=seq(0,4,length.out = 1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,   xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 20)

par(mfrow=c(1, 1))
# Calculate all three vector fields.
VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)



#Quasi-potential with same noise for epsilon of 0.01 (with efficiency of 0.8 - canard)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.8, m = 0.4, p = 0.1)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 0.9131818, y = 2.28127)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 0.9131818, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 2.28127, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
#not quite working
persp3D(z =eq1.local, x=seq(0,4,length.out = 1000),y=seq(0,4,length.out = 1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)



# maybe doesn't matter that same noise for both equations - could put different noise in supplemental
# NO BUT noise on different time scales so could potentially matter!

#TODO try making for different epsilon values across the real/complex divide - AND set out issues
#1. non canard limit cycle -> centre just a flat surface - shouldnt be that
#2. check determinitistic canard skeleton - always same numbers?
#3. canard with stochasticity doesn't follow same straight back to resource isocline (going right) every single time 
#4. coordinate transform to get different noise for C and R

#Quasi-potential with same noise for epsilon of 0.1 (with efficiency of 0.5)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.5, m = 0.4, p = 0.01)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 2.02, y = 1.65)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)
step.number.x <- 1000
step.number.y <- 1000

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 2.02, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 1.65, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
#not quite working
persp3D(z =eq1.local, x=seq(0, 4, length.out=1000), y=seq(0, 4, length.out=1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)

#Quasi-potential with same noise for epsilon of 0.5 (with efficiency of 0.5)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.5, m = 0.4, p = 0.5)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 2.02, y = 1.65)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 2.02, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 1.65, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
#not quite working
persp3D(z =eq1.local,x=seq(0, 4, length.out=1000), y=seq(0, 4, length.out=1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)

#Quasi-potential with same noise for epsilon of 0.9 (with efficiency of 0.5)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.5, m = 0.4, p = 0.9)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 2.02, y = 1.65)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)
step.number.x <- 1000
step.number.y <- 1000

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 2.02, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 1.65, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
#not quite working
persp3D(z =eq1.local, x=seq(0, 4, length.out=1000), y=seq(0, 4, length.out=1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 20, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)



#Quasi-potential with same noise for epsilon of 0.04 (with efficiency of 0.6)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.6, m = 0.4, p = 0.04)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 1.30, y = 2.21)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)
step.number.x <- 1000
step.number.y <- 1000

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 1.30, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 2.21, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)

persp3D(z =eq1.local, x=seq(0,4,length.out = 1000), y=seq(0,4,length.out = 1000), xlim=c(0,3), ylim=c(0,2.5), zlim=c(0,0.0007), col = viridis(n = 100, option = "A"), contour=TRUE, xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 40, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4), arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)


#Quasi-potential with same noise for epsilon of 0.15 (with efficiency of 0.6)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.6, m = 0.4, p = 0.15)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 1.30, y = 2.21)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 1.30, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 2.21, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
#not quite working
persp3D(z =eq1.local,x=seq(0,4,length.out = 1000),y=seq(0,4,length.out = 1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 40, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)


#Quasi-potential with same noise for epsilon of 0.9 (with efficiency of 0.6)
var.eqn.x <- "( r * x * (1 - ( x / k ) ) ) - ( a * x * y ) / ( 1 +  ( a * h * x ) ) "
var.eqn.y <- "p * ( ( e * a * x * y ) / (1 + ( a * h * x ) ) - ( m * y ) )"

model.parms <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.6, m = 0.4, p = 0.9)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms, supress.print = TRUE)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, supress.print = TRUE)

model.state <- c(x = 1.30, y = 2.21)
model.sigma <- 0.01
model.time <- 1000     # we used 12500 in the figures
model.deltat <- 1
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma, lower.bound = 0)
#Problem - if using noise on resource - dynamics can flick over to negative because so close to 0 resources for part of canard
# COuld have been solved with lower.bound = 0

TSPlot(ts.ex1, deltat = model.deltat, ylim = c(0,5), xlim = c(0,5))
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
bounds.x <- c(0, 4)
bounds.y <- c(0, 4)
step.number.x <- 1000
step.number.y <- 1000

eq1.local <- QPotential(x.rhs = parms.eqn.x, x.start = 1.30, x.bound = bounds.x,x.num.steps = step.number.x, y.rhs = parms.eqn.y, y.start = 2.21, y.bound = bounds.y, y.num.steps = step.number.y)

QPContour(surface = eq1.local, dens = c(1000, 1000), x.bound = bounds.x,y.bound = bounds.y, c.parm = 5)
#not quite working
persp3D(z =eq1.local,x=seq(0,4,length.out = 1000),y=seq(0,4,length.out = 1000), xlim=c(0,3), ylim=c(0,2.5), col = viridis(n = 100, option = "A"), contour=TRUE,  xlab="Resource", ylab="Consumer", zlab="Quasipotential", ticktype="detailed", theta = 40, phi = 20)

VDAll <- VecDecomAll(surface = eq1.local, x.rhs = parms.eqn.x, y.rhs = parms.eqn.y,x.bound = bounds.x, y.bound = bounds.y)
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 4), ylim = c(0, 4),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)




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

#Limit cycle
var.eqn.x <- "- (y - beta) + mu * (x - alpha) * (1 - (x - alpha)^2 - (y - beta)^2)"
var.eqn.y <- "(x - alpha) + mu * (y - beta) * (1 - (x - alpha)^2 - (y - beta)^2)"
model.state <- c(x = 3, y = 3)
model.parms <- c(alpha = 4, beta = 5, mu = 0.2)
model.sigma <- 0.1
model.time <- 1000 # we used 2500 in the figures
model.deltat <- 0.005
ts.ex2 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat,x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)
TSPlot(ts.ex2, deltat = model.deltat)                                  # Figure 8
TSPlot(ts.ex2, deltat = model.deltat, dim = 2, line.alpha = 25)        # Figure 9a
TSDensity(ts.ex2, dim = 1)                                             # Histogram
TSDensity(ts.ex2, dim = 2) 

eqn.x <- Model2String(var.eqn.x, parms = model.parms)
eqn.y <- Model2String(var.eqn.y, parms = model.parms)
eq1.qp <- QPotential(x.rhs = eqn.x, x.start = 4.15611, x.bound = c(-0.5, 7.5),x.num.steps = 4000, y.rhs = eqn.y, y.start = 5.98774, y.bound = c(-0.5, 7.5),y.num.steps = 4000)

QPContour(eq1.qp, dens = c(1000, 1000), x.bound = c(-0.5, 7.5),y.bound = c(-0.5, 7.5), c.parm = 10)

#persp3D(z =eq1.qp, xlim=c(0,7), ylim=c(0,7), col = viridis(n = 100, option = "A"), contour=TRUE)

# Calculate all three vector fields.
VDAll <- VecDecomAll(surface = eq1.qp, x.rhs = eqn.x, y.rhs = eqn.y,x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5))
## Plot the deterministic skeleton vector field.
VecDecomPlot(x.field = VDAll[, , 1], y.field = VDAll[, , 2], dens = c(25, 25),x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), xlim = c(0, 7.5), ylim = c(0, 7.5),arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
## Plot the gradient vector field.
VecDecomPlot(x.field = VDAll[, , 3], y.field = VDAll[, , 4], dens = c(25, 25),x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), arrow.type = "proportional",tail.length = 0.25, head.length = 0.025)
## Plot the remainder vector field.
VecDecomPlot(x.field = VDAll[, , 5], y.field = VDAll[, , 6], dens = c(25, 25),x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), arrow.type = "proportional",tail.length = 0.35, head.length = 0.025)


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

times <- seq(0, 600, by = 1)
initstate = c(X = 0.75, Y = 2.26)
parmslc <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.72, m = 0.4)

parmslcep <- c(r = 2.0, k = 3.0, a = 1.1, h = 0.8, e = 0.8, m = 0.4, p = 0.01)
solved <- ode(y = initstate, times = times, func = RM, parms = parmslc)
plot(solved)
