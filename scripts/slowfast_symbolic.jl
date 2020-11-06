#Script for slow fast examination of time delays
using SymPy

# Find equilibria - should be same as normal but also with ε = 0
x, y, r, k, a, m, e, h = symbols("x, y, r, k, a, m, e, h", real = true)

f(x, y) = r * x * (1 - x / k) - a * x * y / (1 + a * h * x)
g(x ,y) = e * a * x * y / (1 + a * h * x) - m * y

SymPy.solve(f(x,y),x)
SymPy.solve(f(x,y),y)
SymPy.solve(g(x,y),x)
SymPy.solve(g(x,y),y)

c(e) = ( 0.4 / (1.1 * (e - 0.8 * 0.4)) ) - 3

SymPy.solve(c(e),e)



#Find isoclines of Type I
@vars R C
@vars r k a m e h

f(R, C) = r * R * (1 - R / k) - a * R * C
g(R ,C) = e * a * R * C - m * R

SymPy.solve(f(R,C),R)
SymPy.solve(f(R,C),C)
SymPy.solve(g(R,C),R)
SymPy.solve(g(R,C),C)


f(R) = r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)

SymPy.diff(f(R), R)

SymPy.solve(SymPy.diff(f(R), R),R)

3.0 / 2 - 1 / ( 2 * 1.1 * 0.8)

2.0 * (1.1 * 0.8 * 3.0 * 0.93181818181819 - 1.1 * 0.8 * 0.93181818181819^2 + 3.0 - 0.93181818181819) / (1.1 * 3.0)
