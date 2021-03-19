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


##proving that when reduce the epsilon value we shift where the complex divide occurs
x, y, r, k, a, m, e, h, ϵ = symbols("x, y, r, k, a, m, e, h, ϵ", real = true)
f(x, y) = r * x * (1 - x / k) - a * x * y / (1 + a * h * x)
g(x ,y) = ϵ * (e * a * x * y / (1 + a * h * x) - m * y )


# first step find the jacobian

f1 = simplify(diff(f(x,y), x))
f2 = simplify(diff(f(x,y), y))
g1 = simplify(diff(g(x,y), x))
g2 = simplify(diff(g(x,y), y))

tr = SymPy.simplify(f1 + 0) #trying kevin's idea of leaving in R* and C*
det = SymPy.simplify((f1 * 0) - (f2 * g1))

trsq = SymPy.simplify(tr^2)
det4 = SymPy.simplify(4 * det)
imag = SymPy.simplify(trsq - det4)

eq_R = m / (a * (e - h * m))
eq_C = SymPy.simplify(r * (a * h * k * (m / (a * (e - h * m))) - a * h * (m / (a * (e - h * m)))^2 + k - (m / (a * (e - h * m)))) / (a * k))

f1e = SymPy.simplify( r - ( (2 * r * eq_R) / k )- ( (a * eq_C) / (1 + a * h * eq_R)^2))
f2e = SymPy.simplify( - a * eq_R / (1 + a * h * eq_R))
g1e = SymPy.simplify( ϵ * e * a * eq_C / (1 + a * h * eq_R)^2)
g2e = SymPy.simplify( ϵ * e * a * eq_R / (1 + a * h * eq_R) - ϵ * m)

tr = SymPy.simplify(f1e + g2e)
det = SymPy.simplify((f1e * g2e) - (f2e * g1e))

trsq = SymPy.simplify(tr^2)
det4 = SymPy.simplify(4 * det)
imag = SymPy.simplify(trsq - det4)

SymPy.solve(imag, ϵ)

SymPy.simplify(SymPy.diff((-m*r*(-a*e*h*k + a*h^2*k*m + e + h*m)^2)/(4*a*e*k*(e - h*m)^2*(-a*e*k + a*h*k*m + m)), e))

#try keeping R* in and do the same method as above to compare bsquared = 4ac
x, y, r, k, a, m, e, h, ϵ, R = symbols("x, y, r, k, a, m, e, h, ϵ, R", real = true)
f(x, y) = r * x * (1 - x / k) - a * x * y / (1 + a * h * x)
g(x ,y) = ϵ * (e * a * x * y / (1 + a * h * x) - m * y )


# first step find the jacobian

f1 = simplify(diff(f(x,y), x))
f2 = simplify(diff(f(x,y), y))
g1 = simplify(diff(g(x,y), x))
g2 = simplify(diff(g(x,y), y))

eq_R = m / (a * (e - h * m))
eq_C = SymPy.simplify(r * (a * h * k * R - a * h * R^2 + k - R) / (a * k))

f1e = SymPy.simplify( r - ( (2 * r * R) / k )- ( (a * eq_C) / (1 + a * h * R)^2))
f2e = SymPy.simplify( - a * R / (1 + a * h * R))
g1e = SymPy.simplify( ϵ * e * a * eq_C / (1 + a * h * R)^2)
g2e = SymPy.simplify( ϵ * e * a * eq_R / (1 + a * h * eq_R) - ϵ * m)

tr = SymPy.simplify(f1e + g2e)
det = SymPy.simplify((f1e * g2e) - (f2e * g1e))

trsq = SymPy.simplify(tr^2)
det4 = SymPy.simplify(4 * det)
imag = SymPy.simplify(trsq - det4)

using Parameters
@with_kw mutable struct RozMacPar
    r = 2.0
    k = 3.0
    a = 1.1
    h = 0.8
    e = 0.7
    m = 0.4
    σ = 0.1
    ε = 0.1
end

function test(p)
    @unpack r, a, k, h, m, e= p
    return h * k - ( (m / (a * (e - h * m))) * (e + h))
end

test(RozMacPar(e=0.45))^2
test(RozMacPar(e=0.5))^2
test(RozMacPar(e=0.52))^2
test(RozMacPar(e=0.75))^2

function test2(p)
    @unpack r, a, k, h, m, e= p
    R = m / (a * (e - h * m))
    return (R * r * a * (h * k - e * R - h * R)^2 ) / (4 * k *e * (k - R))
end


test2(RozMacPar(e=0.45))
test2(RozMacPar(e=0.5))
test2(RozMacPar(e=0.52))
test2(RozMacPar(e=0.75))
# need to simplify - maybe find eigenvalues of normal model first and put epsilon back in - or find part that is complex and deal only with this - or pick parameters (apart from e and ϵ)
# need to input the interior equlibrium

x, y, r, k, a, m, e, h = symbols("x, y, r, k, a, m, e, h", real = true)

i(x, y) = r * x * (1 - x / k) - a * x * y / (1 + a * h * x)
j(x ,y) = e * a * x * y / (1 + a * h * x) - m * y

i1 = diff(i(x,y), x)
i2 = diff(i(x,y), y)
j1 = diff(j(x,y), x)
j2 = diff(j(x,y), y)

simplify(solve((i1 - λ) * (j2 - λ) - (i2 * j1), λ))

#or determinant trace simplification?
simplify(i1 + j2)
simplify((i1 *j2)-(i2*j1))

simplify((i1 + j2)^2 - 4 * ((i1 *j2)-(i2*j1)) )


#attempt where i input parameter values apart from ϵ and efficiency
x, y, r, k, a, m, e, h, ϵ = symbols("x, y, r, k, a, m, e, h, ϵ", real = true)
f(x, y) = 2.0 * x * (1 - x / 3.0) - 1.1 * x * y / (1 + 1.1 * h * x)
g(x ,y) = ϵ * (e * 1.1 * x * y / (1 + 1.1 * h * x) - 0.4 * y )


f1 = diff(f(x,y), x)
f2 = diff(f(x,y), y)
g1 = diff(g(x,y), x)
g2 = diff(g(x,y), y)

λ = symbols("λ", real = true)

simplify(solve((f1 - λ) * (g2 - λ) - (f2 * g1), λ))

eq_R = m / (a * (e - h * m))
eq_C = SymPy.simplify(r * (a * h * k * (m / (a * (e - h * m))) - a * h * (m / (a * (e - h * m)))^2 + k - (m / (a * (e - h * m)))) / (a * k))


f1e = SymPy.simplify( r - ( (2 * r * eq_R) / k )- ( (a * eq_C) / (1 + a * h * eq_R)^2))
f2e = SymPy.simplify( - a * eq_R / (1 + a * h * eq_R))
g1e = SymPy.simplify( ϵ * e * a * eq_C / (1 + a * h * eq_R)^2)
g2e = SymPy.simplify( ϵ * e * a * eq_R / (1 + a * h * eq_R) - ϵ * m)

tr = SymPy.simplify(f1e^2)

det = SymPy.simplify(4 * (- (f2e * g1e)))

imag = SymPy.simplify(tr - det)
.
SymPy.solve(imag, e)



# Non-dimensionalized form
x, y, α, β, c, ϵ = symbols("x, y, α, β, c, ϵ", real = true)

f(x, y) = x * (1 - x ) - x * y / (1 + c * x)
g(x ,y) = ϵ * ( α * x * y / (1 + c * x) - β * y )

SymPy.solve(f(x,y),x)
SymPy.solve(f(x,y),y)
SymPy.solve(g(x,y),x)
SymPy.solve(g(x,y),y)

eq_R = -β/(c*β - α)
eq_C = SymPy.simplify(-c*(-β/(c*β - α))^2 + c*(-β/(c*β - α)) - (-β/(c*β - α)) + 1)
#eq_C_red = -α / (β * c - α) - β * α / ((β * c - α)^2)

f1 = simplify(diff(f(x,y), x))
f2 = simplify(diff(f(x,y), y))
g1 = simplify(diff(g(x,y), x))
g2 = simplify(diff(g(x,y), y))

f1e = SymPy.simplify( - 2 * eq_R - (eq_C / (c * eq_R + 1)^2)+1)
#f1e_2 = SymPy.simplify(- 2 * eq_R - (eq_C_red / (c * eq_R + 1)^2) + 1 )
f2e = SymPy.simplify( - eq_R / (c * eq_R +1))
g1e = SymPy.simplify( ϵ * eq_C * α / ((c * eq_R + 1 )^2))
g2e = SymPy.simplify( ϵ * (eq_R * α - β * ( c * eq_R + 1))/ ( c * eq_R + 1))

tr = SymPy.simplify(f1e + g2e)
det = SymPy.simplify((f1e * g2e) - (f2e * g1e))

trsq = SymPy.simplify(tr^2)
det4 = SymPy.simplify(4 * det)
imag = SymPy.simplify(trsq - det4)

SymPy.solve(imag, ϵ)

SymPy.simplify(SymPy.diff((-β*(c^2*β - c*α + c*β + α)^2)/(4*α*(c*β - α)^2*(c*β - α + β)),α))


SymPy.solve( 4 * α * ϵ * (c * β - α)^2 * (c * β - α + β) + β * (c^2 * β - c * α + c * β + α)^2, α)

#another idea is to use logistic plus type 1 to make proof simpler - still shows what we want to show

x, y, α, β, ϵ = symbols("x, y, α, β, ϵ", real = true)

f(x, y) = x * (1 - x ) - x * y
g(x ,y) = ϵ * ( α * x * y - β * y )

SymPy.solve(f(x,y),x)
SymPy.solve(f(x,y),y)
SymPy.solve(g(x,y),x)
SymPy.solve(g(x,y),y)

eq_R = β/α
eq_C = 1-β/α 

f1 = simplify(diff(f(x,y), x))
f2 = simplify(diff(f(x,y), y))
g1 = simplify(diff(g(x,y), x))
g2 = simplify(diff(g(x,y), y))

f1e = SymPy.simplify( - 2 * eq_R - eq_C +1)
f2e = SymPy.simplify( - eq_R )
g1e = SymPy.simplify( ϵ * eq_C * α )
g2e = SymPy.simplify( ϵ * (eq_R * α - β ))

tr = SymPy.simplify(f1e)
tr^2
det = SymPy.simplify((f1e * g2e) - (f2e * g1e))

imag = SymPy.simplify(tr - det)

SymPy.solve(imag, α)

using PyPlot
using Parameters

function real_complex(α,β,ϵ)
    return β * (- 4 * α * ϵ * (α - β) + β )
end

m = 0.4
r = 2.0
k = 3.0
a = 1.1
e = 0.7

(3.0 * 1.1 * 0.7) / 2.0

0.4/2
real_complex(1.155,0.2,0.01)

let
    β = 0.2
    alpharange = 0.0:0.01:0.9
    epsilonsmall = [real_complex(α, β, 0.1) for α in alpharange]
    epsilonmed = [real_complex(α, β, 0.2) for α in alpharange]
    epsilonlarge = [real_complex(α, β, 0.4) for α in alpharange]
    test = figure()
    plot(alpharange, epsilonsmall, color = "orange")
    plot(alpharange, epsilonmed, color = "red")
    plot(alpharange, epsilonlarge, color = "blue")
    hlines(0.0, 0.0, 0.9)
    xlabel("α")
    ylabel("Function in sqrt")
    return test
end



function real_complex_ndII(α,β,c,ϵ)
    return ( β * (  4 * α * ϵ * (c * β - α)^2 * (c * β - α + β) + β * (c^2 * β - c * α + c * β + α)^2 ) ) / (α^2 * (c * β - α)^2 )
end

m = 0.4
r = 2.0
k = 3.0
a = 1.1
e = 0.7
h = 0.8


0.4/2

1.1*0.8*3.0
real_complex_ndII(1.155,0.2,2.6400000006, 0.01)

let
    β = 0.2
    c = 2.6400000000000006
    alpharange = 0.0:0.01:3.0
    epsilonsmall = [real_complex_ndII(α, β, c, 1) for α in alpharange]
    #epsilonmed = [real_complex(α, β, 0.2) for α in alpharange]
    #epsilonlarge = [real_complex(α, β, 0.4) for α in alpharange]
    test = figure()
    plot(alpharange, epsilonsmall, color = "orange")
    #plot(alpharange, epsilonmed, color = "red")
    #plot(alpharange, epsilonlarge, color = "blue")
    hlines(0.0, 0.0, 3.0)
    xlabel("α")
    ylabel("Function in sqrt")
    return test
end



x, y, α, β, ϵ = symbols("x, y, α, β, ϵ", real = true)

SymPy.simplify(diff(sqrt(β* ϵ * (β* ϵ + 1)), ϵ))
diff(sqrt(β* ϵ * (β* ϵ + 1)), ϵ)

SymPy.simplify(diff(sqrt(β* ϵ * (β* ϵ + 1))/(2 * ϵ), ϵ))