#Script for slow fast examination of time delays

include("packages.jl")
using SymPy

#Non-dimensionalized
@with_kw mutable struct RozMacPar
    r = 2.0
    K = 3.0
    a = 1.1
    h = 0.8
    e = 0.7
    m = 0.4
    σ = 0.1
    ε = 0.1
end

par_rozmac = RozMacPar()
function roz_mac_II!(du, u, p, t,)
    @unpack r, K, a, h, e, m = p
    R, C = u
    du[1] = r * R * (1 - R / K) - a * R * C / (1 + a * h * R)
    du[2] = ε ( e * a * R * C / (1 + a * h * R) - m * C )
    return
end

function roz_mac_II(u, par)
    du = similar(u)
    roz_mac_II!(du, u, par, 0.0)
    return du
end

# Find equilibria - should be same as normal but also with ε = 0
x, y, r, k, a, m, e, h = symbols("x, y, r, k, a, m, e, h", real = true)

f(x, y) = r * x * (1 - x / k) - a * x * y / (1 + a * h * x)
g(x ,y) = e * a * x * y / (1 + a * h * x) - m * y

SymPy.solve(f(x,y),x)
SymPy.solve(f(x,y),y)
SymPy.solve(g(x,y),x)
SymPy.solve(g(x,y),y)
# Find isoclines
function con_iso(p)
    @unpack m, a, e, h = p
    m / (a * (e - h * m))
end

function res_iso(x, p)
    @unpack a, k, r, h = p
    r * (a * h * k * x - a * h * x^2 + k - x) / (a * k)
end

vline(0.4 / (1.1 * (0.7 - 0.8 * 0.4)))

plot(x^2,0,4)
let
    plot(con_iso(par_rozmac), 0, 4, label = "Consumer Isocline")
end

# Create vector fields
# - Before hopf fixed point

# - Before hopf damped oscillations

# - At hopf

# - After hopf (limit cycle)


## Dimensionalized
#Setup
@with_kw mutable struct RozMacPar
    r = 2.0
    K = 3.0
    a = 1.1
    h = 0.8
    e = 0.7
    m = 0.4
    σ = 0.1
end

function roz_mac_II!(du, u, p, t,)
    @unpack r, K, a, h, e, m = p
    R, C = u
    du[1] = r * R * (1 - R / K) - a * R * C / (1 + a * h * R)
    du[2] = e * a * R * C / (1 + a * h * R) - m * C
    return
end

# Find equilibria - should be same as normal but also with ε = 0

# Find isoclines

# Create vector fields
# - Before hopf fixed point

# - Before hopf damped oscillations

# - At hopf

# - After hopf (limit cycle)
