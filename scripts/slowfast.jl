#Script for slow fast examination of time delays

include("packages.jl")

#Non-dimensionalized
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

let
    Plot.plot(x -> res_iso(x ,par_rozmac), 0, 4, label = "Resource Isocline")
    vline!([con_iso(par_rozmac)], label = "Consumer Isocline")
end

# Create vector fields
using PyPlot

minval = 0
maxval = 3
steps = 100
R = repeat(range(minval,stop=maxval,length=steps)',steps)
C = repeat(range(minval,stop=maxval,length=steps),1,steps)
U = 2.0 .* R .* (1 .- R ./ 3.0) .- 1.1 .* R .* C ./ (1 .+ 1.1 * 0.8 .* R)
V = 1 .* ( 0.7 .* 1.1 .* R .* C ./ (1 .+ 1.1 .* 0.8 .* R) .- 0.7 .* C )
speed = sqrt.(U.^2 .+ V.^2)
lw = 5 .* speed ./ maximum(speed) # Line Widths

let
    figure()
    streamplot(R,C,U,V,density=0.6,color="k",linewidth=lw)
    gcf()
end


function roz_mac_res(du, u, p, t,)
    @unpack r, K, a, h, e, m = p
    R, C = u
    du[1] = r * R * (1 - R / K) - a * R * C / (1 + a * h * R)
    du[2] = ε ( e * a * R * C / (1 + a * h * R) - m * C )
    return
end


minval = 0
maxval = 3
steps = 100
R = repeat(range(minval,stop=maxval,length=steps)',steps)
C = repeat(range(minval,stop=maxval,length=steps),1,steps)
U = 2.0 .* R .* (1 .- R ./ 3.0) .- 1.1 .* R .* C ./ (1 .+ 1.1 * 0.8 .* R)
V = 1 .* ( 0.7 .* 1.1 .* R .* C ./ (1 .+ 1.1 .* 0.8 .* R) .- 0.7 .* C )
speed = sqrt.(U.^2 .+ V.^2)
lw = 5 .* speed ./ maximum(speed) # Line Widths

let
    figure()
    streamplot(R,C,U,V,density=0.6,color="k",linewidth=lw)
    gcf()
end
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
