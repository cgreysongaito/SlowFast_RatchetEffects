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
    @unpack r, k, a, h, e, m, ε = p
    R, C = u
    du[1] = r * R * (1 - R / k) - a * R * C / (1 + a * h * R)
    du[2] = ε * ( e * a * R * C / (1 + a * h * R) - m * C )
    return
end

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
# Plot isoclines and vector fields

function roz_mac_res(R, C, p)
    @unpack r, k, h, a, m = p
    return r * R * (1 - R / k) - (a * R * C / (1 + a * h * R) )
end

function roz_mac_con(R, C, eff, ep, p)
    @unpack h, a, m = p
    return ep * ( ( eff * a * R * C ) / (1 + a * h * R) - m * C )
end

function con_iso(eff, p)
    @unpack m, a, h = p
    m / (a * (eff - h * m))
end

function res_iso(R, p)
    @unpack a, k, r, h = p
    r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
end

function roz_mac_plot(eff, ep)
    minval = 0
    maxval = 3
    steps = 100
    resconrange = range(minval,stop=maxval,length=steps)

    U = [roz_mac_res(R, C, par_rozmac) for C in resconrange, R in resconrange]
    V = [roz_mac_con(R, C, eff, ep, par_rozmac) for C in resconrange, R in resconrange]
    speed = sqrt.(U.^2 .+ V.^2)
    lw = 5 .* speed ./ maximum(speed) # Line Widths
    streamplot(collect(resconrange), collect(resconrange), U, V, density = 0.6, color = "k", linewidth = lw)
    PyPlot.plot(collect(resconrange), [res_iso(R, par_rozmac) for R in resconrange])
    return PyPlot.plot(repeat([con_iso(eff, par_rozmac)],100),collect(resconrange))
end

# - Before hopf fixed point
let
    figure(figsize = (8,10))
    subplot(421)
    roz_mac_plot(0.45, 1)
    subplot(422)
    roz_mac_plot(0.45, 0.1)
    subplot(423)
    roz_mac_plot(0.6, 1)
    subplot(424)
    roz_mac_plot(0.6, 0.1)
    subplot(425)
    roz_mac_plot(0.71, 1) #need to differentiate and find max
    subplot(426)
    roz_mac_plot(0.71, 0.1)
    subplot(427)
    roz_mac_plot(0.9, 1)
    subplot(428)
    roz_mac_plot(0.9, 0.1)
    savefig("figs/vectorfieldplot.png")
end

# Plot transients and measure length of transients
function roz_mac_ep_plot(eff,ep)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    u0 = [2.5, 1.5]
    tspan = (0.0, 500.0)

    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    return PyPlot.plot(sol.t, sol.u)
end

let
    figure()
    subplot(421)
    roz_mac_ep_plot(0.45, 1)
    subplot(422)
    roz_mac_ep_plot(0.45, 0.1)
    subplot(423)
    roz_mac_ep_plot(0.6, 1)
    subplot(424)
    roz_mac_ep_plot(0.6, 0.1)
    subplot(425)
    roz_mac_ep_plot(0.73, 1)
    subplot(426)
    roz_mac_ep_plot(0.73, 0.1)
    subplot(427)
    roz_mac_ep_plot(0.9, 1)
    subplot(428)
    roz_mac_ep_plot(0.9, 0.1)
    savefig("figs/transientsplot.png")
end

#0.1 epsilon decreases length of transients after hopf but no perceptible difference before hopf. after hopf less variability in the cycles

#measure length of transients (before hopf) - know equilibria values but then need to check stay at those values for more than one timestep. also starting conditions should be random

function eq_II(eff, p)
    @unpack r, a, k, h, m = p
    eq_II_R = m / (a * (eff - h * m))
    eq_II_C = r * (a * h * k * (m / (a * (eff - h * m))) - a * h * (m / (a * (eff - h * m)))^2 + k - m / (a * (eff - h * m))) / (a * k)
    return vcat(eq_II_R, eq_II_C)
end

function transient_length(effi, ep)
    par = RozMacPar()
    par.e = effi
    par.ε = ep
    eq = eq_II(par.e, par)
    u0 = [2.5, 1.5]
    tspan = (0.0, 100000.0)

    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)

    for i in 1:length(sol.t)
        if isapprox(sol.u[i], eq; atol = 1e-3) & isapprox(sol.u[i+1], eq; atol = 1e-3)
        return sol.t[i]
        end
    end
end

effrange = 0.45:0.01:0.71
eprange = 0.05:0.01:1

con = [transient_length(ei, epi) for ei in effrange, epi in eprange]

let
    figure()
    contourf(eprange, effrange, con,levels = collect(range(0, stop=36000, length = 100)))
    gcf()
    savefig("figs/transientlengthplot.png")
end


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
