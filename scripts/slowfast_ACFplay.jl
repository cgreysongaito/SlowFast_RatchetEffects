let
    using Parameters
    using LinearAlgebra
    using Calculus
    using DifferentialEquations
    using ForwardDiff
    using Distributions
    using StatsBase
    using Random
    using PyPlot
end

include("slowfast_commoncode.jl")

#first one perturbation at start

par = RozMacPar()
par.ε = 1
par.e = 0.5
u0 = [eq_II(par)[1], eq_II(par)[2] * 1.01]
tspan = (0.0, 100.0)
tstart = 0.0
tend = 100.0
tstep = 1
tvals = tstart:tstep:tend

    # cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
prob = ODEProblem(roz_mac_II!, u0, tspan, par)

sol = DifferentialEquations.solve(prob, reltol = 1e-8)

soldis = sol(tvals)

lrange = 0:1:20

acf = autocor(soldis[2, 1:end], collect(lrange))

# conf = 1.96/sqrt(length(solend))
let
    figure()
    subplot(1,2,1)
    bar(collect(lrange), acf)
    # hlines(0 + conf, 0, maximum(lrange))
    # hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    ylabel("ACF")
    subplot(1,2,2)
    plot(soldis.t, soldis.u)
    gcf()
end

# second constant number

acf2 = autocor(fill(5.0,101), collect(lrange))

# conf = 1.96/sqrt(length(solend))
let
    figure()
    subplot(1,2,1)
    bar(collect(lrange), acf2)
    # hlines(0 + conf, 0, maximum(lrange))
    # hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    ylabel("ACF")
    subplot(1,2,2)
    plot(soldis.t, fill(5.0,101))
    gcf()
end

# three three perturbations spaced out

function pert_cb2(integrator)
    integrator.u[2] = integrator.u[2] * ( 1 + rand(Normal(0.0, 0.01)))
end

par = RozMacPar()
par.ε = 1.0
par.e = 0.5
u0 = [eq_II(par)[1], eq_II(par)[2] * ( 1 + rand(Normal(0.0, 0.01)))]
tspan = (0.0, 100.0)
tstart = 0.0
tend = 100.0
tstep = 1
tvals = tstart:tstep:tend

cb = PeriodicCallback(pert_cb2, 1, initial_affect = false)
prob = ODEProblem(roz_mac_II!, u0, tspan, par)


sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)


solpert = sol(tvals)

lrange = 0:1:50

acfpert = autocor(solpert[2, 1:end], collect(lrange))

conf = 1.96/sqrt(length(solpert))
let
    figure()
    subplot(1,2,1)
    bar(collect(lrange), acfpert)
    hlines(0 + conf, 0, maximum(lrange))
    hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    ylabel("ACF")
    subplot(1,2,2)
    plot(solpert.t, solpert.u)
    gcf()
end


### White noise with mean 0.0 and standard deviation 0.01

x = rand(Normal(0.0, 0.01), 101)


acf_wn = autocor(x, collect(lrange))

conf = 1.96/sqrt(101)
let
    figure()
    subplot(1,2,1)
    bar(collect(lrange), acf_wn)
    hlines(0 + conf, 0, maximum(lrange))
    hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    ylabel("ACF")
    subplot(1,2,2)
    plot(collect(0:1:100), x)
    gcf()
end
