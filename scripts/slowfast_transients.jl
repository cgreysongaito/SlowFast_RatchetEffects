#Script for slow fast examination of time delays

include("packages.jl")
include("slowfast_commoncode.jl")

# Plot transients and measure length of transients
# to create starting conditions eq * 1 + rand(Uniform(1e-7, 1e-6))


function roz_mac_ep_plot(eff,ep)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    eq = eq_II(par)
    u0 = randeq.(eq)
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
    roz_mac_ep_plot(0.75, 1)
    subplot(426)
    roz_mac_ep_plot(0.75, 0.1)
    subplot(427)
    roz_mac_ep_plot(0.9, 1)
    subplot(428)
    roz_mac_ep_plot(0.9, 0.1)
    savefig("figs/transientsplot.png")
end

#0.1 epsilon decreases length of transients after hopf but no perceptible difference before hopf. after hopf less variability in the cycles

#measure length of transients (before hopf) - know equilibria values but then need to check stay at those values for more than one timestep. also starting conditions should be random



function transient_length(effi, ep)
    par = RozMacPar()
    par.e = effi
    par.ε = ep
    eq = eq_II(par)
    u0 = randeq.(eq)
    tspan = (0.0, 100000.0)

    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)

    for i in 1:length(sol.t)
        if isapprox(sol.u[i], eq; atol = 1e-7) & isapprox(sol.u[i+1], eq; atol = 1e-7)
        return log(10, sol.t[i])
        end
    end
end

effrange = 0.45:0.01:0.71
eprange = 0.05:0.01:1

con = [transient_length(ei, epi) for ei in effrange, epi in eprange]

let
    figure()
    contourf(eprange, effrange, con)
    xlabel("ε")
    ylabel("Efficiency")
    gcf()
    #savefig("figs/transientlengthplot.png")
end

#need to add code to extend to after the hopf
