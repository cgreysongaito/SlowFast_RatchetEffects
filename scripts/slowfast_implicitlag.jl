#Script for slow fast examination of time delays

include("packages.jl")
include("lag_utils.jl")
include("slowfast_commoncode.jl")

# Compare changing epsilon and implicit lag
function epsilon_lag_plot(eff, lagboth)
    par = RozMacPar()
    par.e = eff
    eq = eq_II(par)
    epvals = 0.006:0.0005:1
    u0 = randeq.(eq)
    tspan = (0.0, 100000.0)
    tstart = 9000
    tend = 100000
    tstep = 0.1
    tvals = tstart:tstep:tend
    lag = fill(0.0, length(epvals))
    per = fill(0.0, length(epvals))
    pha = fill(0.0, length(epvals))

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
        asol = sol(tvals)
        R_peaks = find_peaks(asol[1, :], 1e-8)
        C_peaks = find_peaks(asol[2, :], 1e-8)
        period = diff(R_peaks[1])[1] * tstep
        phase = find_phase(tvals, asol[1, :], asol[2, :])
        delay = phase / period #is the delay a fraction of the full difference between peaks of resource?
        lag[epi] = delay
        per[epi] = period
        pha[epi] = phase
    end

    if lagboth == "lag"
        return PyPlot.plot(collect(epvals), lag)
    else
        PyPlot.plot(collect(epvals), pha, label = "phase")
        return PyPlot.plot(collect(epvals), per, label = "period")
    end
    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end



#PLACEHOLDER not sure it makes sense to graph lag before the hopf - could do with if else function where before the hopf determine time span before reach equlibrium use this to calculate lag - check with GABE but with random tiny perturbation no point in calculating lag etc because the transients is so short (before the hopf)

let
    figure()
    subplot(221)
    epsilon_lag_plot(0.74, "lag")
    # subplot(222)
    # epsilon_lag_plot(0.74, "both")
    # subplot(223)
    # epsilon_lag_plot(0.9, "lag")
    # xlabel("ε")
    # subplot(224)
    # epsilon_lag_plot(0.9, "both")
    # xlabel("ε")
    # gcf()
    savefig("figs/epsilon_lag_plot.png")
end



par = RozMacPar()
par.e = 0.74
eq = eq_II(par)
u0 = randeq.(eq)
tspan = (0.0, 100000.0)
tstart = 9000
tend = 100000
tstep = 0.1
tvals = tstart:tstep:tend


par.ε = 0.009
prob = ODEProblem(roz_mac_II!, u0, tspan, par)
sol = DifferentialEquations.solve(prob, reltol = 1e-8)
asol = sol(tvals)
R_peaks = find_peaks(asol[1, :], 1e-8)
C_peaks = find_peaks(asol[2, :], 1e-8)
period = diff(R_peaks[1])[1] * tstep
phase = find_phase(tvals, asol[1, :], asol[2, :])
delay = phase / period #is the delay a fraction of the full difference between peaks of resource?
