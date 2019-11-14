#Script for slow fast examination of time delays

include("packages.jl")
include("lag_utils.jl")

# Compare changing epsilon and implicit lag
function epsilon_lag_plot(eff, lagboth)
    par = RozMacPar()
    par.e = eff
    eq = eq_II(par.e, par)
    epvals = 0.05:0.01:1
    u0 = randeq.(eq)
    tspan = (0.0, 10000.0)
    tstart = 9000
    tend = 10000
    tstep = 0.1
    tvals = tstart:tstep:tend
    lag = fill(0.0, length(epvals))
    per = fill(0.0, length(epvals))
    pha = fill(0.0, length(epvals))
    htend = 0

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
        if eff < 0.73
            for i in 1:length(sol.t)-1
                if isapprox(sol.u[i], eq; atol = 1e-7) && isapprox(sol.u[i+1], eq; atol = 1e-7)
                    htend = i
                    break
                end
            end
            asol = sol(0:0.1:htend)
        else
            asol = sol(tvals)
        end
        println(htend)
        break asol
        R_peaks = find_peaks(asol[1, :], 1e-20)
        C_peaks = find_peaks(asol[2, :], 1e-20)
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

let
    figure()
    epsilon_lag_plot(0.6, "lag")
    gcf()
end #PLACEHOLDER not sure it makes sense to graph lag before the hopf - could do with if else function where before the hopf determine time span before reach equlibrium use this to calculate lag - check with GABE but with random tiny perturbation no point in calculating lag etc because the transients is so short (before the hopf)

epsilon_lag_plot(0.6, "lag")

let
    figure()

    gcf()
end

let
    figure()
    subplot(221)
    epsilon_lag_plot(0.74, "lag")
    subplot(222)
    epsilon_lag_plot(0.74, "both")
    subplot(223)
    epsilon_lag_plot(0.9, "lag")
    xlabel("ε")
    subplot(224)
    epsilon_lag_plot(0.9, "both")
    xlabel("ε")
    gcf()
    savefig("figs/epsilon_lag_plot.png")
end
