#Script for slow fast examination of time delays
cd("/home/chrisgg/Documents/Guelph/PhD/TimeDelays/")

include("scripts/packages.jl")
include("scripts/slowfast_commoncode.jl")

### NOTE this needs to be fixed - incorrect calc of CV and old way for producing prob


## Compare epsilon and coefficient of variation
# PLACEHOLDER should CV be calculated before the hopf
function epsilon_cv_plot(eff, rescon, cvboth)
    par = RozMacPar()
    par.e = eff
    eq = eq_II(par)
    epvals = 0.05:0.001:1
    u0 = randeq.(eq)
    tspan = (0.0, 100000.0)
    tstart = 90000
    tend = 100000
    tstep = 0.1
    tvals = tstart:tstep:tend
    cv = fill(0.0, length(epvals))
    mn = fill(0.0, length(epvals))
    sd = fill(0.0, length(epvals))

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
        asol = sol(tvals)
        if rescon == "res"
            cv[epi] = std(asol[1, 1:end] ./ mean(asol[1, 1:end]))
            mn[epi] = mean(asol[1, 1:end])
            sd[epi] = std(asol[1, 1:end])
        else
            cv[epi] = std(asol[1, 1:end] ./ mean(asol[2, 1:end]))
            mn[epi] = mean(asol[2, 1:end])
            sd[epi] = std(asol[2, 1:end])
        end
    end

    if cvboth == "cv"
    return PyPlot.plot(collect(epvals), cv)
    else
       PyPlot.plot(collect(epvals), mn, label = "mean")
       return PyPlot.plot(collect(epvals), sd, label = "sd")
    end
end

let
    figure()
    subplot(221)
    PyPlot.title("Resource")
    epsilon_cv_plot(0.74, "res", "cv")
    ylabel("CV")
    subplot(222)
    PyPlot.title("Consumer")
    epsilon_cv_plot(0.74, "con", "cv")
    subplot(223)
    epsilon_cv_plot(0.9, "res", "cv")
    xlabel("ε")
    ylabel("CV")
    subplot(224)
    epsilon_cv_plot(0.9, "con", "cv")
    xlabel("ε")
    gcf()
    #savefig("figs/cv_afterhopf_plot.png")
end
