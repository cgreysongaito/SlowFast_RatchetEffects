include("packages.jl")
include("slowfast_commoncode.jl")


function sto_cv_plot(ep, save)
    par = RozMacPar()
    par.ε = ep
    evals = 0.441:0.005:0.9
    tspan = (0.0, 10000.0)
    tstart = 2000
    tend = 10000
    tstep = 1
    tvals = tstart:tstep:tend
    cv = fill(0.0, length(evals))
    mn = fill(0.0, length(evals))
    sd = fill(0.0, length(evals))
    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)

    for (ei, eval) in enumerate(evals)
        par.e = eval
        u0 = eq_II(par)
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
        solend = sol(tvals)
        cv[ei] = std(solend[2, 1:end] ./ mean(solend[2, 1:end]))
        mn[ei] = mean(solend[2, 1:end])
        sd[ei] = std(solend[2, 1:end])
    end

    maxcv = maximum(cv)
    maxmn = maximum(mn)
    let
        figure()
        subplot(211)
        PyPlot.plot(collect(evals), cv)
        PyPlot.vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maxcv, linestyles = "dashed")
        annotate("TC", (50, 315), xycoords = "figure points", fontsize = 12)
        annotate("R/C", (107, 315), xycoords = "figure points", fontsize = 12)
        annotate("H", (250, 315), xycoords = "figure points", fontsize = 12)
        ylabel("Consumer CV")
        subplot(212)
        PyPlot.plot(collect(evals), mn, label = "mean")
        PyPlot.plot(collect(evals), sd, label = "sd")
        PyPlot.vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maxmn, linestyles = "dashed")
        xlabel("Efficiency (e)")
        ylabel("Mean (blue), SD (orange)")
        if save == "save" && ep < 1
            savefig(joinpath(abpath(), "figs/sto_eptiny_cv_plot.png"))
        elseif save == "save"
            savefig(joinpath(abpath(), "figs/sto_ep1_cv_plot.png"))
        else
        return gcf()
        end
    end
end

sto_cv_plot(1, "save")
sto_cv_plot(1, "show")
sto_cv_plot(0.01, "save")
sto_cv_plot(0.01, "show")

function numsolvplot(u0, tspan, par, ep)
    par.ε = ep
    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return PyPlot.plot(sol.t, sol.u)
end

function stoepboth(eff)
    par = RozMacPar()
    par.e = eff
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)

    let
        figure()
        subplot(211)
        PyPlot.title("ε = 1")
        numsolvplot(u0, tspan, par, 1)
        subplot(212)
        PyPlot.title("ε = 0.01")
        numsolvplot(u0, tspan, par, 0.01)
        #return gcf()
        savefig(joinpath(abpath(), "figs/" * string(eff) * "sto.png"))
    end
end

stoepboth(0.5)

function stoepsingle(eff,ep)
    par = RozMacPar()
    par.e = eff
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)
    numsolvplot(u0, tspan, par, ep)
end

let
    figure()
    stoepsingle(0.45, 0.01)
    gcf()
end




##### Before Hopf , starting values large difference from equilibrium with small epsilon
function randstart(ep, eff)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    tspan = (0.0, 1000000.0)

    uC = [0.5, 1.0, 2.0, 2.9]
    uR = [0.5, 1.0, 2.0, 2.9]
    for i in 1:4, j in 1:4
        u0 = [uC[i], uR[j]]
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
        PyPlot.plot(sol.t, sol[2,1:end])
    end
end

let
    figure()
    randstart(0.0001,0.6)
    #gcf()
    savefig(joinpath(abpath(), "figs/eptiny_beforehopf_sol.png"))
end


par = RozMacPar()
par.ε = 0.01
par.e = 0.6
tspan = (0.0, 1000.0)

u0 = [2.9, 2.9]

prob = ODEProblem(roz_mac_II!, u0, tspan, par)
sol = DifferentialEquations.solve(prob, reltol = 1e-8, maxiters = 1e9)
let
    figure()
    PyPlot.plot(sol.t, sol[2,1:end])
    gcf()
end

#something is not working when epsilon is tiny and consumer is above the "hopf" - consumer reduces to zero and doesn't increase again and do't have enough maxiters (or takes forever)


# Create plots of con-res stochastic model before imag numbers (and before hopf)
function findRCdivide_effx(eff)
    par = RozMacPar()
    par.e = eff
    epvals = 0.00001:0.000001:1

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        equ = eq_II(par)
        eig1 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        if eig1 < 0 || eig1 > 0
            return epval
            break
        end
    end
end

findRCdivide_effx(0.55)

function noiseACF_plot(eff, ep, sto)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    eq = eq_II(par)
    u0 = randeq.(eq)
    tspan = (0.0, 10000.0)
    tstart = 9000
    tend = 10000
    tstep = 1
    tvals = tstart:tstep:tend

    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    if sto == "yes"
        sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    else
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    end
    endsol = sol(tvals)
    PyPlot.plot(endsol[1,1:end], endsol[2,1:end])
    return xlim(0.0,2.8)
end

let
    figure(figsize = (7,10))
    subplots_adjust(right = 0.75)
    subplot(5,1,1)
    noiseACF_plot(0.55, 0.1, "yes")
    title("e = 0.55", fontsize = 15)
    annotate("ε = 0.1", (400, 580), xycoords = "figure points", fontsize = 15)
    subplot(5,1,2)
    noiseACF_plot(0.55, 0.4007, "yes")
    annotate("ε = 0.4007", (400, 460), xycoords = "figure points", fontsize = 15)
    ylabel("Consumer", fontsize = 15)
    subplot(5,1,3)
    noiseACF_plot(0.55, 0.4008, "yes")
    annotate("ε = 0.4008", (400, 350), xycoords = "figure points", fontsize = 15)
    subplot(5,1,4)
    noiseACF_plot(0.55, 0.9, "yes")
    annotate("ε = 0.9", (400, 240), xycoords = "figure points", fontsize = 15)
    subplot(5,1,5)
    noiseACF_plot(1.0, 0.01, "det")
    annotate("e = 1.0\nε = 0.01\nDeterministic", (400, 120), xycoords = "figure points", fontsize = 15)
    xlabel("Resource", fontsize = 15)
    # gcf()
    savefig(joinpath(abpath(), "figs/noiseACF_effbeforehopf_phase.png"))
end
# are we flipping where ACF and white noise should be found - looks like white noise found when eigenvalues have complex - something seems wrong

let
    figure(figsize = (7,10))
    subplots_adjust(right = 0.75)
    subplot(5,1,1)
    noiseACF_plot(0.46, 0.01, "yes")
    title("ε = 0.01", fontsize = 15)
    annotate("e = 0.46", (400, 580), xycoords = "figure points", fontsize = 15)
    subplot(5,1,2)
    noiseACF_plot(0.52, 0.01, "yes")
    annotate("e = 0.52", (400, 460), xycoords = "figure points", fontsize = 15)
    ylabel("Consumer", fontsize = 15)
    subplot(5,1,3)
    noiseACF_plot(0.53, 0.01, "yes")
    annotate("e = 0.53", (400, 350), xycoords = "figure points", fontsize = 15)
    subplot(5,1,4)
    noiseACF_plot(0.71, 0.01, "yes")
    annotate("e = 0.71", (400, 240), xycoords = "figure points", fontsize = 15)
    subplot(5,1,5)
    noiseACF_plot(1.0, 0.01, "det")
    annotate("e = 1.0\nε = 0.01\nDeterministic", (400, 120), xycoords = "figure points", fontsize = 15)
    xlabel("Resource", fontsize = 15)
    gcf()
    #savefig(joinpath(abpath(), "figs/noiseACF_effbeforehopf_phase.png"))
end


##### Autocorrelation analysis as efficiency changes with tiny epsilon - in stochastic

function acf_plot(ep, eff, sto)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)
    tstart = 4000
    tend = 10000
    tstep = 1
    tvals = tstart:tstep:tend

    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)

    if sto == "yes"
        sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    else
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    end

    solend = sol(tvals)
    if ep == 1.0
        lrange = 0:1:200
    else
        lrange = 0:1:2000
    end
    acf = autocor(solend[2, 1:end], collect(lrange))

    PyPlot.bar(collect(lrange), acf)
    ylim(-1,1)
    return ylabel("ACF")
end

let
    figure(figsize = (8,10))
    subplots_adjust(hspace = 0.2, wspace = 0.3)
    subplot(5,2,1)
    acf_plot(1.0,0.45, "yes")
    ax1 = gca()
    setp(ax1.get_xticklabels(),visible=false)
    subplot(5,2,2)
    acf_plot(0.01,0.45, "yes")
    ax5 = gca()
    setp(ax5.get_xticklabels(),visible=false)
    subplot(5,2,3)
    acf_plot(1.0,0.52, "yes")
    ax2 = gca()
    setp(ax2.get_xticklabels(),visible=false)
    subplot(5,2,4)
    acf_plot(0.01,0.52, "yes")
    ax6 = gca()
    setp(ax6.get_xticklabels(),visible=false)
    subplot(5,2,5)
    acf_plot(1.0, 0.53, "yes")
    ax3 = gca()
    setp(ax3.get_xticklabels(),visible=false)
    subplot(5,2,6)
    acf_plot(0.01, 0.53, "yes")
    ax7 = gca()
    setp(ax7.get_xticklabels(),visible=false)
    subplot(5,2,7)
    acf_plot(1.0, 0.71, "yes")
    ax4 = gca()
    setp(ax4.get_xticklabels(),visible=false)
    subplot(5,2,8)
    acf_plot(0.01, 0.71, "yes")
    ax8 = gca()
    setp(ax8.get_xticklabels(),visible=false)
    subplot(5,2,9)
    acf_plot(1.0, 1.0, "no")
    xlabel("Lag")
    subplot(5,2,10)
    acf_plot(0.01, 1.0, "no")
    xlabel("Lag")
    annotate("e = 0.45\nStochastic", (515, 575), xycoords = "figure points", fontsize = 12)
    annotate("e = 0.52\nStochastic", (515, 450), xycoords = "figure points", fontsize = 12)
    annotate("e = 0.53\nStochastic", (515, 350), xycoords = "figure points", fontsize = 12)
    annotate("e = 0.71\nStochastic", (515, 225), xycoords = "figure points", fontsize = 12)
    annotate("e = 1.0\nDeterministic", (515, 125), xycoords = "figure points", fontsize = 12)
    annotate("ε = 1.0", (100, 650), xycoords = "figure points", fontsize = 12)
    annotate("ε = 0.01", (400, 650), xycoords = "figure points", fontsize = 12)
    annotate("R/C", (515, 405), xycoords = "figure points", fontsize = 12)
    annotate("Hopf", (515, 180), xycoords = "figure points", fontsize = 12)
    # gcf()
    savefig(joinpath(abpath(), "figs/ACF.png"))
end
