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

# cb = PeriodicCallback(pert_cb, 1, initial_affect = true)

function numsolvplot(u0, tspan, par, ep, mean)
    par.ε = ep
    par.μ = mean
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

    sto_epboth_fig = figure()
    subplot(211)
    PyPlot.title("ε = 1")
    numsolvplot(u0, tspan, par, 1, 0.001)
    subplot(212)
    PyPlot.title("ε = 0.01")
    numsolvplot(u0, tspan, par, 0.01, 0.001)
    return sto_epboth_fig
    #savefig(joinpath(abpath(), "figs/" * string(eff) * "sto.png"))
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

function noiseACF_plot(eff, ep, sto, seed, mean)
    Random.seed!(seed)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    par.μ = mean
    eq = eq_II(par)
    tspan = (0.0, 10000.0)
    tstart = 1
    tend = 10000
    tstep = 1
    tvals = tstart:tstep:tend

    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
    prob = ODEProblem(roz_mac_II!, eq, tspan, par)
    if sto == "yes"
        sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    else
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    end
    endsol = sol(tvals)
    PyPlot.plot(endsol[1,1:end], endsol[2,1:end])
    iso_plot(range(0, stop = 3, length = 100), par, eff)
    ylim(0.0,3.0)
    return xlim(0.0,3.0)
end

let
    seed = 2
    noiseACF = figure(figsize = (5,12))
    subplot(5,1,1)
    noiseACF_plot(0.55, 0.1, "yes", seed)
    title("(A) e = 0.55, ε = 0.1", fontsize = 15)
    subplot(5,1,2)
    title("(B) e = 0.55, ε = 0.4007", fontsize = 15)
    noiseACF_plot(0.55, 0.4007, "yes", seed)
    ylabel("Consumer", fontsize = 15)
    subplot(5,1,3)
    title("(C) e = 0.55, ε = 0.4008", fontsize = 15)
    noiseACF_plot(0.55, 0.4008, "yes", seed)
    subplot(5,1,4)
    title("(D) e = 0.55, ε = 0.9", fontsize = 15)
    noiseACF_plot(0.55, 0.9, "yes", seed)
    subplot(5,1,5)
    title("(E) e = 1.0, ε = 0.01, Deterministic", fontsize = 15)
    noiseACF_plot(1.0, 0.01, "det", seed)
    xlabel("Resource", fontsize = 15)
    tight_layout()
    return noiseACF
    #savefig(joinpath(abpath(), "figs/noiseACF_effbeforehopf_phase.png"))
end
# are we flipping where ACF and white noise should be found - looks like white noise found when eigenvalues have complex - something seems wrong

let #examining e = 0.52 and ep = 0.01 seeing whether close to zero
    par = RozMacPar()
    par.ε = 0.01
    par.e = 0.52
    test = figure(figsize = (8, 4))
    subplot(1,2,1)
    title("(A) ACF ε = 0.01, e = 0.52, Sto")
    acf_plot(0.01,0.52, "yes", 3)
    subplot(1,2,2)
    title("(B) Phase space ε = 0.01, e = 0.52, Sto")
    noiseACF_plot(0.52, 0.01, "yes", 3, 0.0)
    iso_plot(range(0, stop = 3, length = 100), par, 0.52)
    xlabel("Resource")
    ylabel("Consumer")
    tight_layout()
    # return test
    savefig(joinpath(abpath(), "figs/e52_ep001_isoclines_comparisonwithACF.png"))
end

let
    seed = 2
    mean = 0.001
    sto_canard = figure(figsize = (5,12))
    subplot(5,1,1)
    title("(A) ε = 0.01, e = 0.46", fontsize = 15)
    noiseACF_plot(0.46, 0.01, "yes", seed, mean)
    subplot(5,1,2)
    title("(B) ε = 0.01, e = 0.52", fontsize = 15)
    noiseACF_plot(0.52, 0.01, "yes", seed, mean)
    ylabel("Consumer", fontsize = 15)
    subplot(5,1,3)
    title("(C) ε = 0.01, e = 0.53", fontsize = 15)
    noiseACF_plot(0.53, 0.01, "yes", seed, mean)
    subplot(5,1,4)
    title("(D) ε = 0.01, e = 0.71", fontsize = 15)
    noiseACF_plot(0.71, 0.01, "yes", seed, mean)
    subplot(5,1,5)
    title("(E) ε = 0.01, e = 1.0, Deterministic", fontsize = 15)
    noiseACF_plot(1.0, 0.01, "det", seed, mean)
    xlabel("Resource", fontsize = 15)
    tight_layout()
    return sto_canard
    #savefig(joinpath(abpath(), "figs/noiseACF_effbeforehopf_phase.png"))
end


##### Autocorrelation analysis as efficiency changes with tiny epsilon - in stochastic

function acf_plot(ep, eff, sto, seed)
    Random.seed!(seed)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)
    tstart = 6000
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
        lrange = 0:1:50
    else
        lrange = 0:1:500
    end
    acf = autocor(solend[2, 1:end], collect(lrange))
    conf = 1.96/sqrt(length(solend))
    PyPlot.bar(collect(lrange), acf)
    hlines(0 + conf, 0, maximum(lrange))
    hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    return ylabel("ACF")
end


Random.seed!(seed)
par = RozMacPar()
par.ε = 0.01
par.e = 0.45
u0 = eq_II(par)
tspan = (0.0, 10000.0)
tstart = 6000
tend = 10000
tstep = 1
tvals = tstart:tstep:tend

cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
prob = ODEProblem(roz_mac_II!, u0, tspan, par)

sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)

solend = sol(tvals)
acf = autocor(solend[2, 1:end], collect(0:1:50))
conf = 1.96/sqrt(length(solend))
let
    test = figure()
    subplot(1,2,1)
    plot(solend.t, solend.u)
    subplot(1,2,2)
    bar(collect(0:1:50), acf)
    hlines(0 + conf, 0, maximum(0:1:50))
    hlines(0 - conf, 0, maximum(0:1:50))
    return test
end

let
    seed = 3
    acfplot = figure(figsize = (8,12))
    subplot(5,2,1)
    title("(A) ε = 1.0, e = 0.45, Sto")
    acf_plot(1.0,0.45, "yes", seed)
    subplot(5,2,2)
    title("(B) ε = 0.01, e = 0.45, Sto")
    acf_plot(0.01,0.45, "yes", seed)
    subplot(5,2,3)
    title("(C) ε = 1.0, e = 0.52, Sto")
    acf_plot(1.0,0.52, "yes", seed)
    subplot(5,2,4)
    title("(D) ε = 0.01, e = 0.52, Sto")
    acf_plot(0.01,0.52, "yes", seed)
    subplot(5,2,5)
    title("(E) ε = 1.0, e = 0.53, Sto")
    acf_plot(1.0, 0.53, "yes", seed)
    subplot(5,2,6)
    title("(F) ε = 0.01, e = 0.53, Sto")
    acf_plot(0.01, 0.53, "yes", seed)
    subplot(5,2,7)
    title("(G) ε = 1.0, e = 0.71, Sto")
    acf_plot(1.0, 0.71, "yes", seed)
    subplot(5,2,8)
    title("(H) ε = 0.01, e = 0.71, Sto")
    acf_plot(0.01, 0.71, "yes", seed)
    subplot(5,2,9)
    title("(I) ε = 1.0, e = 1.0, Det")
    acf_plot(1.0, 1.0, "no", seed)
    subplot(5,2,10)
    title("(J) ε = 0.01, e = 1.0, Det")
    acf_plot(0.01, 1.0, "no", seed)
    # annotate("R/C", (515, 405), xycoords = "figure points", fontsize = 12)
    # annotate("Hopf", (515, 180), xycoords = "figure points", fontsize = 12)
    tight_layout()
    return acfplot
    # savefig(joinpath(abpath(), "figs/ACF.png"))
end


# Testing whether high versus low frequency pert increases probablity of canard
# figure can change either efficiency and or frequency - four figures for four values of efficiency , x axis is frequency y axis is proportion

function canard_proportion(eff, ep, repeat, mean)
    freq_vals = 1.0:0.5:10.0
    canard_prop = fill(0.0, length(freq_vals))

    for (freq_i, freq_val) in enumerate(freq_vals)
        par = RozMacPar()
        par.ε = ep
        par.e = eff
        par.μ = mean
        eq = eq_II(par)
        tspan = (0.0, 10000.0)
        tstart = 6000
        tend = 10000
        tstep = 1
        tvals = tstart:tstep:tend
        canard_count = 0
        for j in 1:repeat
            cb = PeriodicCallback(pert_cb, freq_val, initial_affect = true)
            prob = ODEProblem(roz_mac_II!, eq, tspan, par)
            sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
            endsol = sol(tvals)

            for k in 1:length(endsol)
                if 0 < endsol.u[j][2] < 1.6 && 0 < endsol.u[j][1] < 1.6
                    canard_count += 1
                    break
                end
            end
        end
        canard_prop[freq_i] = canard_count / repeat
    end
    return canard_prop
end

let
    canard_prop_plot = figure(figsize = (7, 10))
    subplot(3,1,1)
    title("(A) Perturbation μ = 0.0")
    scatter(1.0:0.5:10.0, canard_proportion(0.46, 0.01, 100, 0.0), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(0.53, 0.01, 100, 0.0), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(0.6, 0.01, 100, 0.0), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(0.71, 0.01, 100, 0.0), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    subplot(3,1,2)
    title("(B) Perturbation μ = 0.0001")
    scatter(1.0:0.5:10.0, canard_proportion(0.46, 0.01, 100, 0.001), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(0.53, 0.01, 100, 0.001), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(0.6, 0.01, 100, 0.001), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(0.71, 0.01, 100, 0.001), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    subplot(3,1,3)
    title("(C) Perturbation μ = 0.0002")
    scatter(1.0:0.5:10.0, canard_proportion(0.46, 0.01, 100, 0.002), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(0.53, 0.01, 100, 0.002), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(0.6, 0.01, 100, 0.002), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(0.71, 0.01, 100, 0.002), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    tight_layout()
    # return canard_prop_plot
    savefig(joinpath(abpath(), "figs/perturbation_canard_proprtion.png"))
end
