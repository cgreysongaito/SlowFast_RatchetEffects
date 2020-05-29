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



function pert_cv_plot(ep)
    evals = 0.441:0.005:0.9
    cv = fill(0.0, length(evals))

    for (ei, eval) in enumerate(evals)
        sol = RozMac_pert(ep, eval, 0.0, 3, 10000.0, 6000.0:1.0:10000.0)
        cv[ei] = std(sol[2, :]) / mean(sol[2, :])
    end

    plot(collect(evals), cv)
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = 0.4, linestyles = "dashed")
    ylabel("Consumer CV")
    return xlabel("Efficiency (e)")
end

function pert_mean_plot(ep)
    evals = 0.441:0.005:0.9
    mn = fill(0.0, length(evals))
    sd = fill(0.0, length(evals))

    for (ei, eval) in enumerate(evals)
        sol = RozMac_pert(ep, eval, 0.0, 3, 10000.0, 6000.0:1.0:10000.0)
        sd[ei] = std(sol[2, :])
        mn[ei] = mean(sol[2, :])
    end

    plot(collect(evals), mn, label = "mean")
    plot(collect(evals), sd, label = "sd")
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maximum(mn), linestyles = "dashed")
    ylabel("Mean (blue), SD (orange)")
    return xlabel("Efficiency (e)")
end

let
    sto_ep1_cv_plot = figure()
    subplot(2,1,1)
    pert_cv_plot(1.0)
    subplot(2,1,2)
    pert_mean_plot(1.0)
    # annotate("TC", (50, 315), xycoords = "figure points", fontsize = 12)
    # annotate("R/C", (107, 315), xycoords = "figure points", fontsize = 12)
    # annotate("H", (250, 315), xycoords = "figure points", fontsize = 12)
    tight_layout()
    return sto_ep1_cv_plot
    # savefig(joinpath(abpath(), "figs/sto_ep1_cv_plot.png"))
end

let
    sto_eptiny_cv_plot = figure()
    subplot(2,1,1)
    pert_cv_plot(0.01)
    subplot(2,1,2)
    pert_mean_plot(0.01)
    # annotate("TC", (50, 315), xycoords = "figure points", fontsize = 12)
    # annotate("R/C", (107, 315), xycoords = "figure points", fontsize = 12)
    # annotate("H", (250, 315), xycoords = "figure points", fontsize = 12)
    tight_layout()
    return sto_eptiny_cv_plot
    # savefig(joinpath(abpath(), "figs/sto_eptiny_cv_plot.png"))
end


function pert_con_minmax_plot(ep, stand)
    evals = 0.441:0.005:0.9
    min_con = fill(0.0, length(evals))
    max_con = fill(0.0, length(evals))

    for (ei, eval) in enumerate(evals)
        sol = RozMac_pert(ep, eval, 0.0, 3, 10000.0, 6000.0:1.0:10000.0)
        if stand == "standardized"
            min_con[ei] = minimum(sol[2,:]) / eq_II(RozMacPar(e = eval))[2]
            max_con[ei] = maximum(sol[2,:]) / eq_II(RozMacPar(e = eval))[2]
        else
            min_con[ei] = minimum(sol[2,:])
            max_con[ei] = maximum(sol[2,:])
        end
    end

    scatter(collect(evals), min_con)
    scatter(collect(evals), max_con)
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maximum(max_con), linestyles = "dashed")
    return xlabel("Efficiency (e)")
end

let
    conminmax_plot = figure()
    subplot(2,1,1)
    pert_con_minmax_plot(0.01, "no")
    ylabel("Consumer Min/Max")
    subplot(2,1,2)
    pert_con_minmax_plot(0.01, "standardized")
    ylabel("Consumer Min/Max \n Standardized by equilibrium values")
    # annotate("TC", (50, 350), xycoords = "figure points", fontsize = 12)
    # annotate("R/C", (107, 350), xycoords = "figure points", fontsize = 12)
    # annotate("H", (250, 350), xycoords = "figure points", fontsize = 12)
    tight_layout()
    return conminmax_plot
end

let
    sto_epboth_fig = figure()
    subplot(2,1,1)
    title("ε = 1")
    pert_timeseries_plot(1.0, 0.5, 0.0, 3, 10000.0, 6000.0:1.0:10000.0)
    subplot(2,1,2)
    title("ε = 0.01")
    pert_timeseries_plot(0.01, 0.5, 0.0, 3, 10000.0, 6000.0:1.0:10000.0)
    tight_layout()
    return sto_epboth_fig
    #savefig(joinpath(abpath(), "figs/" * string(eff) * "sto.png"))
end


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

function phase_epgradient_plot(ep, eff, mean, seed, tsend, tvals)
    pert_phase_plot(ep, eff, mean, seed, tsend, tvals)
    iso_plot(range(0, stop = 3, length = 100), RozMacPar(e = eff))
    xlabel("Resource")
    ylabel("Consumer")
    ylim(0.0,3.0)
    return xlim(0.0,3.0)
end

let
    mean = 0.0
    seed = 2
    tsend = 10000.0
    tvals = 6000.0:1.0:10000.0
    epgradient = figure(figsize = (8,12))
    subplot(4,2,1)
    title("(A) e = 0.55, ε = 0.1", fontsize = 15)
    pert_timeseries_plot(0.1, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,2)
    title("(B) e = 0.55, ε = 0.1", fontsize = 15)
    phase_epgradient_plot(0.1, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,3)
    title("(C) e = 0.55, ε = 0.4007", fontsize = 15)
    pert_timeseries_plot(0.4007, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,4)
    title("(D) e = 0.55, ε = 0.4007", fontsize = 15)
    phase_epgradient_plot(0.4007, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,5)
    title("(E) e = 0.55, ε = 0.4008", fontsize = 15)
    pert_timeseries_plot(0.4008, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,6)
    title("(F) e = 0.55, ε = 0.4008", fontsize = 15)
    phase_epgradient_plot(0.4008, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,7)
    title("(G) e = 0.55, ε = 0.9", fontsize = 15)
    pert_timeseries_plot(0.9, 0.55, mean, seed, tsend, tvals)
    subplot(4,2,8)
    title("(H) e = 0.55, ε = 0.9", fontsize = 15)
    phase_epgradient_plot(0.9, 0.55, mean, seed, tsend, tvals)
    tight_layout()
    return epgradient
    #savefig(joinpath(abpath(), "figs/noiseACF_effbeforehopf_phase.png"))
end
# are we flipping where ACF and white noise should be found - looks like white noise found when eigenvalues have complex - something seems wrong

let
    mean = 0.0
    seed = 2
    tsend = 10000.0
    tvals = 1000.0:1.0:10000.0
    sto_canard = figure(figsize = (5,12))
    subplot(4,1,1)
    title("(A) ε = 0.01, e = 0.46", fontsize = 15)
    phase_epgradient_plot(0.01, 0.46, mean, seed, tsend, tvals)
    subplot(4,1,2)
    title("(B) ε = 0.01, e = 0.52", fontsize = 15)
    phase_epgradient_plot(0.01, 0.52, mean, seed, tsend, tvals)
    subplot(4,1,3)
    title("(C) ε = 0.01, e = 0.53", fontsize = 15)
    phase_epgradient_plot(0.01, 0.53, mean, seed, tsend, tvals)
    subplot(4,1,4)
    title("(D) ε = 0.01, e = 0.71", fontsize = 15)
    phase_epgradient_plot(0.01, 0.71, mean, seed, tsend, tvals)
    tight_layout()
    return sto_canard
    #savefig(joinpath(abpath(), "figs/noiseACF_effbeforehopf_phase.png"))
end


##### Autocorrelation analysis as efficiency changes with tiny epsilon - in stochastic


function acf_plot(lrange, ep, eff, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, seed, tsend, tvals)
    acf = autocor(sol[2, :], collect(lrange))
    conf = 1.96/sqrt(length(sol))
    PyPlot.bar(collect(lrange), acf)
    hlines(0 + conf, 0, maximum(lrange))
    hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    return ylabel("ACF")
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

let
    seed = 3
    acfplot_ep1 = figure(figsize = (12,12))
    subplot(4,3,1)
    title("(A) ε = 1.0, e = 0.45, Sto")
    acf_plot(1.0,0.45, "yes", seed)
    subplot(4,3,2)
    title("(B) ε = 1.0, e = 0.45, Sto")
    timeseries(1.0, 0.45, seed)
    subplot(4,3,3)
    title("(C) ε = 1.0, e = 0.45, Sto")
    phase(1.0, 0.45, seed)
    subplot(4,3,4)
    title("(D) ε = 1.0, e = 0.52, Sto")
    acf_plot(1.0, 0.52, "yes", seed)
    subplot(4,3,5)
    title("(E) ε = 1.0, e = 0.52, Sto")
    timeseries(1.0, 0.52, seed)
    subplot(4,3,6)
    title("(F) ε = 1.0, e = 0.52, Sto")
    phase(1.0, 0.52, seed)
    subplot(4,3,7)
    title("(G) ε = 1.0, e = 0.53, Sto")
    acf_plot(1.0, 0.53, "yes", seed)
    subplot(4,3,8)
    title("(H) ε = 1.0, e = 0.53, Sto")
    timeseries(1.0, 0.53, seed)
    subplot(4,3,9)
    title("(I) ε = 1.0, e = 0.53, Sto")
    phase(1.0, 0.53, seed)
    subplot(4,3,10)
    title("(J) ε = 1.0, e = 0.71, Sto")
    acf_plot(1.0, 0.71, "yes", seed)
    subplot(4,3,11)
    title("(K) ε = 1.0, e = 0.71, Sto")
    timeseries(1.0, 0.71, seed)
    subplot(4,3,12)
    title("(L) ε = 1.0, e = 0.71, Sto")
    phase(1.0, 0.71, seed)
    tight_layout()
    # return acfplot_ep1
    savefig(joinpath(abpath(), "figs/ACFplot_ep1.png"))
end

let
    seed = 3
    tsend = 10000.0
    tvals = 5000.0:100.0:10000.0
    lrange = 0:1:50
    acfplot_ep001 = figure(figsize = (12,12))
    subplot(4,3,1)
    title("(A) ε = 0.01, e = 0.45, Sto")
    acf_plot(lrange, 0.01, 0.45, seed, tsend, tvals)
    subplot(4,3,2)
    title("(B) ε = 0.01, e = 0.45, Sto")
    pert_timeseries(0.01, 0.45, seed, tsend, tvals)
    subplot(4,3,3)
    title("(C) ε = 0.01, e = 0.45, Sto")
    pert_phase(0.01, 0.45, seed, tsend, tvals)
    subplot(4,3,4)
    title("(D) ε = 0.01, e = 0.52, Sto")
    acf_plot(lrange, 0.01, 0.52, seed, tsend, tvals)
    subplot(4,3,5)
    title("(E) ε = 0.01, e = 0.52, Sto")
    pert_timeseries(0.01, 0.52, seed, tsend, tvals)
    subplot(4,3,6)
    title("(F) ε = 0.01, e = 0.52, Sto")
    pert_phase(0.01, 0.52, seed, tsend, tvals)
    subplot(4,3,7)
    title("(G) ε = 0.01, e = 0.53, Sto")
    acf_plot(lrange, 0.01, 0.53, seed, tsend, tvals)
    subplot(4,3,8)
    title("(H) ε = 0.01, e = 0.53, Sto")
    pert_timeseries(0.01, 0.53, seed, tsend, tvals)
    subplot(4,3,9)
    title("(I) ε = 0.01, e = 0.53, Sto")
    pert_phase(0.01, 0.53, seed, tsend, tvals)
    subplot(4,3,10)
    title("(J) ε = 0.01, e = 0.71, Sto")
    acf_plot(lrange, 0.01, 0.71, seed, tsend, tvals)
    subplot(4,3,11)
    title("(K) ε = 0.01, e = 0.71, Sto")
    pert_timeseries(0.01, 0.71, seed, tsend, tvals)
    subplot(4,3,12)
    title("(L) ε = 0.01, e = 0.71, Sto")
    pert_phase(0.01, 0.71, seed, tsend, tvals)
    tight_layout()
    return acfplot_ep001
    # savefig(joinpath(abpath(), "figs/ACFplot_ep001.png"))
end

# Testing whether high versus low frequency pert increases probablity of canard
# figure can change either efficiency and or frequency - four figures for four values of efficiency , x axis is frequency y axis is proportion

function canard_proportion(repeat, ep, eff, mean, tsend, tvals)
    freq_vals = 1.0:0.5:10.0
    canard_prop = fill(0.0, length(freq_vals))

    for (freq_i, freq_val) in enumerate(freq_vals)
        canard_count = 0
        for j in 1:repeat
            sol = RozMac_pert(ep, eff, mean, rand(1:100000), tsend, tvals)

            for k in 1:length(sol)
                if 0 < sol.u[j][2] < 1.6 && 0 < sol.u[j][1] < 1.6
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
    tsend = 10000.0
    tvals = 6000.0:1.0:10000.0
    canard_prop_plot = figure(figsize = (7, 10))
    subplot(3,1,1)
    title("(A) Perturbation μ = 0.0")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.46, 0.0, tsend, tvals), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.53, 0.0, tsend, tvals), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.6, 0.0, tsend, tvals), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.71, 0.0, tsend, tvals), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    subplot(3,1,2)
    title("(B) Perturbation μ = 0.0001")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.46, 0.001, tsend, tvals), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.53, 0.001, tsend, tvals), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.6, 0.001, tsend, tvals), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.71, 0.001, tsend, tvals), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    subplot(3,1,3)
    title("(C) Perturbation μ = 0.0002")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.46, 0.002, tsend, tvals), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.53, 0.001, tsend, tvals), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.6, 0.001, tsend, tvals), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.71, 0.001, tsend, tvals), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    tight_layout()
    return canard_prop_plot
    # savefig(joinpath(abpath(), "figs/perturbation_canard_proprtion.png"))
end



#### Old code - not sure if want
# let #examining e = 0.52 and ep = 0.01 seeing whether close to zero
#     par = RozMacPar()
#     par.ε = 0.01
#     par.e = 0.52
#     test = figure(figsize = (8, 4))
#     subplot(1,2,1)
#     title("(A) ACF ε = 0.01, e = 0.52, Sto")
#     acf_plot(0.01,0.52, "yes", 3)
#     subplot(1,2,2)
#     title("(B) Phase space ε = 0.01, e = 0.52, Sto")
#     noiseACF_plot(0.52, 0.01, "yes", 3, 0.0)
#     iso_plot(range(0, stop = 3, length = 100), par, 0.52)
#     xlabel("Resource")
#     ylabel("Consumer")
#     tight_layout()
#     # return test
#     savefig(joinpath(abpath(), "figs/e52_ep001_isoclines_comparisonwithACF.png"))
# end



##### Before Hopf , starting values large difference from equilibrium with small epsilon
# function randstart(ep, eff)
#     par = RozMacPar()
#     par.ε = ep
#     par.e = eff
#     tspan = (0.0, 1000000.0)
#
#     uC = [0.5, 1.0, 2.0, 2.9]
#     uR = [0.5, 1.0, 2.0, 2.9]
#     for i in 1:4, j in 1:4
#         u0 = [uC[i], uR[j]]
#         prob = ODEProblem(roz_mac_II!, u0, tspan, par)
#         sol = DifferentialEquations.solve(prob, reltol = 1e-8)
#         PyPlot.plot(sol.t, sol[2,1:end])
#     end
# end
#
# let
#     figure()
#     randstart(0.0001,0.6)
#     #gcf()
#     savefig(joinpath(abpath(), "figs/eptiny_beforehopf_sol.png"))
# end

#something is not working when epsilon is tiny and consumer is above the "hopf" - consumer reduces to zero and doesn't increase again and do't have enough maxiters (or takes forever)
# due to stiffness of problem
