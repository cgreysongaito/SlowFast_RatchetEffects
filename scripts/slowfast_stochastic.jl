include("packages.jl")
include("slowfast_commoncode.jl")
include("slowfast_canardfinder.jl")


## Figure 2 (proportion quasi-canard with different consumer isoclines)
function prop_canard(ep, eff)
    res_start = 0.0:0.1:3.0
    con_start = 0.0:0.1:2.5
    count_true = 0
    for (resi, resvalue) in enumerate(res_start)
        for (coni, convalue) in enumerate(con_start)
            if cf_returnmap(ep, eff, 1, 0.0, 1234, resvalue, convalue, 5000.0, 2000.0:1.0:5000.0) == true
                count_true += 1
            end
        end
    end
    return count_true / (length(res_start)*length(con_start))
end

cf_returnmap(0.1, 0.6, 1, 0.0, 1234, 0.1, 0.1, 5000.0, 2000.0:1.0:5000.0)

prop_canard(0.01,0.6)

let
    test = figure()
    pert_phase_plot(0.1, 0.6, 1, 0.0, rand, 0.1, 0.1, 5000.0, 2000.0:1.0:5000.0)
    return test
end

# pert_phase_plot(0.1, 0.6, 1, 0.0, 7, 0.1, 0.1, 5000.0, 2000.0:1.0:5000.0) - why shows nothing

function prop_canard_data(eff)
    epsilon_range = 0.0:0.01:1.0 #shouldn't include 0.0
    propcanard = fill(0.0,length(epsilon_range))
    for (epi, epvalue) in enumerate(epsilon_range)
        propcanard[epi] = prop_canard(epvalue, eff)
    end
    return propcanard
end

fulldataset = [prop_canard_data(0.5), prop_canard_data(0.6), prop_canard_data(0.7), prop_canard_data(0.8)]
fulldataset[1]
let
    figure2 = figure()
    subplot(1,4,1)
    plot(collect(0.0:0.1:1.0), fulldataset[1])
    subplot(1,4,2)
    plot(collect(0.0:0.1:1.0), fulldataset[2])
    subplot(1,4,3)
    plot(collect(0.0:0.1:1.0), fulldataset[3])
    subplot(1,4,4)
    plot(collect(0.0:0.1:1.0), fulldataset[4])
    return figure2
end





function pert_cvsdmean(ep)
    evals = 0.441:0.005:0.9
    cv = fill(0.0, length(evals))
    mn = fill(0.0, length(evals), 1)
    sd = fill(0.0, length(evals), 1)

    for (ei, eval) in enumerate(evals)
        sol = RozMac_pert(ep, eval, 0.0, 1, 3, 10000.0, 6000.0:1.0:10000.0)
        cv[ei] = std(sol[2, :]) / mean(sol[2, :])
        sd[ei] = std(sol[2, :])
        mn[ei] = mean(sol[2, :])
    end
    cv[cv .==0.0] .= NaN
    sd[sd .==0.0] .= NaN
    mn[mn .==0.0] .= NaN
    return hcat(evals, cv, sd, mn)
end

let
    data = pert_cvsdmean(0.01)
    sto_eptiny_cv_plot = figure()
    subplot(3,1,1)
    plot(data[:, 1], data[:, 2])
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = 0.4, linestyles = "dashed")
    ylabel("Coefficient of \n Variation (CV)")
    xlabel("Efficiency (e)")

    subplot(3,1,2)
    plot(data[:,1], data[:, 4])
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maximum(filter(!isnan, data[:,4])), linestyles = "dashed")
    ylabel("Mean (μ)")
    xlabel("Efficiency (e)")

    subplot(3,1,3)
    plot(data[:,1], data[:, 3])
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maximum(filter(!isnan, data[:,4])), linestyles = "dashed")
    ylabel("Standard \n Deviation (σ) ")
    xlabel("Efficiency (e)")
    ylim(0.0, 0.35)
    annotate("TC", (55, 315), xycoords = "figure points", fontsize = 12)
    annotate("R/C", (112, 315), xycoords = "figure points", fontsize = 12)
    annotate("H", (255, 315), xycoords = "figure points", fontsize = 12)
    tight_layout()
    # return sto_eptiny_cv_plot
    savefig(joinpath(abpath(), "figs/sto_eptiny_cv_plot.png"))
end


function pert_con_minmax(ep)
    evals = 0.441:0.005:0.9
    min_con = fill(0.0, length(evals))
    max_con = fill(0.0, length(evals))
    min_con_stand = fill(0.0, length(evals))
    max_con_stand = fill(0.0, length(evals))

    for (ei, eval) in enumerate(evals)
        sol = RozMac_pert(ep, eval, 0.0, 3, 10000.0, 6000.0:1.0:10000.0)
        min_con[ei] = minimum(sol[2,:])
        max_con[ei] = maximum(sol[2,:])
        min_con_stand[ei] = minimum(sol[2,:]) / eq_II(RozMacPar(e = eval))[2]
        max_con_stand[ei] = maximum(sol[2,:]) / eq_II(RozMacPar(e = eval))[2]
    end

    min_con[min_con .== 0.0] .= NaN
    max_con[max_con .== 0.0] .= NaN
    min_con_stand[min_con_stand .== 0.0] .= NaN
    max_con_stand[max_con_stand .== 0.0] .= NaN
    return hcat(evals, min_con, max_con, min_con_stand, max_con_stand)
end

let
    data = pert_con_minmax(0.01)
    conminmax_plot = figure()
    subplot(2,1,1)
    scatter(data[:,1], data[:, 2])
    scatter(data[:,1], data[:, 3])
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maximum(filter(!isnan, data[:,3])), linestyles = "dashed")
    xlabel("Efficiency (e)")
    ylabel("Consumer Min/Max")

    subplot(2,1,2)
    scatter(data[:,1], data[:, 4])
    scatter(data[:,1], data[:, 5])
    vlines([0.441,0.5225, 0.710], ymin = 0.0, ymax = maximum(filter(!isnan, data[:,5])), linestyles = "dashed")
    xlabel("Efficiency (e)")
    ylabel("Consumer Min/Max \n Standardized by equilibrium values")
    # annotate("TC", (50, 350), xycoords = "figure points", fontsize = 12)
    # annotate("R/C", (107, 350), xycoords = "figure points", fontsize = 12)
    # annotate("H", (250, 350), xycoords = "figure points", fontsize = 12)
    tight_layout()
    # return conminmax_plot
    savefig(joinpath(abpath(), "figs/consumer_minmax_pert.png"))
end


# Create plots of con-res stochastic model before/after imag numbers (and before hopf)
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

# Figure showing stochastic canards can be created before the Hopf when ep = 0.01
let
    mean = 0.0
    seed = 3
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


function acf_plot(lrange, ep, eff, mean, freq, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, mean, freq, seed, tsend, tvals)
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
    mean = 0.0
    freq = 1
    seed = 3
    tsend = 10000.0
    tvals = 5000.0:100.0:10000.0
    lrange = 0:1:50
    acfplot_ep001 = figure(figsize = (12,12))
    subplot(4,3,1)
    title("(A) ε = 0.01, e = 0.46, Sto")
    acf_plot(lrange, 0.01, 0.46, mean, freq, seed, tsend, tvals)
    subplot(4,3,2)
    title("(B) ε = 0.01, e = 0.46, Sto")
    pert_consumer_timeseries_plot(0.01, 0.46, mean, freq, seed, tsend, tvals)
    subplot(4,3,3)
    title("(C) ε = 0.01, e = 0.46, Sto")
    pert_phase_plot(0.01, 0.46, mean, freq,  seed, tsend, tvals)
    subplot(4,3,4)
    title("(D) ε = 0.01, e = 0.52, Sto")
    acf_plot(lrange, 0.01, 0.52, mean, freq,  seed, tsend, tvals)
    subplot(4,3,5)
    title("(E) ε = 0.01, e = 0.52, Sto")
    pert_consumer_timeseries_plot(0.01, 0.52, mean, freq,  seed, tsend, tvals)
    subplot(4,3,6)
    title("(F) ε = 0.01, e = 0.52, Sto")
    pert_phase_plot(0.01, 0.52, mean, freq, seed, tsend, tvals)
    subplot(4,3,7)
    title("(G) ε = 0.01, e = 0.53, Sto")
    acf_plot(lrange, 0.01, 0.53, mean, freq, seed, tsend, tvals)
    subplot(4,3,8)
    title("(H) ε = 0.01, e = 0.53, Sto")
    pert_consumer_timeseries_plot(0.01, 0.53, mean, freq, seed, tsend, tvals)
    subplot(4,3,9)
    title("(I) ε = 0.01, e = 0.53, Sto")
    pert_phase_plot(0.01, 0.53, mean, freq, seed, tsend, tvals)
    subplot(4,3,10)
    title("(J) ε = 0.01, e = 0.71, Sto")
    acf_plot(lrange, 0.01, 0.71, mean, freq, seed, tsend, tvals)
    subplot(4,3,11)
    title("(K) ε = 0.01, e = 0.71, Sto")
    pert_consumer_timeseries_plot(0.01, 0.71, mean, freq, seed, tsend, tvals)
    subplot(4,3,12)
    title("(L) ε = 0.01, e = 0.71, Sto")
    pert_phase_plot(0.01, 0.71, mean, freq, seed, tsend, tvals)
    tight_layout()
    # return acfplot_ep001
    savefig(joinpath(abpath(), "figs/ACFplot_ep001.png"))
end

# Testing whether see ACF structure in resource or closer to white noise - prediction closer to white noise because on faster time scale
# I was wrong - looks pretty similar between resources and consumers - WHY? but is time scale of ACF presently obscuring short time scale of resources?
function acf_res_plot(lrange, ep, eff, mean, freq, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, mean, freq, seed, tsend, tvals)
    acf = autocor(sol[1, :], collect(lrange))
    conf = 1.96/sqrt(length(sol))
    PyPlot.bar(collect(lrange), acf)
    hlines(0 + conf, 0, maximum(lrange))
    hlines(0 - conf, 0, maximum(lrange))
    ylim(-1,1)
    xlabel("Lag")
    return ylabel("ACF")
end

let
    mean = 0.0
    freq = 1
    seed = 3
    tsend = 10000.0
    tvals = 9900.0:1.0:tsend
    lrange = 0:1:50
    eff = 0.51
    test = figure()
    subplot(4, 1, 1)
    acf_res_plot(lrange, 0.01, eff, mean, freq, seed, tsend, tvals)
    subplot(4, 1, 2)
    acf_plot(lrange, 0.01, eff, mean, freq, seed, tsend, tvals)
    subplot(4, 1, 3)
    pert_timeseries_plot(0.01, eff, mean, freq, seed, tsend, tvals)
    subplot(4, 1, 4)
    pert_phase_plot(0.01, eff, mean, freq, seed, tsend, tvals)
    return test
end

# Testing whether high versus low frequency pert increases probablity of canard
# figure can change either efficiency and or frequency - four figures for four values of efficiency , x axis is frequency y axis is proportion

function canard_proportion(repeat, ep, eff, mean, tsend, tvals)
    freq_vals = 1.0:0.5:10.0
    canard_prop = fill(0.0, length(freq_vals))

    for (freq_i, freq_val) in enumerate(freq_vals)
        canard_count = 0
        for j in 1:repeat
            sol = RozMac_pert(ep, eff, mean, freq_val, rand(1:100000), tsend, tvals)

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
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.46, 0.0001, tsend, tvals), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.53, 0.0001, tsend, tvals), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.6, 0.0001, tsend, tvals), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.71, 0.0001, tsend, tvals), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    subplot(3,1,3)
    title("(C) Perturbation μ = 0.0002")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.46, 0.0002, tsend, tvals), label = "e = 0.46")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.53, 0.0002, tsend, tvals), label = "e = 0.53")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.6, 0.0002, tsend, tvals), label = "e = 0.6")
    scatter(1.0:0.5:10.0, canard_proportion(100, 0.01, 0.71, 0.0002, tsend, tvals), label = "e = 0.71")
    legend()
    ylim(0.0,0.22)
    xlabel("Perturbation frequency")
    ylabel("Proportion")
    tight_layout()
    # return canard_prop_plot
    savefig(joinpath(abpath(), "figs/perturbation_canard_proportion.png"))
end


## Testing whether ep = 0.01 with initial perturbation (and no futher pert) returns to equilibrium

function RozMac_initial_pert(ep, eff, mean, seed, tsend, tvals)
    Random.seed!(seed)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    par.μ = mean
    u0 = [eq_II(par)[1], eq_II(par)[2] + rand(Normal(mean, 0.01))]
    tspan = (0, tsend)
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    println(eq_II(par))
    println(u0)
    return solend = sol(tvals)
end

let
    tsend = 75.0
    sol = RozMac_initial_pert(1, 0.59, 0.0, 3, tsend, 0.0:1.0:tsend)
    test = figure()
    plot(sol.t, sol[2, :])
    return test
end

## Could you see canard (oscilattory structure) in ACF for e = 0.46 and ep = 0.01
# Answer yes when mean of pert large enough
let
    mean = 0.002
    seed = 3
    tsend = 10000.0
    tvals = 0.0:100.0:tsend
    lrange = 0:1:50
    freq = 1
    test = figure()
    subplot(3,1,1)
    acf_plot(lrange, 0.01, 0.46, mean, freq, seed, tsend, tvals)
    subplot(3,1,2)
    pert_timeseries_plot(0.01, 0.46, mean, freq, seed, tsend, tvals)
    subplot(3, 1, 3)
    pert_phase_plot(0.01, 0.46, mean, freq, seed, tsend, tvals)
    tight_layout()
    # return test
    savefig(joinpath(abpath(), "figs/perturbation_eff046ep001mean0002.png"))
end

# Can you see oscillations and autocorelation outside of initial small lags even when eff is small
# answer yes
let
    mean = 0.0
    seed = 3
    tsend = 30000.0
    tvals = 6000.0:500.0:tsend
    lrange = 0:1:25
    freq = 1
    test = figure()
    subplot(3,1,1)
    acf_plot(lrange, 0.01, 0.48, mean, freq, seed, tsend, tvals)
    subplot(3,1,2)
    pert_timeseries_plot(0.01, 0.48, mean, freq, seed, tsend, tvals)
    subplot(3, 1, 3)
    pert_phase_plot(0.01, 0.48, mean, freq, seed, tsend, tvals)
    tight_layout()
    # return test
    savefig(joinpath(abpath(), "figs/perturbation_acfoscillations_eff048.png"))
end


## Exploring whether can increase signal of white noise if match pertubation frequency to length of time takes to return to equilibrium (and ACF time steps)
let
    mean = 0.0
    seed = 1234
    tsend = 500.0
    lrange = 0:1:15
    test = figure()
    subplot(4,1,1)
    acf_plot(lrange, 1.0, 0.52, mean, 30, seed, tsend, 0.0:31.0:tsend)
    subplot(4,1,2)
    acf_plot(lrange, 1.0, 0.52, mean, 30, seed, tsend, 0.0:1.0:tsend)
    subplot(4, 1, 3)
    pert_consumer_timeseries_plot(1.0, 0.52, mean, 30, seed, tsend, 0.0:31.0:tsend)
    subplot(4, 1, 4)
    pert_consumer_timeseries_plot(1.0, 0.52, mean, 30, seed, tsend, 0.0:1.0:tsend)
    tight_layout()
    return test
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
