include("packages.jl")
include("slowfast_commoncode.jl")

#Figure 1
include("slowfast_eigen.jl")

let
    prop_real = figure()
    data = findRCdivide_epx_data()
    plot(data[:,1], data[:,2], color = "black")
    fill_between(data[:,1], data[:,2], color = "blue", alpha=0.3)
    fill_between(data[:,1], fill(1.0, length(data[:,2])), data[:,2], color = "orange", alpha=0.3)
    xlim(-0.05,2.05)
    #ylim(0.410, 1)
    ylabel("Proportion Real", fontsize = 15)
    xlabel("ε", fontsize = 15)
    return prop_real
    # savefig(joinpath(abpath(), "figs/epsilonxaxis_propReal.png"))
end


let
    data = eff_maxeigen_data(0.01)
    maxeigen_ep001 = figure()
    plot(data[:,1], data[:,2], color = "black")
    title("ϵ = 0.01")
    ylabel("Re(λₘₐₓ)", fontsize = 15)
    xlim(0.4, 0.8)
    ylim(-0.35, 0.1)
    xlabel("Efficiency", fontsize = 15)
    hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
    vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
    fill([0.441,0.65415, 0.65415, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
    fill([0.65415, 0.710, 0.710, 0.65415], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
    annotate("TC", (90, 298), xycoords = "figure points", fontsize = 12)
    annotate("H", (335, 298), xycoords = "figure points", fontsize = 12)
    # return maxeigen_ep001
    savefig(joinpath(abpath(), "figs/maxeigen_ep001.png"))
end

let
    data = eff_maxeigen_data(0.8)
    maxeigen_ep08 = figure()
    plot(data[:,1], data[:,2], color = "black")
    title("ϵ = 0.8")
    ylabel("Re(λₘₐₓ)", fontsize = 15)
    xlabel("Efficiency", fontsize = 15)
    xlim(0.4, 0.8)
    ylim(-0.35, 0.1)
    hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
    vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
    fill([0.441,0.529, 0.529, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
    fill([0.529,0.710, 0.710, 0.529], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
    annotate("TC", (86, 298), xycoords = "figure points", fontsize = 12)
    annotate("H", (330, 298), xycoords = "figure points", fontsize = 12)
    return maxeigen_ep08
    # savefig(joinpath(abpath(), "figs/maxeigen_ep08.png"))
end

# Figure 2 (white noise and isoclines)
using DataFrames
using CSV

data05 = CSV.read("/home/chrisgg/julia/TimeDelays/canard05eq.csv", DataFrame)
data06 = CSV.read("/home/chrisgg/julia/TimeDelays/canard06eq.csv", DataFrame)
data07 = CSV.read("/home/chrisgg/julia/TimeDelays/canard07eq.csv", DataFrame)
data08 = CSV.read("/home/chrisgg/julia/TimeDelays/canard08eq.csv", DataFrame)

using Plots
pyplot()

short_data = CSV.read("/home/chrisgg/julia/TimeDelays/canard_whitenoise_short.csv", DataFrame)
long_data = CSV.read("/home/chrisgg/julia/TimeDelays/canard_whitenoise_long.csv", DataFrame)
longer_data = CSV.read("/home/chrisgg/julia/TimeDelays/canard_whitenoise_longer.csv", DataFrame)

begin
    plot(collect(0.001:0.001:0.12), data05.prop[1:120], color = "black")
    ylims!(0,1)
    lens!([0,0.008], [0,0.003], inset = (1, bbox(0.5, 0.0, 0.4, 0.4)))
end

let
    figure2 = figure(figsize = (10,5))
    subplot(2,4,1)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.5))
    title("Non-excitable\n(Real λ)")
    ylim(0,2.5)
    xlim(0.0,3.0)
    ylabel("Consumer")
    xlabel("Resource")
    subplot(2,4,2)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.6))
    title("Excitable\n(Complex λ)")
    xlabel("Resource")
    ylim(0,2.5)
    xlim(0.0,3.0)
    subplot(2,4,3)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.7))
    title("Excitable\n(Near Hopf)")
    xlabel("Resource")
    ylim(0,2.5)
    xlim(0.0,3.0)
    subplot(2,4,4)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.8))
    ylim(0,2.5)
    xlim(0.0,3.0)
    title("Canard")
    xlabel("Resource")
    subplot(2,4,5)
    plot(collect(0.001:0.001:0.15), short_data.prop05, color = "black", linestyle = "dashed")
    # plot(collect(0.001:0.001:0.15), long_data.prop05, color = "black", linestyle = "dashed")
    plot(collect(0.001:0.001:0.15), longer_data.prop05, color = "black", linestyle = "dotted")
    # plot(collect(0.2:0.1:1.0), data05.prop[121:end], color = "black")
    ylabel("Proportion with\nquasi-canard")
    xlabel("ϵ")
    ylim(0,1)
    subplot(2,4,6)
    plot(collect(0.001:0.001:0.15), short_data.prop06, color = "black", linestyle = "dashed")
    # plot(collect(0.001:0.001:0.15), long_data.prop06, color = "black", linestyle = "dashed")
    plot(collect(0.001:0.001:0.15), longer_data.prop06, color = "black", linestyle = "dotted")
    # plot(collect(0.2:0.1:1.0), data06.prop[121:end], color = "black")
    ylim(0,1)
    xlabel("ϵ")
    subplot(2,4,7)
    plot(collect(0.001:0.001:0.15), short_data.prop07, color = "black", linestyle = "dashed")
    # plot(collect(0.001:0.001:0.15), long_data.prop07, color = "black", linestyle = "dashed")
    plot(collect(0.001:0.001:0.15), longer_data.prop07, color = "black", linestyle = "dotted")
    xlabel("ϵ")
    # plot(collect(0.2:0.1:1.0), data07.prop[121:end], color = "black")
    ylim(0,1)
    subplot(2,4,8)
    plot(collect(0.001:0.001:0.15), short_data.prop08, color = "black", linestyle = "dashed")
    # plot(collect(0.001:0.001:0.15), long_data.prop08, color = "black", linestyle = "dashed")
    plot(collect(0.001:0.001:0.15), longer_data.prop08, color = "black", linestyle = "dotted")
    xlabel("ϵ")
    # plot(collect(0.2:0.1:1.0), data08.prop[121:end], color = "black")
    ylim(0,1)
    tight_layout()
    return figure2
    # savefig(joinpath(abpath(), "figs/canard_whitenoise_prop.png"))
end


# Figure 3 (time series of canards)
let
    sol_05_ep001 = RozMac_pert(0.01, 0.5, 1, 0.0, 1234, 5000.0, 3500.0:2.0:5000.0)
    sol_06_ep001 = RozMac_pert(0.01, 0.6, 1, 0.0, 1234, 5000.0, 3500.0:2.0:5000.0)
    sol_07_ep001 = RozMac_pert(0.01, 0.7, 1, 0.0, 1234, 5000.0, 3500.0:2.0:5000.0)
    timeseries = figure(figsize=(10,4))
    subplot(2,3,1)
    plot(sol_05_ep001.t, sol_05_ep001.u)
    title("e = 0.5, ϵ = 0.01")
    ylim(0,2.5)
    ylabel("Biomass")
    subplot(2,3,2)
    plot(sol_06_ep001.t, sol_06_ep001.u)
    title("e = 0.6, ϵ = 0.01")
    ylim(0,2.5)
    xlabel("Time")
    subplot(2,3,3)
    plot(sol_07_ep001.t, sol_07_ep001.u)
    title("e = 0.7, ϵ = 0.01")
    ylim(0,2.5)
    subplot(2,3,4)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.5))
    plot(sol_05_ep001[1, :], sol_05_ep001[2, :], color = "black")
    xlabel("Resource")
    ylabel("Consumer")
    ylim(0,2.5)
    xlim(0,3)
    subplot(2,3,5)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.6))
    plot(sol_06_ep001[1, :], sol_06_ep001[2, :], color = "black")
    xlabel("Resource")
    ylabel("Consumer")
    ylim(0,2.5)
    xlim(0,3)
    tight_layout()
    return timeseries
    # savefig(joinpath(abpath(), "figs/figure2b.png"))
end

# Figure 4

let
    figure2 = figure(figsize = (6,3))
    subplot(1,3,1)
    plot(collect(0.0:0.1:0.9), rn_fulldataset[1], color="black")
    fill_between(collect(0.0:0.1:0.9), fill(1.0,10), rn_fulldataset[1], color="#E66100", alpha = 0.3)
    fill_between(collect(0.0:0.1:0.9), rn_fulldataset[1], color="#5D3A9B", alpha = 0.3)
    ylabel("Proportion with quasi-canards")
    ylim(0,1)
    title("e = 0.5, ϵ = 0.079")
    subplot(1,3,2)
    plot(collect(0.0:0.1:0.9), rn_fulldataset[2], color="black")
    fill_between(collect(0.0:0.1:0.9), fill(1.0,10), rn_fulldataset[2], color="#E66100", alpha = 0.3)
    fill_between(collect(0.0:0.1:0.9), rn_fulldataset[2], color="#5D3A9B", alpha = 0.3)
    xlabel("Noise correlation (t = -1)")
    ylim(0,1)
    title("e = 0.6, ϵ = 0.079")
    subplot(1,3,3)
    plot(collect(0.0:0.1:0.9), rn_fulldataset[3], color="black")
    fill_between(collect(0.0:0.1:0.9), fill(1.0,10), rn_fulldataset[3], color="#E66100", alpha = 0.3)
    fill_between(collect(0.0:0.1:0.9), rn_fulldataset[3], color="#5D3A9B", alpha = 0.3)
    ylim(0,1)
    title("e = 0.7, ϵ = 0.079")
    tight_layout()
    # return figure2
    savefig(joinpath(abpath(), "figs/canard_rednoise_prop.png"))
end

let
    data = eff_maxeigen_data(0.079)
    maxeigen_ep0079 = figure()
    plot(data[:,1], data[:,2], color = "black")
    ylabel("Re(λₘₐₓ)", fontsize = 15)
    xlim(0.4, 0.8)
    ylim(-0.35, 0.1)
    xlabel("Efficiency", fontsize = 15)
    hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
    vlines([0.5,0.6,0.7], ymin = -0.35, ymax = 0.1, color = "green")
    vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
    fill([0.441,0.60055, 0.60055, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
    fill([0.60055,0.71, 0.71, 0.60055], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
    annotate("TC", (90, 298), xycoords = "figure points", fontsize = 12)
    annotate("H", (335, 298), xycoords = "figure points", fontsize = 12)
    return maxeigen_ep0079
    #savefig(joinpath(abpath(), "figs/maxeigen_ep0079.png"))
end




let
    data = eff_maxeigen_data(0.079)
    maxeigen_ep0079 = figure()
    plot(data[:,1], data[:,2], color = "black")
    ylabel("Re(λₘₐₓ)", fontsize = 15)
    xlim(0.4, 0.8)
    ylim(-0.35, 0.1)
    xlabel("Efficiency", fontsize = 15)
    hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
    vlines([0.5,0.6,0.7], ymin = -0.35, ymax = 0.1, color = "green")
    vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
    fill([0.441,0.60055, 0.60055, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
    fill([0.60055,0.71, 0.71, 0.60055], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
    annotate("TC", (90, 298), xycoords = "figure points", fontsize = 12)
    annotate("H", (335, 298), xycoords = "figure points", fontsize = 12)
    # return maxeigen_ep0079
    savefig(joinpath(abpath(), "figs/maxeigen_ep0079.png"))
end
