include("packages.jl")
include("slowfast_commoncode.jl")

# Figure 2 (primer of different trajectories)

let
    sol_axial = RozMac_pert(0.01, 0.48, 1, 0.8, 9, 1500.0, 0.0:2.0:1500.0)
    sol_qc = RozMac_pert(0.01, 0.48, 1, 0.8, 123, 2500.0, 1750.0:2.0:2500.0)
    sol_equil = RozMac_pert(0.01, 0.48, 1, 0.0, 1234, 5000.0, 0.0:2.0:5000.0)
    figure2 = figure(figsize=(3,10))
    subplot(5,1,1)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.48))  
    plot(sol_qc[1, :], sol_qc[2, :], color = "#440154FF")
    plot(sol_axial[1, :], sol_axial[2, :], color = "#73D055FF")
    xlabel("Resource")
    ylabel("Consumer")
    ylim(0,2.5)
    xlim(0,3)
    subplot(5,1,2)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.48))
    plot(sol_equil[1, :], sol_equil[2, :], color = "#FDE725FF")
    xlabel("Resource")
    ylabel("Consumer")
    ylim(0,2.5)
    xlim(0,3)
    subplot(5,1,3)
    plot(sol_qc.t, sol_qc.u)
    ylim(0,3.0)
    xlabel("Time")
    ylabel("Biomass")
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(5,1,4)
    plot(sol_axial.t, sol_axial.u)
    ylim(0,3.0)
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    xlabel("Time")
    ylabel("Biomass")
    subplot(5,1,5)
    plot(sol_equil.t, sol_equil.u)
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    xlabel("Time")
    ylabel("Biomass")
    ylim(0,3.0)
    tight_layout()
    # return figure2
    savefig(joinpath(abpath(), "figs/phase_timeseries_examples.pdf"))
end


# Figure 3 (white noise and isoclines)
using DataFrames
using CSV

wn_data05_short = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff05_short.csv", DataFrame)
wn_data07_short = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff07_short.csv", DataFrame)
wn_data06_short = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff06_short.csv", DataFrame)
wn_data08_short = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff08_short.csv", DataFrame)
wn_data05_long = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff05_long.csv", DataFrame)
wn_data06_long = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff06_long.csv", DataFrame)
wn_data07_long = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff07_long.csv", DataFrame)
wn_data08_long = CSV.read("/home/chrisgg/julia/TimeDelays/data/wn_eff08_long.csv", DataFrame)

let
    figure3 = figure(figsize = (8,5))
    subplot(2,3,1)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.5))
    title("Non-excitable\n(Real λ)")
    ylim(0,2.5)
    xlim(0.0,3.0)
    ylabel("Consumer")
    xlabel("Resource")
    subplot(2,3,2)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.6))
    title("Excitable\n(Complex λ)")
    xlabel("Resource")
    ylim(0,2.5)
    xlim(0.0,3.0)
    subplot(2,3,3)
    iso_plot(0.0:0.1:3.0, RozMacPar(e = 0.7))
    title("Excitable\n(Near Hopf)")
    xlabel("Resource")
    ylim(0,2.5)
    xlim(0.0,3.0)
    subplot(2,3,4)
    plot(1 ./ wn_data05_short.xrange, wn_data05_short.canard, color = "black", linestyle = "dashed")
    plot(1 ./ wn_data05_long.xrange, wn_data05_long.canard, color = "black", linestyle = "dotted")
    ylabel("Proportion with\nquasi-canard")
    xlabel("1/ϵ")
    ylim(0,1)
    subplot(2,3,5)
    plot(1 ./ wn_data06_short.xrange, wn_data06_short.canard, color = "black", linestyle = "dashed")
    plot(1 ./ wn_data06_long.xrange, wn_data06_long.canard, color = "black", linestyle = "dotted")
    ylim(0,1)
    xlabel("1/ϵ")
    subplot(2,3,6)
    plot(1 ./ wn_data07_short.xrange, wn_data07_short.canard, color = "black", linestyle = "dashed")
    plot(1 ./ wn_data07_long.xrange, wn_data07_long.canard, color = "black", linestyle = "dotted")
    xlabel("1/ϵ")
    ylim(0,1)
    tight_layout()
    annotate("a)", (32, 330), xycoords = "figure points", fontsize = 15)
    annotate("b)", (200, 330), xycoords = "figure points", fontsize = 15)
    annotate("c)", (370, 330), xycoords = "figure points", fontsize = 15)
    annotate("d)", (32, 165), xycoords = "figure points", fontsize = 15)
    annotate("e)", (200, 165), xycoords = "figure points", fontsize = 15)
    annotate("f)", (370, 165), xycoords = "figure points", fontsize = 15)

    return figure3
    # savefig(joinpath(abpath(), "figs/canard_whitenoise_prop.pdf"))
end


# Figure 4

rn_data_ep0004eff05 = CSV.read("/home/chrisgg/julia/TimeDelays/data/rn_ep0004eff05.csv", DataFrame)
rn_data_ep0004eff06 = CSV.read("/home/chrisgg/julia/TimeDelays/data/rn_ep0004eff06.csv", DataFrame)
rn_data_ep0004eff07 = CSV.read("/home/chrisgg/julia/TimeDelays/data/rn_ep0004eff07.csv", DataFrame)
rn_data_ep0079eff05 = CSV.read("/home/chrisgg/julia/TimeDelays/data/rn_ep0079eff05.csv", DataFrame)
rn_data_ep0079eff06 = CSV.read("/home/chrisgg/julia/TimeDelays/data/rn_ep0079eff06.csv", DataFrame)
rn_data_ep0079eff07 = CSV.read("/home/chrisgg/julia/TimeDelays/data/rn_ep0079eff07.csv", DataFrame)


let
    figure4 = figure(figsize = (8,4))
   

    subplot(2,3,1)
    fill_between(rn_data_ep0079eff05.xrange, rn_data_ep0079eff05.canard, color="#440154FF")
    fill_between(rn_data_ep0079eff05.xrange, rn_data_ep0079eff05.canard, rn_data_ep0079eff05.canard_plus_axial, color="#73D055FF")
    fill_between(rn_data_ep0079eff05.xrange, fill(1.0,91), rn_data_ep0079eff05.canard_plus_axial, color="#FDE725FF")
    ylabel("Proportion")
    ylim(0,1)
    subplot(2,3,2)
    fill_between(rn_data_ep0079eff06.xrange, rn_data_ep0079eff06.canard, color="#440154FF")
    fill_between(rn_data_ep0079eff06.xrange, rn_data_ep0079eff06.canard, rn_data_ep0079eff06.canard_plus_axial, color="#73D055FF")
    fill_between(rn_data_ep0079eff06.xrange, fill(1.0,91), rn_data_ep0079eff06.canard_plus_axial, color="#FDE725FF")
    xlabel("Noise correlation (t = -1)")
    ylim(0,1)
    subplot(2,3,3)
    fill_between(rn_data_ep0079eff07.xrange, rn_data_ep0079eff07.canard, color="#440154FF")
    fill_between(rn_data_ep0079eff07.xrange, rn_data_ep0079eff07.canard, rn_data_ep0079eff07.canard_plus_axial, color="#73D055FF")
    fill_between(rn_data_ep0079eff07.xrange, fill(1.0,91), rn_data_ep0079eff07.canard_plus_axial, color="#FDE725FF")
    ylim(0,1)

    subplot(2,3,4)
    fill_between(rn_data_ep0004eff05.xrange, rn_data_ep0004eff05.canard, color="#440154FF")
    fill_between(rn_data_ep0004eff05.xrange, rn_data_ep0004eff05.canard, rn_data_ep0004eff05.canard_plus_axial, color="#73D055FF")
    fill_between(rn_data_ep0004eff05.xrange, fill(1.0,91), rn_data_ep0004eff05.canard_plus_axial, color="#FDE725FF")
    ylabel("Proportion")
    ylim(0,1)
    subplot(2,3,5)
    fill_between(rn_data_ep0004eff06.xrange, rn_data_ep0004eff06.canard, color="#440154FF")
    fill_between(rn_data_ep0004eff06.xrange, rn_data_ep0004eff06.canard, rn_data_ep0004eff06.canard_plus_axial, color="#73D055FF")
    fill_between(rn_data_ep0004eff06.xrange, fill(1.0,91), rn_data_ep0004eff06.canard_plus_axial, color="#FDE725FF")
    xlabel("Noise correlation (t = -1)")
    ylim(0,1)
    subplot(2,3,6)
    fill_between(rn_data_ep0004eff07.xrange, rn_data_ep0004eff07.canard, color="#440154FF")
    fill_between(rn_data_ep0004eff07.xrange, rn_data_ep0004eff07.canard, rn_data_ep0004eff07.canard_plus_axial, color="#73D055FF")
    fill_between(rn_data_ep0004eff07.xrange, fill(1.0,91), rn_data_ep0004eff07.canard_plus_axial, color="#FDE725FF")
    ylim(0,1)
    tight_layout()
    annotate("a)", (20, 280), xycoords = "figure points", fontsize = 15)
    annotate("b)", (205, 280), xycoords = "figure points", fontsize = 15)
    annotate("c)", (385, 280), xycoords = "figure points", fontsize = 15)
    annotate("d)", (20, 140), xycoords = "figure points", fontsize = 15)
    annotate("e)", (205, 140), xycoords = "figure points", fontsize = 15)
    annotate("f)", (385, 140), xycoords = "figure points", fontsize = 15)
    
    return figure4
    # savefig(joinpath(abpath(), "figs/canard_rednoise_prop.pdf"))
end



# Figure 5 (schematic B)

let 
    figure5 = figure(figsize = (3,4))
    subplot(2,1,1)
    roz_mac_plot(0.01, 0.48, 2, 0.3)
    subplot(2,1,2)
    roz_mac_plot(0.01, 0.48, 2, 0.3)
    tight_layout()
    # return figure5
    savefig(joinpath(abpath(), "figs/figure5.pdf"))
end


let 
    
end

#Original Figure 1 (proportion real) - replaced with proof
# include("slowfast_eigen.jl")

# let
#     prop_real = figure()
#     data = findRCdivide_epx_data()
#     plot(data[:,1], data[:,2], color = "black")
#     fill_between(data[:,1], data[:,2], color = "blue", alpha=0.3)
#     fill_between(data[:,1], fill(1.0, length(data[:,2])), data[:,2], color = "orange", alpha=0.3)
#     xlim(-0.05,2.05)
#     #ylim(0.410, 1)
#     ylabel("Proportion Real", fontsize = 15)
#     xlabel("ε", fontsize = 15)
#     return prop_real
#     # savefig(joinpath(abpath(), "figs/epsilonxaxis_propReal.png"))
# end


# let
#     data = eff_maxeigen_data(0.01)
#     maxeigen_ep001 = figure()
#     plot(data[:,1], data[:,2], color = "black")
#     title("ϵ = 0.01")
#     ylabel("Re(λₘₐₓ)", fontsize = 15)
#     xlim(0.4, 0.8)
#     ylim(-0.35, 0.1)
#     xlabel("Efficiency", fontsize = 15)
#     hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
#     vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
#     fill([0.441,0.65415, 0.65415, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
#     fill([0.65415, 0.710, 0.710, 0.65415], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
#     annotate("TC", (90, 298), xycoords = "figure points", fontsize = 12)
#     annotate("H", (335, 298), xycoords = "figure points", fontsize = 12)
#     # return maxeigen_ep001
#     savefig(joinpath(abpath(), "figs/maxeigen_ep001.png"))
# end

# let
#     data = eff_maxeigen_data(0.8)
#     maxeigen_ep08 = figure()
#     plot(data[:,1], data[:,2], color = "black")
#     title("ϵ = 0.8")
#     ylabel("Re(λₘₐₓ)", fontsize = 15)
#     xlabel("Efficiency", fontsize = 15)
#     xlim(0.4, 0.8)
#     ylim(-0.35, 0.1)
#     hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
#     vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
#     fill([0.441,0.529, 0.529, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
#     fill([0.529,0.710, 0.710, 0.529], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
#     annotate("TC", (86, 298), xycoords = "figure points", fontsize = 12)
#     annotate("H", (330, 298), xycoords = "figure points", fontsize = 12)
#     return maxeigen_ep08
#     # savefig(joinpath(abpath(), "figs/maxeigen_ep08.png"))
# end


# let
#     data = eff_maxeigen_data(0.079)
#     maxeigen_ep0079 = figure()
#     plot(data[:,1], data[:,2], color = "black")
#     ylabel("Re(λₘₐₓ)", fontsize = 15)
#     xlim(0.4, 0.8)
#     ylim(-0.35, 0.1)
#     xlabel("Efficiency", fontsize = 15)
#     hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
#     vlines([0.5,0.6,0.7], ymin = -0.35, ymax = 0.1, color = "green")
#     vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
#     fill([0.441,0.60055, 0.60055, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
#     fill([0.60055,0.71, 0.71, 0.60055], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
#     annotate("TC", (90, 298), xycoords = "figure points", fontsize = 12)
#     annotate("H", (335, 298), xycoords = "figure points", fontsize = 12)
#     # return maxeigen_ep0079
#     savefig(joinpath(abpath(), "figs/maxeigen_ep0079.png"))
# end


# let
#     data = eff_maxeigen_data(0.004)
#     maxeigen_ep0079 = figure()
#     plot(data[:,1], data[:,2], color = "black")
#     ylabel("Re(λₘₐₓ)", fontsize = 15)
#     xlim(0.4, 0.8)
#     ylim(-0.35, 0.1)
#     xlabel("Efficiency", fontsize = 15)
#     hlines(0.0, 0.4,0.8, linestyles = "dashed", linewidth = 0.5)
#     vlines([0.5,0.6,0.7], ymin = -0.35, ymax = 0.1, color = "green")
#     vlines([0.441, 0.710], ymin = -0.35, ymax = 0.1, linestyles = "dashed")
#     fill([0.441,0.60055, 0.60055, 0.441], [0.1, 0.1, -0.35, -0.35], "blue", alpha=0.3)
#     fill([0.60055,0.71, 0.71, 0.60055], [0.1, 0.1, -0.35, -0.35], "orange", alpha=0.3)
#     annotate("TC", (90, 298), xycoords = "figure points", fontsize = 12)
#     annotate("H", (335, 298), xycoords = "figure points", fontsize = 12)
#     return maxeigen_ep0079
#     #savefig(joinpath(abpath(), "figs/maxeigen_ep0079.png"))
# end