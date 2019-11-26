#Script for slow fast examination of time delays

include("packages.jl")
include("slowfast_commoncode.jl")
## Compare epsilon and reactivity
function epsilon_reac_plot(eff)
    par = RozMacPar()
    par.e = eff
    equ = eq_II(par)
    epvals = 0.0005:0.0001:1.5
    reac = fill(0.0, length(epvals))

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        reac[epi] = ν_stability(jacmat(roz_mac_II, equ, par))
    end

    return PyPlot.plot(collect(epvals), reac)
    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end

let
    figure(figsize = (8,10))
    subplot(411)
    epsilon_reac_plot(0.45)
    ylabel("Reactivity")
    subplot(412)
    epsilon_reac_plot(0.6)
    ylabel("Reactivity")
    subplot(413)
    epsilon_reac_plot(0.74)
    ylabel("Reactivity")
    subplot(414)
    epsilon_reac_plot(0.9)
    ylabel("Reactivity")
    xlabel("ε")
    gcf()
    savefig("figs/epsilon_reac_plot.png")
end
