#Script for slow fast examination of time delays

include("packages.jl")
include("slowfast_commoncode.jl")
## Compare epsilon and eigenvalue
function epsilon_maxeigen_plot(eff)
    par = RozMacPar()
    par.e = eff
    equ = eq_II(par.e, par)
    epvals = 0.05:0.01:1
    max_eig = fill(0.0, length(epvals))

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        max_eig[epi] = λ_stability(jacmat(roz_mac_II, equ, par))
    end

    return PyPlot.plot(collect(epvals), max_eig)
    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end


let
    figure(figsize = (8,10))
    subplot(411)
    epsilon_maxeigen_plot(0.45)
    ylabel("Dominant λ")
    subplot(412)
    epsilon_maxeigen_plot(0.6)
    ylabel("Dominant λ")
    subplot(413)
    epsilon_maxeigen_plot(0.74)
    ylabel("Dominant λ")
    subplot(414)
    epsilon_maxeigen_plot(0.9)
    ylabel("Dominant λ")
    xlabel("ε")
    gcf()
    savefig("figs/epsilon_eigen_plot.png")
end


function epsilon_comeigen_plot(eff)
    par = RozMacPar()
    par.e = eff
    equ = eq_II(par.e, par)
    epvals = 0.0001:0.0001:1
    eigr1 = fill(0.0, length(epvals))
    eigi1 = fill(0.0, length(epvals))
    eigr2 = fill(0.0, length(epvals))
    eigi2 = fill(0.0, length(epvals))

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        eigr1[epi] = real.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        eigi1[epi] = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        eigr2[epi] = real.(eigvals(jacmat(roz_mac_II, equ, par))[2])
        eigi2[epi] = imag.(eigvals(jacmat(roz_mac_II, equ, par))[2])
    end

    PyPlot.plot(eigr1, eigi1)
    PyPlot.plot(eigr2, eigi2)
    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end

let
    figure()
    epsilon_comeigen_plot(0.78)
    gcf()
end
