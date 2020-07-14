#Script for slow fast examination of time delays

include("packages.jl")
include("slowfast_commoncode.jl")
## Compare epsilon and eigenvalue
function epsilon_maxeigen_plot(eff)
    par = RozMacPar()
    par.e = eff
    epvals = 0.05:0.01:1
    max_eig = fill(0.0, length(epvals))

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        equ = eq_II(par)
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
    #savefig("figs/epsilon_eigen_plot.png")
end


function epsilon_comeigen_plot(eff, realcom)
    par = RozMacPar()
    par.e = eff
    epvals = 0.00001:0.000001:1
    eig1 = fill(0.0, length(epvals))
    eig2 = fill(0.0, length(epvals))


    if realcom == "real"
        for (epi, epval) in enumerate(epvals)
            par.ε = epval
            equ = eq_II(par)
            eig1[epi] = real.(eigvals(jacmat(roz_mac_II, equ, par))[1])
            eig2[epi] = real.(eigvals(jacmat(roz_mac_II, equ, par))[2])
        end
        PyPlot.plot(collect(epvals), eig1)
        PyPlot.plot(collect(epvals), eig2)
    else
        for (epi, epval) in enumerate(epvals)
            par.ε = epval
            equ = eq_II(par)
            eig1[epi] = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
            eig2[epi] = imag.(eigvals(jacmat(roz_mac_II, equ, par))[2])
        end
        PyPlot.plot(collect(epvals), eig1)
        PyPlot.plot(collect(epvals), eig2)
    end

    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end

let
    figure(figsize = (8,10))
    subplot(911)
    PyPlot.title("Real")
    epsilon_comeigen_plot(0.45, "real")
    subplot(912)
    epsilon_comeigen_plot(0.55, "real")
    subplot(913)
    epsilon_comeigen_plot(0.60, "real")
    subplot(914)
    epsilon_comeigen_plot(0.65, "real")
    subplot(915)
    epsilon_comeigen_plot(0.70, "real")
    subplot(916)
    epsilon_comeigen_plot(0.75, "real")
    subplot(917)
    epsilon_comeigen_plot(0.80, "real")
    subplot(918)
    epsilon_comeigen_plot(0.85, "real")
    subplot(919)
    epsilon_comeigen_plot(0.90, "real")
    gcf()
    #savefig(joinpath(abpath(), "figs/ep_eigen_real.png"))
end

let
    figure(figsize = (8,10))
    subplot(911)
    PyPlot.title("Imaginary")
    epsilon_comeigen_plot(0.45, "com")
    subplot(912)
    epsilon_comeigen_plot(0.55, "com")
    subplot(913)
    epsilon_comeigen_plot(0.60, "com")
    subplot(914)
    epsilon_comeigen_plot(0.65, "com")
    subplot(915)
    epsilon_comeigen_plot(0.70, "com")
    subplot(916)
    epsilon_comeigen_plot(0.75, "com")
    subplot(917)
    epsilon_comeigen_plot(0.80, "com")
    subplot(918)
    epsilon_comeigen_plot(0.85, "com")
    subplot(919)
    epsilon_comeigen_plot(0.90, "com")
    #gcf()
    savefig(joinpath(abpath(), "figs/ep_eigen_imag.png"))
end

function eigen_epeff(effi, ep, realcom, eigen)
    par = RozMacPar()
    par.e = effi
    par.ε = ep
    equ = eq_II(par)
    if realcom == "real"
        eig1 = real.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        eig2 = real.(eigvals(jacmat(roz_mac_II, equ, par))[2])
    else
        eig1 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        eig2 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[2])
    end
    if eigen == 1
        return eig1
    else
        return eig2
    end
end

effrange = 0.45:0.01:0.9
eprange = 0.001:0.001:1

function eigen_contour(realcom, eigen)
    return [eigen_epeff(ei, epi, realcom,eigen) for ei in effrange, epi in eprange]
end

maxrint = [-0.36, -0.30, -0.24, -0.18, -0.12, -0.06, -0.0001, 0.0001, 0.06, 0.12, 0.18]
subrint = [-2.00, -1.75, -1.5, -1.25, -1.00, -0.75, -0.5, -0.25, -0.0001, 0.0001, 0.25]
negcint = sort([0.00, -0.0001, -0.08, -0.16, -0.24, -0.32, -0.40, -0.48, -0.56, -0.64])
poscint = [0.00, 0.0001, 0.08, 0.16, 0.24, 0.32, 0.40, 0.48, 0.56, 0.64]

function eigen_contour_plot(realcom, eigen, levs)
    contourf(eprange, effrange, eigen_contour(realcom, eigen), levels = levs)
    colorbar()
end

let
    figure()
    subplot(211)
    PyPlot.title("Real")
    eigen_contour_plot("real", 1, subrint)
    ylabel("Efficiency")
    subplot(212)
    eigen_contour_plot("real", 2, maxrint)
    xlabel("ε")
    ylabel("Efficiency")
    #gcf()
    savefig(joinpath(abpath(), "figs/eigencontour_real.png"))
end

let
    figure()
    subplot(211)
    PyPlot.title("Imaginary")
    eigen_contour_plot("imag", 1, negcint)
    ylabel("Efficiency")
    subplot(212)
    eigen_contour_plot("imag", 2, poscint)
    ylabel("Efficiency")
    xlabel("ε")
    #gcf()
    savefig(joinpath(abpath(), "figs/eigencontour_imag.png"))
end


#When efficiency is 0.55, find epsilon value where imaginary numbers start


par = RozMacPar()
par.e = 0.55
equ = eq_II(par)
epvals = 0.00001:0.000001:1

for (epi, epval) in enumerate(epvals)
    par.ε = epval
    eig1 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
    #println(eig1)
    if eig1 < 0 || eig1 > 0
         return epval
         break
         #println(epval)
    end
end

# When epp is 1, find the excitable and non excitable phases for e - deterministic
function findRCdivide(ep)
    par = RozMacPar()
    par.ε = ep
    evals = 0.441:0.00005:0.9

    for (ei, eval) in enumerate(evals)
        par.e = eval
        equ = eq_II(par)
        eig1 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        if eig1 < 0 || eig1 > 0
            return eval
            break
        end
    end
end

findRCdivide(1.0)

# # create graph of epsilon on x and efficiency value where RC divide by transposing graph of efficiency on x and epsilon value that creates RC divide
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

findRCdivide_effx(0.5225)

function findRCdivide_effx_plot()
    evals = 0.5225:0.0005:1
    epRC = fill(0.0, length(evals))
    for (ei, eval) in enumerate(evals)
        epRC[ei] = findRCdivide_effx(eval)
    end
    return hcat(reverse(epRC), reverse(collect(evals)))
end

let
    figure()
    data = findRCdivide_effx_plot()
    plot(data[:,1], data[:,2])
    xlim(-0.05,1.05)
    ylim(0.410, 1)
    ylabel("Efficiency where R/C divide occurs", fontsize = 15)
    xlabel("ε", fontsize = 15)
    hlines([0.441,0.5225, 0.710], xmin = -0.05, xmax = 1.05, linestyles = "dashed")
    annotate("TC", (408, 54), xycoords = "figure points", fontsize = 12)
    annotate("R/C", (408, 89), xycoords = "figure points", fontsize = 12)
    annotate("H", (408, 175), xycoords = "figure points", fontsize = 12)
    # gcf()
    savefig(joinpath(abpath(), "figs/epsilonxaxis_RCdivide.png"))
end


# graph real and imag for deterministis ep = 1.0 for my own understanding
function eff_maxeigen_plot()
    par = RozMacPar()
    par.ε = 1.0
    evals = 0.441:0.005:0.9
    max_eig = fill(0.0, length(evals))

    for (ei, eval) in enumerate(evals)
        par.e = eval
        equ = eq_II(par)
        max_eig[ei] = λ_stability(jacmat(roz_mac_II, equ, par))
    end

    return PyPlot.plot(collect(evals), max_eig)
    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end

let
    figure()
    eff_maxeigen_plot()
    gcf()
end

function eff_comeigen_plot(realcom)
    par = RozMacPar()
    par.ε = 1.0
    evals = 0.441:0.005:0.9
    eig1 = fill(0.0, length(evals))
    eig2 = fill(0.0, length(evals))


    if realcom == "real"
        for (ei, eval) in enumerate(evals)
            par.e = eval
            equ = eq_II(par)
            eig1[ei] = real.(eigvals(jacmat(roz_mac_II, equ, par))[1])
            eig2[ei] = real.(eigvals(jacmat(roz_mac_II, equ, par))[2])
        end
        PyPlot.plot(collect(evals), eig1)
        PyPlot.plot(collect(evals), eig2)
    else
        for (ei, eval) in enumerate(evals)
            par.e = eval
            equ = eq_II(par)
            eig1[ei] = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
            eig2[ei] = imag.(eigvals(jacmat(roz_mac_II, equ, par))[2])
        end
        PyPlot.plot(collect(evals), eig1)
        PyPlot.plot(collect(evals), eig2)
    end

    #ylabel("Implicit Lag", fontsize = 15)
    # ylim(-0.01, 0.51)
    #xlabel("ε", fontsize = 15)
end

let
    figure()
    subplot(211)
    eff_comeigen_plot("real")
    subplot(212)
    eff_comeigen_plot("complex")
    gcf()
end
