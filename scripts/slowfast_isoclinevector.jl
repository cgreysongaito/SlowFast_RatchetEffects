#Script for slow fast examination of time delays

include("packages.jl")
include("slowfast_commoncode.jl")
# Plot isoclines and vector fields

function roz_mac_res(R, C, p)
    @unpack r, k, h, a, m = p
    return r * R * (1 - R / k) - (a * R * C / (1 + a * h * R) )
end

function roz_mac_con(R, C, eff, ep, p)
    @unpack h, a, m = p
    return ep * ( ( eff * a * R * C ) / (1 + a * h * R) - m * C )
end

function con_iso(eff, p)
    @unpack m, a, h = p
    m / (a * (eff - h * m))
end

function res_iso(R, p)
    @unpack a, k, r, h = p
    r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
end

function roz_mac_plot(eff, ep)
    minval = 0
    maxval = 3
    steps = 100
    resconrange = range(minval,stop=maxval,length=steps)

    U = [roz_mac_res(R, C, par_rozmac) for C in resconrange, R in resconrange]
    V = [roz_mac_con(R, C, eff, ep, par_rozmac) for C in resconrange, R in resconrange]
    speed = sqrt.(U.^2 .+ V.^2)
    lw = 5 .* speed ./ maximum(speed) # Line Widths
    streamplot(collect(resconrange), collect(resconrange), U, V, density = 0.6, color = "k", linewidth = lw)
    PyPlot.plot(collect(resconrange), [res_iso(R, par_rozmac) for R in resconrange])
    return PyPlot.plot(repeat([con_iso(eff, par_rozmac)],100),collect(resconrange))
end

# - Streamlot figure
let
    figure(figsize = (8,10))
    subplot(421)
    roz_mac_plot(0.45, 1)
    subplot(422)
    roz_mac_plot(0.45, 0.01)
    subplot(423)
    roz_mac_plot(0.6, 1)
    subplot(424)
    roz_mac_plot(0.6, 0.01)
    subplot(425)
    roz_mac_plot(0.71, 1) #need to differentiate and find max
    subplot(426)
    roz_mac_plot(0.71, 0.01)
    subplot(427)
    roz_mac_plot(1.1, 1)
    subplot(428)
    roz_mac_plot(1.1, 0.01)
    savefig("figs/vectorfieldplot.png")
    #gcf()
end

 
