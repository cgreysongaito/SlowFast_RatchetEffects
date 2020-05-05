#Script for slow fast examination of time delays

include("packages.jl")
include("slowfast_commoncode.jl")

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
