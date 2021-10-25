#Increase the number of threads by going to Julia Extension settings in VS Code or setting the environment variable in bash
include("packages.jl")
include("slowfast_commoncode.jl")
include("slowfast_metabolicmodel.jl")
include("slowfast_canardfinder.jl")

function iso_plot_yodinn(resrange, par)
    data = [res_iso_yodinn(R, par) for R in resrange]
    plot(collect(resrange), data)
    return vlines(con_iso_yodinn(par), 0.0, maximum(data), colors="orange")
end


function pert_YodInn_phase_plot(ep, R0, freq, r, seed, tsend, tvals)
    sol = YodInn_pert(ep, R0, freq, r, seed, tsend, tvals)
    plot(sol[1, :], sol[2, :], color = "green")
    xlabel("Resource")
    return ylabel("Consumer")
end

let 
    ep = 0.0000001 #max ep for YodInn 0.0000001
    R0 = 0.8
    corr = 0.0
    ran = 18
    tsend = 12000.0
    println(cf_returnmap("YodInn", ep, R0, 1, corr, ran, tsend, 0.0:1.0:tsend))
    test = figure()
    iso_plot_yodinn(0.0:0.1:3.0, YodInnScalePar(R₀ = R0))
    pert_YodInn_phase_plot(ep, R0, 1, corr, ran, tsend, 0.0:1.0:tsend)
    return test
end


# let 
#         ep = 0.01
#         eff = 0.6
#         corr = 0.0
#         ran = 187
#         tsend = 6000.0
#         println(cf_returnmap("RozMac", ep, eff, 1, corr, ran, tsend, 0.0:1.0:tsend))
#         test = figure()
#         iso_plot(0.0:0.1:3.0, RozMacPar(e = eff))
#         pert_phase_plot(ep, eff, 1, corr, ran, tsend, 0.0:1.0:tsend)
#         return test
# end



function prop_canard_whitenoise(model, ep, effR0, reps, tsend)
    count_canard = 0
    count_axial = 0
    count_nothing = 0
    for i in 1:reps
        pattern = cf_returnmap(model, ep, effR0, 1, 0.0, i, tsend, 0.0:1.0:tsend)
        if pattern == "canard"
            count_canard += 1
        elseif pattern == "axial"
            count_axial += 1
        else
            count_nothing += 1
        end
    end
    return vcat(count_canard / reps, count_axial / reps, count_nothing / reps)
end

function prop_canard_whitenoise_data(model, effR0, iter, tsend)
    if model == "RozMac"
        epsilon_range = 0.001:0.001:0.15
    elseif model =="YodInn"
        epsilon_range = 0.0000001:0.00001:0.001
    else
        error("model must be either RozMac or YodInn")
    end
    data = Vector{Vector{Float64}}(undef,length(epsilon_range))
    @threads for i in eachindex(epsilon_range)
        @inbounds data[i] =  prop_canard_whitenoise(model, epsilon_range[i], effR0, iter, tsend)
    end
    return prep_data(data, epsilon_range)
end

### RozMac model
##NOTE maybe change 0.6 to 0.65 (easier to see in prop real graph)
begin
    wn_eff05_short_RozMac = prop_canard_whitenoise_data("RozMac",0.5, 1000, 6000.0)
    CSV.write(joinpath(abpath(), "data/wn_eff05_short_RozMac.csv"), wn_eff05_short_RozMac)
    wn_eff065_short_RozMac = prop_canard_whitenoise_data("RozMac", 0.65, 1000, 6000.0)
    CSV.write(joinpath(abpath(), "data/wn_eff065_short_RozMac.csv"), wn_eff065_short_RozMac)
    wn_eff07_short_RozMac = prop_canard_whitenoise_data("RozMac", 0.7, 1000, 6000.0)
    CSV.write(joinpath(abpath(), "data/wn_eff07_short_RozMac.csv"), wn_eff07_short_RozMac)
end

begin
    wn_eff05_long_RozMac = prop_canard_whitenoise_data("RozMac", 0.5, 1000, 24000.0)
    CSV.write(joinpath(abpath(), "data/wn_eff05_long_RozMac.csv"), wn_eff05_long_RozMac)
    wn_eff06_long_RozMac = prop_canard_whitenoise_data("RozMac", 0.6, 1000, 24000.0)
    CSV.write(joinpath(abpath(), "data/wn_eff06_long_RozMac.csv"), wn_eff06_long_RozMac)
    wn_eff07_long_RozMac = prop_canard_whitenoise_data("RozMac", 0.7, 1000, 24000.0)
    CSV.write(joinpath(abpath(), "data/wn_eff07_long_RozMac.csv"), wn_eff07_long_RozMac)
end

### YodInn model
begin
    wn_R12_short_YodInn = prop_canard_whitenoise_data("YodInn", 1.2, 1000, 6000.0)
    CSV.write(joinpath(abpath(), "data/wn_R12_short_YodInn.csv"), wn_R12_short_YodInn)
    wn_R10_short_YodInn = prop_canard_whitenoise_data("YodInn", 1.0, 1000, 6000.0)
    CSV.write(joinpath(abpath(), "data/wn_R10_short_YodInn.csv"), wn_R10_short_YodInn)
    wn_R08_short_YodInn = prop_canard_whitenoise_data("YodInn", 0.8, 1000, 6000.0)
    CSV.write(joinpath(abpath(), "data/wn_R08_short_YodInn.csv"), wn_R08_short_YodInn)
end

begin
    wn_R12_long_YodInn = prop_canard_whitenoise_data("YodInn", 1.2, 1000, 24000.0)
    CSV.write(joinpath(abpath(), "data/wn_R12_long_YodInn.csv"), wn_R12_long_YodInn)
    wn_R10_long_YodInn = prop_canard_whitenoise_data("YodInn", 1.0, 1000, 24000.0)
    CSV.write(joinpath(abpath(), "data/wn_R10_long_YodInn.csv"), wn_R10_long_YodInn)
    wn_R08_long_YodInn = prop_canard_whitenoise_data("YodInn", 0.8, 1000, 24000.0)
    CSV.write(joinpath(abpath(), "data/wn_R08_long_YodInn.csv"), wn_R08_long_YodInn)
end
