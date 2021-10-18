using .Threads
using CSV
using DataFrames

include("packages.jl")
include("slowfast_commoncode.jl")
include("slowfast_metabolicmodel.jl")
include("slowfast_canardfinder.jl")

function prop_canard_rednoise(model, ep, effR0, r, reps)
    count_canard = 0
    count_axial = 0
    count_nothing = 0
    tsend = 6000.0
    for i in 1:reps
        pattern = cf_returnmap(model, ep, effR0, 1, r, i, tsend, 0.0:1.0:tsend)
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

function prop_canard_rednoise_data(model, ep, effR0, reps)
    r_range = 0.0:0.01:0.9
    data = Vector{Vector{Float64}}(undef,length(r_range))
    @threads for i in eachindex(r_range)
        @inbounds data[i] =  prop_canard_rednoise(model, ep, effR0, r_range[i], reps)
    end

    return prep_data(data, r_range)
end

### RozMac model
begin
    rn_ep0079eff05_RozMac = prop_canard_rednoise_data("RozMac", 0.079, 0.5, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0079eff05_RozMac.csv"), rn_ep0079eff05_RozMac)
    rn_ep0079eff06_RozMac = prop_canard_rednoise_data("RozMac",0.079, 0.6, 1000)
    CSV.write(joinpath(abpath(),"data/rn_ep0079eff06_RozMac.csv"), rn_ep0079eff06_RozMac)
    rn_ep0079eff07_RozMac = prop_canard_rednoise_data("RozMac",0.079, 0.7, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0079eff07_RozMac.csv"), rn_ep0079eff07_RozMac)
end

begin
    rn_ep0004eff05_RozMac = prop_canard_rednoise_data("RozMac",0.004, 0.5, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0004eff05_RozMac.csv"), rn_ep0004eff05_RozMac)
    rn_ep0004eff06_RozMac = prop_canard_rednoise_data("RozMac",0.004, 0.6, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0004eff06_RozMac.csv"), rn_ep0004eff06_RozMac)
    rn_ep0004eff07_RozMac = prop_canard_rednoise_data("RozMac",0.004, 0.7, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0004eff07_RozMac.csv"), rn_ep0004eff07_RozMac)
end

### YodInn model
begin
    rn_ep00079eff12_YodInn = prop_canard_rednoise_data("YodInn", 0.0079, 1.2, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep00079eff12_YodInn.csv"), rn_ep00079eff12_YodInn)
    rn_ep00079eff10_YodInn = prop_canard_rednoise_data("YodInn", 0.0079, 1.0, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep00079eff10_YodInn.csv"), rn_ep00079eff10_YodInn)
    rn_ep00079eff08_YodInn = prop_canard_rednoise_data("YodInn", 0.079, 0.8, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep00079eff08_YodInn.csv"), rn_ep00079eff08_YodInn)
end

begin
    rn_ep0004eff12_YodInn = prop_canard_rednoise_data("YodInn", 0.004, 1.2, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0004eff12_YodInn.csv"), rn_ep0004eff12_YodInn)
    rn_ep0004eff10_YodInn = prop_canard_rednoise_data("YodInn", 0.004, 1.0, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0004eff10_YodInn.csv"), rn_ep0004eff10_YodInn)
    rn_ep0004eff08_YodInn = prop_canard_rednoise_data("YodInn", 0.004, 0.8, 1000)
    CSV.write(joinpath(abpath(), "data/rn_ep0004eff08_YodInn.csv"), rn_ep0004eff08_YodInn)
end
