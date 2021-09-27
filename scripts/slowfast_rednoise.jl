using .Threads
using CSV
using DataFrames

include("/home/chrisgg/julia/TimeDelays/scripts/packages.jl")
include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_commoncode.jl")
include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_canardfinder.jl")

function prop_canard_rednoise(ep, eff, r, reps)
    count_canard = 0
    count_axial = 0
    count_nothing = 0
    tsend = 6000.0
    for i in 1:reps
        pattern = cf_returnmap(ep, eff, 1, r, i, tsend, 0.0:1.0:tsend)
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

function prop_canard_rednoise_data(ep, eff, reps)
    r_range = 0.0:0.01:0.9
    data = Vector{Vector{Float64}}(undef,length(r_range))
    @threads for i in eachindex(r_range)
        @inbounds data[i] =  prop_canard_rednoise(ep, eff, r_range[i], reps)
    end

    return prep_data(data, r_range)
end

begin
    rn_ep0079eff05 = prop_canard_rednoise_data(0.079, 0.5, 1000)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/rn_ep0079eff05.csv", rn_ep0079eff05)
    rn_ep0079eff06 = prop_canard_rednoise_data(0.079, 0.6, 1000)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/rn_ep0079eff06.csv", rn_ep0079eff06)
    rn_ep0079eff07 = prop_canard_rednoise_data(0.079, 0.7, 1000)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/rn_ep0079eff07.csv", rn_ep0079eff07)
end

begin
    rn_ep0004eff05 = prop_canard_rednoise_data(0.004, 0.5, 1000)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/rn_ep0004eff05.csv", rn_ep0004eff05)
    rn_ep0004eff06 = prop_canard_rednoise_data(0.004, 0.6, 1000)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/rn_ep0004eff06.csv", rn_ep0004eff06)
    rn_ep0004eff07 = prop_canard_rednoise_data(0.004, 0.7, 1000)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/rn_ep0004eff07.csv", rn_ep0004eff07)
end

