using Hwloc
Hwloc.num_physical_cores()
using .Threads
nthreads()
#Increase the number of threads by going to Julia Extension settings in VS Code or setting the environment variable in bash

using CSV
using DataFrames


include("/home/chrisgg/julia/TimeDelays/scripts/packages.jl")
include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_commoncode.jl")
include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_canardfinder.jl")

function prop_canard_whitenoise(ep, eff, reps, tsend)
    count_canard = 0
    count_axial = 0
    count_nothing = 0
    for i in 1:reps
        pattern = cf_returnmap(ep, eff, 1, 0.0, i, tsend, 0.0:1.0:tsend)
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

function prop_canard_whitenoise_data(eff, iter, tsend)
    epsilon_range = 0.001:0.001:0.15
    data = Vector{Vector{Float64}}(undef,length(epsilon_range))
    @threads for i in eachindex(epsilon_range)
        @inbounds data[i] =  prop_canard_whitenoise(epsilon_range[i], eff, iter, tsend)
    end
    return prep_data(data, epsilon_range)
end

begin
    wn_eff05_short = prop_canard_whitenoise_data(0.5, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff05_short.csv", wn_eff05_short)
    wn_eff06_short = prop_canard_whitenoise_data(0.6, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff06_short.csv", wn_eff06_short)
    wn_eff07_short = prop_canard_whitenoise_data(0.7, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff07_short.csv", wn_eff07_short)
end

begin
    wn_eff05_long = prop_canard_whitenoise_data(0.5, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff05_long.csv", wn_eff05_long)
    wn_eff06_long = prop_canard_whitenoise_data(0.6, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff06_long.csv", wn_eff06_long)
    wn_eff07_long = prop_canard_whitenoise_data(0.7, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff07_long.csv", wn_eff07_long)
end