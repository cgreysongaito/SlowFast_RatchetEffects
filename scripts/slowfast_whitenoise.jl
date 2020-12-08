using Distributed
using CSV
using DataFrames
addprocs(length(Sys.cpu_info())-1)

@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/packages.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_commoncode.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_canardfinder.jl")


## Figure 2 (proportion quasi-canard with different consumer isoclines)
@everywhere function prop_canard_whitenoise(ep, eff, reps, tsend)
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

let
    test = figure()
    pert_phase_plot(0.079, 0.5, 1, 0.9, 5, 20000.0, 2000.0:1.0:20000.0)
    return test
end

function prop_canard_whitenoise_data(eff, iter, tsend)
    epsilon_range = 0.001:0.001:0.15
    data = pmap(ep -> prop_canard_whitenoise(ep, eff, iter, tsend), epsilon_range)
    return prep_data(data, epsilon_range)
end

begin
    wn_eff05_short = prop_canard_whitenoise_data(0.5, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff05_short.csv", wn_eff05_short)
    wn_eff06_short = prop_canard_whitenoise_data(0.6, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff06_short.csv", wn_eff06_short)
    wn_eff07_short = prop_canard_whitenoise_data(0.7, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff07_short.csv", wn_eff07_short)
    wn_eff08_short = prop_canard_whitenoise_data(0.8, 1000, 6000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff08_short.csv", wn_eff08_short)
end

begin
    wn_eff05_long = prop_canard_whitenoise_data(0.5, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff05_long.csv", wn_eff05_long)
    wn_eff06_long = prop_canard_whitenoise_data(0.6, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff06_long.csv", wn_eff06_long)
    wn_eff07_long = prop_canard_whitenoise_data(0.7, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff07_long.csv", wn_eff07_long)
    wn_eff08_long = prop_canard_whitenoise_data(0.8, 1000, 24000.0)
    CSV.write("/home/chrisgg/julia/TimeDelays/data/wn_eff08_long.csv", wn_eff08_long)
end
