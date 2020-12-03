using Distributed
using CSV
using DataFrames
addprocs(length(Sys.cpu_info())-1)

@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/packages.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_commoncode.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_canardfinder.jl")


## Figure 2 (proportion quasi-canard with different consumer isoclines)
@everywhere function prop_canard(ep, eff, iter, tsend)
    count_true = 0
    for i in 1:iter
        if cf_returnmap(ep, eff, 1, 0.0, i, tsend, 2000.0:1.0:tsend) == true
            count_true += 1
        end
    end
    return count_true / iter
end

let
    test = figure()
    pert_phase_plot(0.079, 0.5, 1, 0.9, 5, 20000.0, 2000.0:1.0:20000.0)
    return test
end

function prop_canard_data(eff, iter, tsend)
    epsilon_range = 0.001:0.001:0.15
    data = pmap(ep -> prop_canard(ep, eff, iter, tsend), epsilon_range)
    return data
end

wn_short_dataset = [prop_canard_data(0.5, 1000, 6000.0), prop_canard_data(0.6, 1000, 6000.0), prop_canard_data(0.7, 1000, 6000.0), prop_canard_data(0.8, 1000, 6000.0)]
wn_long_dataset = [prop_canard_data(0.5, 1000, 12000.0), prop_canard_data(0.6, 1000, 12000.0), prop_canard_data(0.7, 1000, 12000.0), prop_canard_data(0.8, 1000, 12000.0)]
wn_longer_dataset = [prop_canard_data(0.5, 1000, 24000.0), prop_canard_data(0.6, 1000, 24000.0), prop_canard_data(0.7, 1000, 24000.0), prop_canard_data(0.8, 1000, 24000.0)]

CSV.write("/home/chrisgg/julia/TimeDelays/canard_whitenoise_short.csv", DataFrame(prop05 = wn_short_dataset[1], prop06 = wn_short_dataset[2], prop07 = wn_short_dataset[3], prop08 = wn_short_dataset[4]))
CSV.write("/home/chrisgg/julia/TimeDelays/canard_whitenoise_long.csv", DataFrame(prop05 = wn_long_dataset[1], prop06 = wn_long_dataset[2], prop07 = wn_long_dataset[3], prop08 = wn_long_dataset[4]))
CSV.write("/home/chrisgg/julia/TimeDelays/canard_whitenoise_longer.csv", DataFrame(prop05 = wn_longer_dataset[1], prop06 = wn_longer_dataset[2], prop07 = wn_longer_dataset[3], prop08 = wn_longer_dataset[4]))
