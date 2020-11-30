using Distributed
using CSV
using DataFrames
addprocs(length(Sys.cpu_info())-1)

@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/packages.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_commoncode.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_canardfinder.jl")


## Figure 2 (proportion quasi-canard with different consumer isoclines)
@everywhere function prop_canard(ep, eff, iter)
    count_true = 0
    tsend = 6000.0
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

function prop_canard_data(eff, iter)
    epsilon_range1 = 0.001:0.001:0.2
    epsilon_range2 = 0.2:0.1:1.0
    data1 = pmap(ep -> prop_canard(ep, eff, iter), epsilon_range1)
    data2 = pmap(ep -> prop_canard(ep, eff, iter), epsilon_range2)
    return vcat(data1, data2)
end

wn_dataset = [prop_canard_data(0.75, 1000), prop_canard_data(0.85, 1000), prop_canard_data(0.9, 1000), prop_canard_data(0.95, 1000)]

CSV.write("/home/chrisgg/julia/TimeDelays/canard_whitenoise.csv", DataFrame(prop05 = fulldataset[1], prop06 = fulldataset[2], prop07 = fulldataset[3], prop08 = fulldataset[4]))
