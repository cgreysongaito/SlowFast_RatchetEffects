using Distributed
using CSV
using DataFrames
addprocs(length(Sys.cpu_info())-1)

@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/packages.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_commoncode.jl")
@everywhere include("/home/chrisgg/julia/TimeDelays/scripts/slowfast_canardfinder.jl")

@everywhere function prop_canard_rednoise(ep, eff, r, reps)
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
    r_range = 0.0:0.1:0.9
    data = pmap(r -> prop_canard_rednoise(ep, eff, r, reps), r_range)
    # for (corri, corrvalue) in enumerate(corr_range)
    #     propcanard[corri] = prop_canard_rednoise(ep, corrvalue)
    # end
    return prep_data(data, r_range)
end

test = prop_canard_rednoise_data(0.004, 0.5, 1000)



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

let
    test = figure()
    pert_phase_plot(0.004, 0.7, 1, 0.9, 11, 6000.0, 0.0:1.0:6000.0)
    return test
end

@everywhere function prop_axialsol_rednoise(ep, eff, r, reps)
    count_true = 0
    tsend = 6000.0
    for i in 1:reps
        data = RozMac_pert(ep, eff, 1, r, i, tsend, 2000.0:1.0:tsend)
        if data.u[end][1] == 3.00 && data.u[end][2] == 0.00
            count_true += 1
        end
    end
    return count_true / reps
end

function prop_axialsol_rednoise_data(ep, eff, reps)
    r_range = 0.0:0.1:0.9
    data = pmap(r -> prop_axialsol_rednoise(ep, eff, r, reps), r_range)
    # for (corri, corrvalue) in enumerate(corr_range)
    #     propcanard[corri] = prop_canard_rednoise(ep, corrvalue)
    # end
    return data
end
#
axial_rn_dataset = [prop_axialsol_rednoise_data(0.004, 0.5,1000), prop_axialsol_rednoise_data(0.004, 0.6,1000), prop_axialsol_rednoise_data(0.004, 0.7,1000)]
#
let
    figure3b = figure()
    subplot(1,3,1)
    plot(collect(0.0:0.1:0.9), axial_rn_dataset[1], color="black")
    fill_between(collect(0.0:0.1:0.9), fill(1.0,10), axial_rn_dataset[1], color="#E66100")
    fill_between(collect(0.0:0.1:0.9), axial_rn_dataset[1], color="#5D3A9B")
    ylabel("Proportion with quasi-canards")
    ylim(0,1)
    title("e = 0.5, ϵ = 0.079")
    subplot(1,3,2)
    plot(collect(0.0:0.1:0.9), axial_rn_dataset[2], color="black")
    fill_between(collect(0.0:0.1:0.9), fill(1.0,10), axial_rn_dataset[2], color="#E66100")
    fill_between(collect(0.0:0.1:0.9), axial_rn_dataset[2], color="#5D3A9B")
    xlabel("Noise correlation (t = -1)")
    ylim(0,1)
    title("e = 0.6, ϵ = 0.079")
    subplot(1,3,3)
    plot(collect(0.0:0.1:0.9), axial_rn_dataset[3], color="black")
    fill_between(collect(0.0:0.1:0.9), fill(1.0,10), axial_rn_dataset[3], color="#E66100")
    fill_between(collect(0.0:0.1:0.9), axial_rn_dataset[3], color="#5D3A9B")
    ylim(0,1)
    title("e = 0.7, ϵ = 0.079")
    tight_layout()
    return figure3b
    # savefig(joinpath(abpath(), "figs/canard_rednoise_prop.png"))
end

#testing red noise perturbation
# using DSP
#
# let
#     test = figure()
#     plot(0.0:1.0:49, noise(0.8, 50))
#     return test
# end
#
# function redtest!(du, u, p, t,)
#     R, C = u
#     du[1] = 0
#     du[2] = 0
#     return
# end
#
# function redtest_pert(freq, r, seed, tsend, tvals)
#     Random.seed!(seed)
#     par = RozMacPar()
#     noise = noise(r, tsend / freq)
#     count = 1
#     u0 = [0.0, noise[1]]
#     tspan = (0, tsend)
#
#     function pert_cb2(integrator)
#         count += 1
#         integrator.u[2] = noise[count]
#     end
#
#     cb = PeriodicCallback(pert_cb2, freq, initial_affect = false) #as of may 29th - initial_affect does not actually do the affect on the first time point
#     prob = ODEProblem(redtest!, u0, tspan, par)
#     sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
#     return solend = sol(tvals)
# end
#
# function redtest_pert_plot(freq, r, seed, tsend, tvals)
#     sol = redtest_pert(freq, r, seed, tsend, tvals)
#     plot(sol.t, sol.u)
#     return ylabel("Resource & \n Consumer Biomass")
# end
#
# let
#     test = figure()
#     pert_timeseries_plot(1.0, -0.8, 1234, 500.0, 200.0:1.0:500.0)
#     return test
# end
#
# ps = periodogram(redtest_pert(1.0, 0.0, 1234, 500.0, 200.0:1.0:500.0)[2,:])
#
# ps.freq
#
# let
#     test = figure()
#     plot(ps.freq,ps.power)
#     return test
# end
