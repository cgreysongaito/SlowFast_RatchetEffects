include("packages.jl")
include("slowfast_commoncode.jl")
include("slowfast_canardfinder.jl")

# Figure 3
function prop_canard_rednoise(ep, r)
    res_start = 0.0:0.1:3.0
    con_start = 0.0:0.1:2.5
    count_true = 0
    for (resi, resvalue) in enumerate(res_start)
        for (coni, convalue) in enumerate(con_start)
            if cf_returnmap(ep, 0.55, 1, r, 1234, resvalue, convalue, 5000.0, 2000.0:1.0:5000.0) == true
                count_true += 1
            end
        end
    end
    return count_true / (length(res_start)*length(con_start))
end

function prop_canard_rednoise_data(ep)
    corr_range = 0.0:0.1:0.9
    propcanard = fill(0.0,length(epsilon_range))
    for (corri, corrvalue) in enumerate(corr_range)
        propcanard[corri] = prop_canard_rednoise(ep, corrvalue)
    end
    return propcanard
end

fulldataset = [prop_canard_rednoise_data(0.1), prop_canard_rednoise_data(0.5), prop_canard_rednoise_data(0.8)]

let
    figure2 = figure()
    subplot(1,4,1)
    plot(collect(0.0:0.1:0.9), fulldataset[1])
    subplot(1,4,2)
    plot(collect(0.0:0.1:0.9), fulldataset[2])
    subplot(1,4,3)
    plot(collect(0.0:0.1:0.9), fulldataset[3])
    subplot(1,4,4)
    plot(collect(0.0:0.1:0.9), fulldataset[4])
    return figure2
end




#testing red noise perturbation
using DSP

let
    test = figure()
    plot(0.0:1.0:49, noise(0.8, 50))
    return test
end

function redtest!(du, u, p, t,)
    R, C = u
    du[1] = 0
    du[2] = 0
    return
end

function redtest_pert(freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = RozMacPar()
    noise = noise(r, tsend / freq)
    count = 1
    u0 = [0.0, noise[1]]
    tspan = (0, tsend)

    function pert_cb2(integrator)
        count += 1
        integrator.u[2] = noise[count]
    end

    cb = PeriodicCallback(pert_cb2, freq, initial_affect = false) #as of may 29th - initial_affect does not actually do the affect on the first time point
    prob = ODEProblem(redtest!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return solend = sol(tvals)
end

function redtest_pert_plot(freq, r, seed, tsend, tvals)
    sol = redtest_pert(freq, r, seed, tsend, tvals)
    plot(sol.t, sol.u)
    return ylabel("Resource & \n Consumer Biomass")
end

let
    test = figure()
    pert_timeseries_plot(1.0, -0.8, 1234, 500.0, 200.0:1.0:500.0)
    return test
end

ps = periodogram(redtest_pert(1.0, 0.0, 1234, 500.0, 200.0:1.0:500.0)[2,:])

ps.freq

let
    test = figure()
    plot(ps.freq,ps.power)
    return test
end
