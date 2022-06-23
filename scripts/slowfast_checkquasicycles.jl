include("packages.jl")
include("slowfast_commoncode.jl")

#maybe i need to increase noise to 0.045 SD because gabes parameters made a smaller resource isocline so the noise was a larger proportion of the resource isocline

let
    data = RozMac_pert(1.0, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    test = figure()
    plot(data.t, data.u)
    return test
end

let
    data = RozMac_pert(0.6, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    test = figure()
    plot(data.t, data.u)
    return test
end

function RozMac_pert_gellner(ep, eff, freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = RozMacPar(r=1.0, k=1.7, h=0.6, m=0.5, a=2.1)
    par.ε = ep
    par.e = eff
    noise = noise_creation(r, tsend / freq)
    count = 1
    u0 = [eq_II(par)[1], eq_II(par)[2] + noise[1]]
    tspan = (0.0, tsend)

    function pert_cb2(integrator)
        count += 1
        if isapprox(integrator.u[2], 0.00000000; atol = 1e-8)
            integrator.u[2] = 0.00000000
        else
            integrator.u[2] = maximum([integrator.u[2] + noise[count], 0.0])
        end
    end

    cb = PeriodicCallback(pert_cb2, freq, initial_affect = false) #as of may 29th - initial_affect does not actually do the affect on the first time point
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return solend = sol(tvals)
end

let
    data = RozMac_pert_gellner(1.0, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    test = figure()
    plot(data.t, data.u)
    return test
end

let
    data = RozMac_pert_gellner(0.6, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    test = figure()
    plot(data.t, data.u)
    return test
end

minimum(RozMac_pert_gellner(1.0, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000))
minimum(RozMac_pert_gellner(0.6, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000))

minimum(RozMac_pert(1.0, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)[2,1:end])
minimum(RozMac_pert(0.5, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)[2,1:end])

#ACF
let
    lrange = 0:1:40
    sol = RozMac_pert(1/1.6, 0.6, 0.1, 0.0, 1, 10000.0, 9000.0:1.0:10000.0)
    ACFdata = autocor(sol[2, 1:end], collect(lrange))
    test = figure()
    plot(lrange, ACFdata)
    return test
end


function quasicycle_data(ep, reps)
    lrange = 0:1:40
    ACFdata = Vector{Vector{Float64}}(undef,reps)
    avprep = zeros(reps)
    avfinal = zeros(length(lrange))
    @threads for i in 1:reps
        sol = RozMac_pert(ep, 0.6, 1.0, 0.0, i, 10000.0, 9000.0:1.0:10000.0)
        @inbounds ACFdata[i] = autocor(sol[2, 1:end], collect(lrange))
    end
    for j in 1:length(lrange)
        for l in 1:reps
        avprep[l] = ACFdata[l][j]
        end
        avfinal[j] = mean(avprep)
    end
    return avfinal
end


let 
    isogellner = figure()
    iso_plot(0.0:0.1:3.0, RozMacPar(r=1.0, k=1.7, h=0.6, e=0.65, m=0.5, a=2.1))
    return isogellner
end


let 
    isogg = figure()
    iso_plot(0.0:0.1:3.0, RozMacPar(e=0.65))
    return isogg
end

#graphing of time series to check quasi canards (when slowing consumer and with real and complex)
#GG parameters
let 
    complexfast = RozMac_pert(1.0, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    complexmedfast = RozMac_pert(0.8, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    complexslow = RozMac_pert(0.6, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    realfast = RozMac_pert(1.0, 0.5, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    realmedfast = RozMac_pert(0.8, 0.5, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    realslow = RozMac_pert(0.6, 0.5, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)

    timeseries = figure()
    subplot(2,3,1)
    plot(realfast.t, realfast[2, 1:end])
    ylim(1.565, 1.71)
    subplot(2,3,2)
    plot(realmedfast.t, realmedfast[2, 1:end])
    ylim(1.565, 1.71)
    subplot(2,3,3)
    plot(realslow.t, realslow[2, 1:end])
    ylim(1.565, 1.71)
    subplot(2,3,4)
    plot(complexfast.t, complexfast[2, 1:end])
    ylim(2.19, 2.335)
    subplot(2,3,5)
    plot(complexmedfast.t, complexmedfast[2, 1:end])
    ylim(2.19, 2.335)
    subplot(2,3,6)
    plot(complexslow.t, complexslow[2, 1:end])
    ylim(2.19, 2.335)
    tight_layout()
    return timeseries
end

let 
    complexfast = RozMac_pert(1.0, 0.65, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    complexmedfast = RozMac_pert(0.8, 0.65, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    complexslow = RozMac_pert(0.6, 0.65, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    realfast = RozMac_pert(1.0, 0.5, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    realmedfast = RozMac_pert(0.8, 0.5, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    realslow = RozMac_pert(0.6, 0.5, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)

    timeseries = figure()
    subplot(2,3,1)
    plot(realfast.t, realfast[2, 1:end])
    ylim(1.565, 1.71)
    subplot(2,3,2)
    plot(realmedfast.t, realmedfast[2, 1:end])
    ylim(1.565, 1.71)
    subplot(2,3,3)
    plot(realslow.t, realslow[2, 1:end])
    ylim(1.565, 1.71)
    subplot(2,3,4)
    plot(complexfast.t, complexfast[2, 1:end])
    ylim(2.19, 2.335)
    subplot(2,3,5)
    plot(complexmedfast.t, complexmedfast[2, 1:end])
    ylim(2.19, 2.335)
    subplot(2,3,6)
    plot(complexslow.t, complexslow[2, 1:end])
    ylim(2.19, 2.335)
    tight_layout()
    return timeseries
end

#Gellner parameters
let 
    complexfast = RozMac_pert_gellner(1.0, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    complexmedfast = RozMac_pert_gellner(0.8, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    complexslow = RozMac_pert_gellner(0.6, 0.65, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    realfast = RozMac_pert_gellner(1.0, 0.5, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    realmedfast = RozMac_pert_gellner(0.8, 0.5, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)
    realslow = RozMac_pert_gellner(0.6, 0.5, 1.0, 0.0, 1, 10000.0, 9000.0:1.0:10000)

    timeseries = figure()
    subplot(2,3,1)
    plot(realfast.t, realfast[2, 1:end])
    ylim(0.265, 0.42)
    subplot(2,3,2)
    plot(realmedfast.t, realmedfast[2, 1:end])
    ylim(0.265, 0.42)
    subplot(2,3,3)
    plot(realslow.t, realslow[2, 1:end])
    ylim(0.265, 0.42)
    subplot(2,3,4)
    plot(complexfast.t, complexfast[2, 1:end])
    ylim(0.455, 0.6)
    subplot(2,3,5)
    plot(complexmedfast.t, complexmedfast[2, 1:end])
    ylim(0.455, 0.6)
    subplot(2,3,6)
    plot(complexslow.t, complexslow[2, 1:end])
    ylim(0.455, 0.6)
    tight_layout()
    return timeseries
end

let 
    complexfast = RozMac_pert_gellner(1.0, 0.65, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    complexmedfast = RozMac_pert_gellner(0.8, 0.65, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    complexslow = RozMac_pert_gellner(0.6, 0.65, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    realfast = RozMac_pert_gellner(1.0, 0.5, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    realmedfast = RozMac_pert_gellner(0.8, 0.5, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)
    realslow = RozMac_pert_gellner(0.6, 0.5, 1.0, 0.0, 125, 10000.0, 9000.0:1.0:10000)

    timeseries = figure()
    subplot(2,3,1)
    plot(realfast.t, realfast[2, 1:end])
    ylim(0.265, 0.42)
    subplot(2,3,2)
    plot(realmedfast.t, realmedfast[2, 1:end])
    ylim(0.265, 0.42)
    subplot(2,3,3)
    plot(realslow.t, realslow[2, 1:end])
    ylim(0.265, 0.42)
    subplot(2,3,4)
    plot(complexfast.t, complexfast[2, 1:end])
    ylim(0.455, 0.6)
    subplot(2,3,5)
    plot(complexmedfast.t, complexmedfast[2, 1:end])
    ylim(0.455, 0.6)
    subplot(2,3,6)
    plot(complexslow.t, complexslow[2, 1:end])
    ylim(0.455, 0.6)
    tight_layout()
    return timeseries
end