include("packages.jl")

map = 0.0:1.0:10
phase  = 2*pi*rand()
s = sin.((2 * pi * map/maximum(map)) + phase)

s = fill(0.0, length(map))
for i in 1:length(map)
    s[i] = sin(2 * pi * map[i]/maximum(map) + phase)
end



function sine_prep(freq::Int64, map)
    phase = 2*pi*rand()
    s = fill(0.0, length(map))
    for i in 1:length(map)
        s[i] = sin(2 * pi * freq * map[i]/maximum(map) + phase)
    end
    return s
end


let
    map = 0.0:1.0:100
    test = figure()
    plot(map, sine_prep(1, map))
    return test
end

function weighted_sum(amplitudes, noises, map)
    output = fill(0.0, length(map))
    for k in 1:length(noises)
        for x in 1:length(map)
            output[x] += amplitudes[k] * noises[k][x]
        end
    end
    return output
end


amplitudes = [0.1, 0.1, 0.2, 0.3, 0.5, 1.0]
frequencies = [1, 2, 4, 8, 16, 32]

noisestest = [sine_prep(f, 0.0:1.0:100) for f in frequencies]
sumofnoises = weighted_sum(amplitudes, noisestest, 0.0:1.0:100)

let
    map = 0.0:1.0:100
    test = figure()
    plot(map, sumofnoises)
    return test
end



frequencies = 1.0:1.0:30

amplitudes = [f^1 for f in frequencies]

function noise(f_exp, map, freq)
    amplitudes = [f^f_exp for f in freq]
    noises = [sine_prep(f, map) for f in freq]
    sum_of_noises = weighted_sum(amplitudes, noises, map)
    return sum_of_noises
end

let
    map = 0.0:1.0:100
    freq = 1.0:1.0:30
    test = figure()
    plot(map, noise(0, map, freq))
    return test
end


# problem looks like we need to standardise variation


#https://atmos.washington.edu/~breth/classes/AM582/lect/lect8-notes.pdf


function noise2(r, len)
    white = rand(Normal(0.0, 0.01), Int64(len))
    finalnoise = [white[1]]
    for i in 2:Int64(len)
        finalnoise = append!(finalnoise, r * finalnoise[i-1] + white[i] * ( 1 - r^2 )^(1/2))
    end
    return finalnoise
end

noise2(0.8, 50.0)

let
    test = figure()
    plot(0.0:1.0:49, noise2(0.8, 50))
    return test
end




@with_kw mutable struct RozMacPar
    r = 2.0
    k = 3.0
    a = 1.1
    h = 0.8
    e = 0.7
    m = 0.4
    σ = 0.1
    ε = 0.1
end

par_rozmac = RozMacPar()

function roz_mac_II!(du, u, p, t,)
    @unpack r, k, a, h, e, m, ε = p
    R, C = u
    du[1] = r * R * (1 - R / k) - a * R * C / (1 + a * h * R)
    du[2] = ε * ( e * a * R * C / (1 + a * h * R) - m * C )
    return
end

function roz_mac_II(u, par)
    du = similar(u)
    roz_mac_II!(du, u, par, 0.0)
    return du
end

function eq_II(p)
    @unpack r, a, k, h, m, e = p
    eq_II_R = m / (a * (e - h * m))
    eq_II_C = r * (a * h * k * eq_II_R - a * h * eq_II_R^2 + k - eq_II_R) / (a * k)
    return vcat(eq_II_R, eq_II_C)
end




function pert_cb(integrator)
    count += 1
    if isapprox(integrator.u[2], 0.00000000; atol = 1e-8)
        integrator.u[2] = 0.00000000
    else
        integrator.u[2] = maximum([integrator.u[2] + noise[count], 0.0])
    end
end


function RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    noise = noise2(r, tsend / freq)
    count = 1
    u0 = [eq_II(par)[1], eq_II(par)[2] + noise[1]]
    tspan = (0, tsend)

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

RozMac_pert(1.0, 0.6,1.0, 1234, 5000.0, 2000.0:1.0:5000.0)

function pert_timeseries_plot(ep, eff, freq, r, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    plot(sol.t, sol.u)
    return ylabel("Resource & \n Consumer Biomass")
end

let
    test = figure()
    pert_timeseries_plot(1.0, 0.6, 1.0, 0.0, 1234, 5000.0, 2000.0:1.0:5000.0)
    return test
end
