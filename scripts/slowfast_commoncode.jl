# Slow fast common code

function abpath()
    replace(@__DIR__, "scripts" => "")
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
    μ  = 0.0
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


randeq(x) = x * ( 1 + rand(Uniform(1e-7, 1e-6)))

jacmat(model, eq, par) = ForwardDiff.jacobian(eq -> model(eq, par), eq)

λ_stability(M) = maximum(real.(eigvals(M)))
ν_stability(M) = λ_stability((M + M') / 2)

function roz_mac_res(R, C, p)
    @unpack r, k, h, a, m = p
    return r * R * (1 - R / k) - (a * R * C / (1 + a * h * R) )
end

function roz_mac_con(R, C, eff, ep, p)
    @unpack h, a, m = p
    return ep * ( ( eff * a * R * C ) / (1 + a * h * R) - m * C )
end

function con_iso(p)
    @unpack m, a, h, e = p
    m / (a * (e - h * m))
end

function res_iso(R, p)
    @unpack a, k, r, h = p
    r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
end

function iso_plot(resrange, par)
    plot(collect(resrange), [res_iso(R, par) for R in resrange])
    return plot(repeat([con_iso(par)], length(resrange)),collect(resrange))
end

function roz_mac_plot(ep, eff)
    resconrange = range(0, stop = 3, length=100)
    U = [roz_mac_res(R, C, par_rozmac) for C in resconrange, R in resconrange]
    V = [roz_mac_con(R, C, eff, ep, par_rozmac) for C in resconrange, R in resconrange]
    speed = sqrt.(U.^2 .+ V.^2)
    lw = 5 .* speed ./ maximum(speed) # Line Widths
    streamplot(collect(resconrange), collect(resconrange), U, V, density = 0.6, color = "k", linewidth = lw)
    return iso_plot(resconrange, RozMacPar(e = eff))
end

function pert_cb(integrator)
    if isapprox(integrator.u[2], 0.00000000; atol = 1e-8)
        integrator.u[2] = 0.00000000
    else
        integrator.u[2] = maximum([integrator.u[2] + rand(Normal(integrator.p.μ, 0.01)), 0.0])
    end
end


function RozMac_pert(ep, eff, mean, freq, seed, tsend, tvals)
    Random.seed!(seed)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    par.μ = mean
    u0 = [eq_II(par)[1], eq_II(par)[2] + rand(Normal(mean, 0.01))]
    tspan = (0, tsend)
    cb = PeriodicCallback(pert_cb, freq, initial_affect = false) #as of may 29th - initial_affect does not actually do the affect on the first time point
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return solend = sol(tvals)
end


function pert_timeseries_plot(ep, eff, mean, freq, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, mean, freq, seed, tsend, tvals)
    plot(sol.t, sol.u)
    return ylabel("Resource & \n Consumer Biomass")
end

function pert_consumer_timeseries_plot(ep, eff, mean, freq, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, mean, freq, seed, tsend, tvals)
    plot(sol.t, sol[2, :])
    return ylabel("Consumer biomass")
end

function pert_phase_plot(ep, eff, mean, freq, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, mean, freq, seed, tsend, tvals)
    plot(sol[1, :], sol[2, :])
    xlabel("Resource")
    return ylabel("Consumer")
end
