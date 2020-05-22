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
    eq_II_C = r * (a * h * k * (m / (a * (e - h * m))) - a * h * (m / (a * (e - h * m)))^2 + k - m / (a * (e - h * m))) / (a * k)
    return vcat(eq_II_R, eq_II_C)
end


randeq(x) = x * ( 1 + rand(Uniform(1e-7, 1e-6)))

jacmat(model, eq, par) = ForwardDiff.jacobian(eq -> model(eq, par), eq)

λ_stability(M) = maximum(real.(eigvals(M)))
ν_stability(M) = λ_stability((M + M') / 2)

function pert_cb(int)
    int.u[2] = int.u[2] * ( 1 + rand(Normal(0.0, 0.01)))
end


function roz_mac_res(R, C, p)
    @unpack r, k, h, a, m = p
    return r * R * (1 - R / k) - (a * R * C / (1 + a * h * R) )
end

function roz_mac_con(R, C, eff, ep, p)
    @unpack h, a, m = p
    return ep * ( ( eff * a * R * C ) / (1 + a * h * R) - m * C )
end

function con_iso(eff, p)
    @unpack m, a, h = p
    m / (a * (eff - h * m))
end

function res_iso(R, p)
    @unpack a, k, r, h = p
    r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
end

function iso_plot(resrange, par, eff)
    PyPlot.plot(collect(resrange), [res_iso(R, par) for R in resrange])
    return PyPlot.plot(repeat([con_iso(eff, par)], length(resrange)),collect(resrange))
end

function roz_mac_plot(ep, eff)
    resconrange = range(0, stop = 3, length=100)
    U = [roz_mac_res(R, C, par_rozmac) for C in resconrange, R in resconrange]
    V = [roz_mac_con(R, C, eff, ep, par_rozmac) for C in resconrange, R in resconrange]
    speed = sqrt.(U.^2 .+ V.^2)
    lw = 5 .* speed ./ maximum(speed) # Line Widths
    streamplot(collect(resconrange), collect(resconrange), U, V, density = 0.6, color = "k", linewidth = lw)
    return iso_plot(resconrange, par_rozmac, eff)
end
