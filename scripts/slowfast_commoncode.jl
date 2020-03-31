# Slow fast common code
cd("/home/chrisgg/Documents/Guelph/PhD/TimeDelays/")

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


randeq(x) = x * 1 + rand(Uniform(1e-7, 1e-6))

jacmat(model, eq, par) = ForwardDiff.jacobian(eq -> model(eq, par), eq)

λ_stability(M) = maximum(real.(eigvals(M)))
ν_stability(M) = λ_stability((M + M') / 2)
