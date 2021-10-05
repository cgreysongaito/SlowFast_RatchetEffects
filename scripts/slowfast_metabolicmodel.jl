include("packages.jl")


@with_kw mutable struct YodInnPar
    aₜ = 54.9
    aᵣ = 34.3
    rat = 10^-6
    δ = 0.55
    fⱼ = 0.99
    aⱼ = 89.2
    k = 3.0
    R₀ = 1.0
    x = ( aₜ / aᵣ ) * rat^0.25
    y = fⱼ * aⱼ / aₜ
end

par_yodinn = YodInnPar()

function yod_inn_II!(du, u, p, t)
    @unpack x, y, δ, k, R₀ = p
    R, C = u
    du[1] = R * (1 - R / k) -  ( x * y / (1 - δ) )  * R * C  / (R₀ + R)
    du[2] = x * y * R * C  / (R₀ + R) - x * C
    return
end

function eq_yodinn_II(p)
    @unpack x, y, δ, k, R₀ = p
    eq_II_R = R₀ / (y - 1)
    eq_II_C = ( ( 1 - δ ) / x ) * eq_II_R * (1 - eq_II_R / k)
    return vcat(eq_II_R, eq_II_C)
end

eq_yodinn_II(par_yodinn)

function con_iso(p)
    @unpack R₀, y = p
    R₀/(y - 1)
end

function res_iso(R, p)
    @unpack x, y, δ, R₀, k  = p
    (R*(R*δ - R + R₀*δ - R₀) + k*(-R*δ + R - R₀*δ + R₀))/(k*x*y)
end

function iso_plot(resrange, par)
    plot(collect(resrange), [res_iso(R, par) for R in resrange])
    return plot(repeat([con_iso(par)], length(resrange)),collect(resrange))
end

let
    test = figure()
    iso_plot(0.0:0.1:3.0, YodInnPar(R₀ = 0.8))
    return test
end



par = YodInnSimPar(rat = 10^-6)
u0 = [0.1, 0.1]
tspan = (0.0, 5000)
tvals = 0.0:1.0:5000
prob = ODEProblem(yod_inn_sim_II!, u0, tspan, par)
sol = DifferentialEquations.solve(prob, reltol = 1e-8)
solt = sol(tvals)


let
    par = YodInnPar(R₀ = 0.71)
    u0 = [0.1, 0.1]
    tspan = (0.0, 2000)
    tvals = 0.0:1.0:2000
    prob = ODEProblem(yod_inn_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    solt = sol(tvals)
    fixed = figure()
    plot(tvals, solt.u)
    return fixed
end

#write code that scales rat by epsilon and then scales fe accordingly
@with_kw mutable struct YodInnScalePar
    aₜ = 54.9
    aᵣ = 34.3
    rat = 10^-6
    δ = 0.55
    fⱼ = 0.99
    aⱼ = 89.2
    k = 3.0
    R₀ = 1.0
    ϵ = 0.1
    fe = 
    x = ( aₜ / aᵣ ) * rat^0.25
    y = fⱼ * aⱼ / aₜ
end

par_yodinnscale = YodInnScalePar()

function yod_inn_scale_II!(du, u, p, t)
    @unpack x, y, δ, k, R₀ = p
    R, C = u
    du[1] = R * (1 - R / k) -  ( x * y / (1 - δ) )  * R * C  / (R₀ + R)
    du[2] = x * y * R * C  / (R₀ + R) - x * C
    return
end

#test with stochasticity if get quasi-canards


#isoclines
using SymPy
@vars R C
@vars x y R₀ k δ

f(R, C) = R * (1 - R / k) - ( x * y / (1 - δ) ) * R * C / ( R₀ + R )
g(R ,C) = x * y * R * C / (R + R₀) - x * C

SymPy.solve(f(R,C),R)
SymPy.solve(f(R,C),C)
SymPy.solve(g(R,C),R)
SymPy.solve(g(R,C),C)






