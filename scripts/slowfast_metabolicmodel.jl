# include("packages.jl")
# include("slowfast_commoncode.jl")

@with_kw mutable struct YodInnScalePar
    aₜ = 54.9
    aᵣ = 34.3
    rat = 10^-6
    δ = 0.55
    fⱼ = 0.99
    aⱼ = 89.2
    k = 3.0
    R₀ = 0.8
    ε = 0.1
    x = ( aₜ / aᵣ ) * rat^0.25
    y = fⱼ * aⱼ / aₜ
end

function yod_inn_scale_II!(du, u, p, t)
    @unpack aₜ, aᵣ, y, δ, k, R₀, ε, rat = p
    fₑ = ( (rat * ε) ^ 0.25 ) / ( rat ^ 0.25 )
    xₛ = ( aₜ / aᵣ ) * (rat * ε)^0.25
    R, C = u
    du[1] = R * (1 - R / k) -  ( xₛ * y / ( (1 - δ) * fₑ ) )  * R * C  / (R₀ + R)
    du[2] = ( xₛ * y * R * C  / (R₀ + R) ) - xₛ * C
    return
end

function eq_yodinn_II(p)
    @unpack x, y, δ, k, R₀ = p
    eq_II_R = R₀ / (y - 1)
    eq_II_C = ( ( 1 - δ ) / x ) * eq_II_R * (1 - eq_II_R / k)
    return vcat(eq_II_R, eq_II_C)
end

function con_iso_yodinn(p)
    @unpack R₀, y = p
    R₀/(y - 1)
end

function res_iso_yodinn(R, p)
    @unpack x, y, δ, R₀, k  = p
    (R*(R*δ - R + R₀*δ - R₀) + k*(-R*δ + R - R₀*δ + R₀))/(k*x*y)
end

function YodInn_pert(ep, R0, freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = YodInnScalePar()
    par.ε = ep
    par.R₀ = R0
    noise = noise_creation(r, tsend / freq)
    count = 1
    u0 = [eq_yodinn_II(par)[1], eq_yodinn_II(par)[2] + noise[1]]
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
    prob = ODEProblem(yod_inn_scale_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return solend = sol(tvals)
end

# #Feasible R₀
# function feasibleR0(p)
#     @unpack y, k = p
#     lower = ( (y-1)/(y+1) ) * k
#     upper = (y-1) * k
#     return vcat(lower, upper)
# end

# feasibleR0(YodInnScalePar()) #0.699849 1.82557

# function iso_plot(resrange, par)
#     plot(collect(resrange), [res_iso_yodinn(R, par) for R in resrange])
#     return plot(repeat([con_iso_yodinn(par)], length(resrange)),collect(resrange))
# end

# let
#     isoclines = figure()
#     iso_plot(0.0:0.1:3.0, YodInnScalePar(R₀ = 0.72))
#     return isoclines
# end

# function pert_YodInn_phase_plot(ep, R0, freq, r, seed, tsend, tvals)
#     sol = YodInn_pert(ep, R0, freq, r, seed, tsend, tvals)
#     plot(sol[1, :], sol[2, :])
#     xlabel("Resource")
#     return ylabel("Consumer")
# end

## isoclines
# using SymPy
# @vars R C
# @vars x y R₀ k δ

# f(R, C) = R * (1 - R / k) - ( x * y / (1 - δ) ) * R * C / ( R₀ + R )
# g(R ,C) = x * y * R * C / (R + R₀) - x * C

# SymPy.solve(f(R,C),R)
# SymPy.solve(f(R,C),C)
# SymPy.solve(g(R,C),R)
# SymPy.solve(g(R,C),C)

# #find maximum of resource isocline
# using Calculus
# Calculus.simplify(Calculus.differentiate("(R*(R*δ - R + R₀*δ - R₀) + k*(-R*δ + R - R₀*δ + R₀))/(k*x*y)", :R))
# @vars R C
# @vars x y R₀ k δ a
# r(R) = (((((((R * δ - R) + R₀ * δ) - R₀) + R * (δ - 1)) + k * (1 + -δ)) * (k * x * y)) / (k * x * y) ^ 2)

# SymPy.solve(r(R), R)

# function hopf(p)
#     @unpack k, R₀, δ, x, y= p
#     hopf_R = -R₀/2 + k/2
#     hopf_C = (hopf_R*(hopf_R*δ - hopf_R + R₀*δ - R₀) + k*(-hopf_R*δ + hopf_R - R₀*δ + R₀))/(k*x*y)
#     return vcat(hopf_R,hopf_C)
# end

# hopf(YodInnScalePar()) #R=1.1, C=6.651085519687272

