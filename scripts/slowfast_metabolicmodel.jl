#### Code and functions specific to the Yodzis-Innes metabolic model
#### "Slow organisms exhibit sudden population disappearances in a reddened world" by Greyson-Gaito, Gellner, & McCann.

@with_kw mutable struct YodInnScalePar
    aₜ = 54.9
    aᵣ = 34.3
    rat = 10^-1
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
end #setting up the Yodzis-Innes model

function yod_inn_scale_II(u, par)
    du = similar(u)
    yod_inn_scale_II!(du, u, par, 0.0)
    return du
end #setting up the Yodzis-Innes model

function eq_yodinn_II(p)
    @unpack x, y, δ, k, R₀ = p
    eq_II_R = R₀ / (y - 1)
    eq_II_C = ( ( 1 - δ ) / x ) * eq_II_R * (1 - eq_II_R / k)
    return vcat(eq_II_R, eq_II_C)
end #calculating the interior equilibrium of the Yodzis-Innes model

function con_iso_yodinn(p)
    @unpack R₀, y = p
    R₀/(y - 1)
end #calculating the consumer isocline of the Yodzis-Innes model

function res_iso_yodinn(R, p)
    @unpack x, y, δ, R₀, k  = p
    (R*(R*δ - R + R₀*δ - R₀) + k*(-R*δ + R - R₀*δ + R₀))/(k*x*y)
end #function to help calculate the resource isocline of the Yodzis-Innes model

function iso_plot_YodInn(resrange, par)
    data = [res_iso_yodinn(R, par) for R in resrange]
    plot(collect(resrange), data)
    return vlines(con_iso_yodinn(par), 0, 0.5, colors="orange")
end #plots the isoclines of the Yodzis-Innes model

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
end #function to numerically solve the Yodzis-Innes consumer-resource model with added noise to the consumer

function pert_YodInn_phase_plot(ep, R0, freq, r, seed, tsend, tvals)
    sol = YodInn_pert(ep, R0, freq, r, seed, tsend, tvals)
    plot(sol[1, :], sol[2, :], color = "green")
    xlabel("Resource")
    return ylabel("Consumer")
end #plot a phase space time series for the Yodzis-Innes model

### Below are functions to help calculate the feasible R₀ values, the isocline functions, and the Hopf bifurcation values for the Yodzis-Innes model
# #Feasible R₀
# function feasibleR0(p)
#     @unpack y, k = p
#     lower = ( (y-1)/(y+1) ) * k
#     upper = (y-1) * k
#     return vcat(lower, upper)
# end

# feasibleR0(YodInnScalePar()) #0.699849 1.82557

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

# function hopf_YodInn(p)
#     @unpack k, R₀, δ, x, y= p
#     hopf_R = -R₀/2 + k/2
#     hopf_C = (hopf_R*(hopf_R*δ - hopf_R + R₀*δ - R₀) + k*(-hopf_R*δ + hopf_R - R₀*δ + R₀))/(k*x*y)
#     return vcat(hopf_R,hopf_C)
# end

# hopf_YodInn(YodInnScalePar()) #R=1.1, C=6.651085519687272

