include("packages.jl")
include("slowfast_commoncode.jl")

noise_creation(0.9, 1000)

let 
    white = figure()
    plot(0:1:999, noise_creation(0.0, 1000), color = "grey" )
    # return white
    savefig(joinpath(abpath(), "figs/whitenoise_csee.png"))
end

let 
    red = figure()
    plot(0:1:999, noise_creation(0.9, 1000), color = "red" )
    # return red
    savefig(joinpath(abpath(), "figs/rednoise_csee.png"))
end

let
    par = RozMacPar(ε = 1, e = 0.46)
    u0 = [1.5, 1.5]
    tspan = (0.0, 200)
    tvals = 0.0:1.0:200
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    solt = sol(tvals)
    fixed = figure()
    plot(tvals, solt.u)
    # return fixed
    savefig(joinpath(abpath(), "figs/fixedpoint_csee.png"))

end

let
    par = RozMacPar(ε = 1, e = 0.66)
    u0 = [1.5, 1.5]
    tspan = (0.0, 200)
    tvals = 0.0:1.0:200
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    solt = sol(tvals)
    fixed = figure()
    plot(tvals, solt.u)
    # return fixed
    savefig(joinpath(abpath(), "figs/damped_csee.png"))

end

let
    par = RozMacPar(ε = 1, e = 0.8)
    u0 = [eq_II(par)[1]+0.01, eq_II(par)[2]+0.01]
    tspan = (0.0, 200)
    tvals = 0.0:1.0:200
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, reltol = 1e-8)
    solt = sol(tvals)
    fixed = figure()
    plot(tvals, solt.u)
    # return fixed
    savefig(joinpath(abpath(), "figs/cycles_csee.png"))

end