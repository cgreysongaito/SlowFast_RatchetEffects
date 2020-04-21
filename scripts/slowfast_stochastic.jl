include("packages.jl")
include("slowfast_commoncode.jl")

##
function pert_cb(int)
    int.u[2] = int.u[2] * ( 1 + rand(Normal(0.0, 0.01)))
end

function sto_cv_plot(ep, save)
    par = RozMacPar()
    par.ε = ep
    evals = 0.441:0.005:0.9
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)
    cv = fill(0.0, length(evals))
    mn = fill(0.0, length(evals))
    sd = fill(0.0, length(evals))
    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)

    for (ei, eval) in enumerate(evals)
        par.e = eval
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)

        cv[ei] = std(sol[1, 1:end] ./ mean(sol[2, 1:end]))
        mn[ei] = mean(sol[2, 1:end])
        sd[ei] = std(sol[2, 1:end])
    end

    let
        figure()
        subplot(211)
        PyPlot.plot(collect(evals), cv)
        ylabel("Consumer CV")
        subplot(212)
        PyPlot.plot(collect(evals), mn, label = "mean")
        PyPlot.plot(collect(evals), sd, label = "sd")
        xlabel("Efficiency (e)")
        ylabel("Mean (blue), SD (orange)")
        if save == "save" && ep < 1
            savefig(joinpath(abpath(), "figs/sto_eptiny_cv_plot.png"))
        elseif save == "save"
            savefig(joinpath(abpath(), "figs/sto_ep1_cv_plot.png"))
        else
        return gcf()
        end
    end
end

sto_cv_plot(1, "save")
sto_cv_plot(1, "show")
sto_cv_plot(0.01, "save")
sto_cv_plot(0.01, "show")

function numsolvplot(u0, tspan, par, ep)
    par.ε = ep
    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return PyPlot.plot(sol.t, sol.u)
end

function stoepboth(eff)
    par = RozMacPar()
    par.e = eff
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)

    let
        figure()
        subplot(211)
        numsolvplot(u0, tspan, par, 1)
        subplot(212)
        numsolvplot(u0, tspan, par, 0.01)
        return gcf()
    end
end

stoepboth(0.5)

function stoepsingle(eff,ep)
    par = RozMacPar()
    par.e = eff
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)
    numsolvplot(u0, tspan, par, ep)
end

let
    figure()
    stoepsingle(0.45, 0.01)
    gcf()
end




##### Before Hopf , starting values large difference from equilibrium with small epsilon
function randstart(ep, eff)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    tspan = (0.0, 10000000.0)

    uC = [0.5, 1.0, 2.0, 2.9]
    uR = [0.5, 1.0, 2.0, 2.9]
    for i in 1:4, j in 1:4
        u0 = [uC[i], uR[j]]
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
        PyPlot.plot(sol.t, sol[2,1:end])
    end
end

let
    figure()
    randstart(0.0001,0.6)
    gcf()
end


par = RozMacPar()
par.ε = 0.01
par.e = 0.6
tspan = (0.0, 1000.0)

u0 = [2.9, 2.9]

prob = ODEProblem(roz_mac_II!, u0, tspan, par)
sol = DifferentialEquations.solve(prob, reltol = 1e-8, maxiters = 1e9)
let
    figure()
    PyPlot.plot(sol.t, sol[2,1:end])
    gcf()
end

#something is not working when epsilon is tiny and consumer is above the "hopf" - consumer reduces to zero and doesn't increase again and do't have enough maxiters (or takes forever)
