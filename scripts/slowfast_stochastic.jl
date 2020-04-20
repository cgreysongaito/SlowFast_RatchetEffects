include("packages.jl")
include("slowfast_commoncode.jl")

function sto_cv_plot(ep, cvboth)
    par = RozMacPar()
    par.ε = ep
    evals = 0.4:0.005:0.9
    u0 = eq_II(par)
    tspan = (0.0, 10000.0)
    cv = fill(0.0, length(evals))
    mn = fill(0.0, length(evals))
    sd = fill(0.0, length(evals))
    function pert_cb(int)
        int.u[2] = int.u[2] * ( 1 + rand(Normal(0.0, 0.01)))
    end
    cb = PeriodicCallback(pert_cb, 1, initial_affect = true)

    for (ei, eval) in enumerate(evals)
        par.e = eval
        prob = ODEProblem(roz_mac_II!, u0, tspan, par)
        sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)

        cv[ei] = std(sol[1, 1:end] ./ mean(sol[2, 1:end]))
        mn[ei] = mean(sol[2, 1:end])
        sd[ei] = std(sol[2, 1:end])
    end

    if cvboth == "cv"
        return PyPlot.plot(collect(evals), cv)
    else
        PyPlot.plot(collect(evals), mn, label = "mean")
        return PyPlot.plot(collect(evals), sd, label = "sd")
    end
end

let
    figure()
    subplot(211)
    sto_cv_plot(1, "cv")
    ylabel("Consumer CV")
    subplot(212)
    sto_cv_plot(1, "sd")
    xlabel("Efficiency (e)")
    ylabel("Mean (blue), SD (orange)")
    gcf()
    #savefig("figs/cv_afterhopf_plot.png")
end

let
    figure()
    subplot(211)
    sto_cv_plot(0.01, "cv")
    ylabel("Consumer CV")
    subplot(212)
    sto_cv_plot(0.01, "sd")
    xlabel("Efficiency (e)")
    ylabel("Mean (blue), SD (orange)")
    gcf()
    #savefig("figs/cv_afterhopf_plot.png")
end

#joinpath(abpath(), "figs/figure1.png")
