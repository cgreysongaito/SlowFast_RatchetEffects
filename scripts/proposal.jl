#script to produce figures for time delay section of proposal presentation


include("packages.jl")

@with_kw mutable struct LogisticPar
    r = 2.0
    K = 1.0
end

function logistic!(du, u, p, t)
    @unpack r, K = p
    R = u[1]
    du[1] = r * R * (1 - R / K)
    return
end

let
    t_end = 40
    t_span = (0.0, t_end)

    par = LogisticPar()
    par.r = 2.5

    u0 = [0.1]
    probf = ODEProblem(logistic!, u0, t_span, par)
    solf = DifferentialEquations.solve(probf)

    par.r = 0.5
    probs = ODEProblem(logistic!, u0, t_span, par)
    sols = DifferentialEquations.solve(probs)

    figure()
    plot(solf.t, solf.u)
    plot(sols.t,sols.u)
    savefig("figs/proposallogistic.png")
end

#Lag

@with_kw mutable struct RickerPar
    r = 2.5
    K = 1.0
end

function ricker(N0, n_steps, p)
    @unpack r, K = p
    N = fill(0.0, n_steps)
    N[1] = N0
    for t in 1:(length(N) - 1)
        N[t + 1] = N[t] * exp(r * N[t] * (1 - N[t] / K))
    end
    return N
end

let
    par = RickerPar()
    par.r = 0.5
    outs = ricker(0.1, 100, par)
    par.r = 2.5
    outf = ricker(0.1, 100, par)
    figure()
    plot(range(1,length=100), outf)
    plot(range(1,length=100), outs)
    savefig("figs/proposalricker.png")
end
