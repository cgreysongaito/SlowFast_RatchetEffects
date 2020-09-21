#script to produce figures for time delay section of proposal presentation


include("packages.jl")

function abpath()
    replace(@__DIR__, "scripts" => "")
end

@with_kw mutable struct LogisticPar
    r = 2.0
    K = 1.0
    τ = 0.75
end

function logistic!(du, u, p, t)
    @unpack r, K = p
    R = u[1]
    du[1] = r * R * (1 - R / K)
    return
end

function logistic_lag!(du, u, h, p, t)
    @unpack r, K, τ = p
    hist = h(p, t-τ)[1]
    du[1] = r * u[1] * (1 - hist / K)
    return
end





let
    t_end = 20
    t_span = (0.0, t_end)

    par = LogisticPar()
    par.r = 1.5

    u0 = [0.1]
    probf = ODEProblem(logistic!, u0, t_span, par)
    solf = DifferentialEquations.solve(probf)

    par.r = 0.5
    probs = ODEProblem(logistic!, u0, t_span, par)
    sols = DifferentialEquations.solve(probs)

    test = figure()
    plot(solf.t, solf.u, label = "r = 1.5")
    plot(sols.t,sols.u, label = "r = 0.5")
    xlabel("Time", fontsize = 15)
    ylabel("Biomass", fontsize = 15)
    xticks(fontsize = 12)
    yticks(fontsize = 12)
    ylim(0.0,1.5)
    legend(fontsize = 15)
    legend(fontsize = 15)
    # return test
    savefig(joinpath(abpath(), "figs/proposallogistic.png"))
end

#Lag

let
    t_end = 20
    t_span = (0.0, t_end)
    h(p, t) = 0.0

    par = LogisticPar()
    par.r = 1.5
    lags = [par.τ]
    u0 = [0.1]
    probf = DDEProblem(logistic_lag!, u0, h, t_span, par; constant_lags = lags)
    solf = DifferentialEquations.solve(probf)

    par.r = 0.5
    probs = DDEProblem(logistic_lag!, u0, h, t_span, par; constant_lags = lags)
    sols = DifferentialEquations.solve(probs)

    test = figure()
    plot(solf.t, solf.u, label = "r = 1.5")
    plot(sols.t,sols.u, label = "r = 0.5")
    xlabel("Time", fontsize = 15)
    ylabel("Biomass", fontsize = 15)
    xticks(fontsize = 12)
    yticks(fontsize = 12)
    ylim(0.0,1.5)
    legend(fontsize = 15)
    # return test
    savefig(joinpath(abpath(),"figs/proposallogisticlag.png"))
end



#Discrete


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


## Functional responses

@with_kw mutable struct funresPar
    a = 1.1
    R0 = 0.5
    h = 0.8
end

parfunres = funresPar()
function fun_resIIa(p, rvals)
    @unpack a, h = p
    rate = [ ( a * i ) / ( 1 + a * h * i )  for i in rvals]
    return rate
end


function denom_fun_resIIa(p, rvals)
    @unpack a, h = p
    rate = [  1 + a * h * i   for i in rvals]
    return rate
end

let
    rvals = 0.0:0.1:4
    fun_res = figure()
    plot(rvals, denom_fun_resIIa(funresPar(h = 0.5), rvals), label = "0.5")
    plot(rvals, denom_fun_resIIa(funresPar(h = 0.8), rvals), label = "0.8")
    plot(rvals, denom_fun_resIIa(funresPar(h = 1.0), rvals), label = "1.0")
    legend()
    return fun_res
end

function fun_resIIb(p, rvals)
    @unpack a, R0 = p
    rate = [ ( a * i ) / ( R0 + i )  for i in rvals]
    return rate
end

function denom_fun_resIIb(p, rvals)
    @unpack a, R0 = p
    rate = [  R0 + i   for i in rvals]
    return rate
end


denom_fun_resIIb(funresPar(R0 = 0.5), 0.0:0.1:8)

let
    rvals = 0.0:0.1:8
    fun_res = figure()
    plot(rvals, [  0.5 * i   for i in rvals], label = "0.5")
    plot(rvals, [  0.8 * i   for i in rvals], label = "0.8")
    plot(rvals, [  1.0 * i   for i in rvals], label = "1.0")
    plot(rvals, denom_fun_resIIb(funresPar(R0 = 0.5), rvals), label = "denom")
    # plot(rvals, denom_fun_resIIb(funresPar(R0 = 0.8), rvals), label = "0.8")
    # plot(rvals, denom_fun_resIIb(funresPar(R0 = 1.0), rvals), label = "1.0")
    legend()
    return fun_res
end

function fun_resIII(p, rvals)
    @unpack a, R0 = p
    rate = [ ( a * i^2 ) / ( R0^2 + i^2 )  for i in rvals]
    return rate
end

let
    rvals = 0.0:0.1:10
    fun_res = figure()
    plot(rvals, fun_resIIa(funresPar(a = 0.5), rvals), label = "0.5")
    plot(rvals, fun_resIIa(funresPar(a = 0.8), rvals), label = "0.8")
    plot(rvals, fun_resIIa(funresPar(a = 1.0), rvals), label = "1.0")
    legend()
    return fun_res
end

let
    rvals = 0.0:0.1:15
    fun_res = figure()
    plot(rvals, fun_resIIb(funresPar(a = 0.5), rvals), label = "0.5")
    plot(rvals, fun_resIIb(funresPar(a = 0.8), rvals), label = "0.8")
    plot(rvals, fun_resIIb(funresPar(a = 1.0), rvals), label = "1.0")
    legend()
    return fun_res
end

let
    rvals = 0.0:0.1:3
    fun_res = figure()
    plot(rvals, fun_resIII(funresPar(R0 = 0.5), rvals), label = "0.5")
    plot(rvals, fun_resIII(funresPar(R0 = 0.8), rvals), label = "0.8")
    plot(rvals, fun_resIII(funresPar(R0 = 1.0), rvals), label = "1.0")
    legend()
    return fun_res
end

@vars R a h

SymPy.solve((2 * a * R) / (1 + a * h * R) - 1 / h, R)
