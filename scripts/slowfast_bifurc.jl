#slow fast bifurcation

include("packages.jl")
include("slowfast_commoncode.jl")

# without changing epsilon

par = RozMacPar()
par.ε = 1.0
equ = eq_II(par.e, par)

function con_iso(eff, p)
    @unpack m, a, h = p
    m / (a * (eff - h * m))
end

function res_iso(R, p)
    @unpack a, k, r, h = p
    r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
end

R, C, r, k, a, m, e, h = symbols("R, C, r, k, a, m, e, h", real = true)

f(R, C) = r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)

#Transcritical
#when res_iso equals zero - and what e value creates this
SymPy.solve(f(R,C),R)

g(e) = m / (a * (e - h * m)) - k

SymPy.solve(g(e),e)


function transcrit(p)
    @unpack a, k, h, m = p
    h * m + m / (a * k)
end

transcrit(par)

#when res_iso has a differential of 0 (max)
Calculus.simplify(Calculus.differentiate("r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)", :R))

y(R) = (((r * ((a * h * k - a * h * (2 * R)) - 1)) * (a * k)) / (a * k) ^ 2)

SymPy.solve(y(R), R)

function hopf(p)
    @unpack k, a, h, r, m= p
    R = k / 2 - 1 / (2 * a * h)
    m / (R * a) + ( h * m)
end

hopf(par)
