include("packages.jl")
include("slowfast_commoncode.jl")

@vars R C r a h e m k ε

function percapR(equ)
    return eval(Expr(:call, :(/), equ, :R))
end

function percapC(equ)
    return eval(Expr(:call, :(/), equ, :C))
end

# Type II

#With ε
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :R)))
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :C)))
Calculus.simplify(percapC(differentiate("ε * ( e * a * R * C / (1 + a * h * R) - m * C )", :R)))
Calculus.simplify(percapC(differentiate("ε * ( e * a * R * C / (1 + a * h * R) - m * C )", :C)))

#Without ε
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :R)))
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :C)))
Calculus.simplify(percapC(differentiate("e * a * R * C / (1 + a * h * R) - m * C", :R)))
Calculus.simplify(percapC(differentiate("e * a * R * C / (1 + a * h * R) - m * C", :C)))

#Without ep and e
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :R)))
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :C)))
Calculus.simplify(percapC(differentiate("a * R * C / (1 + a * h * R) - m * C", :R)))
Calculus.simplify(percapC(differentiate("a * R * C / (1 + a * h * R) - m * C", :C)))

function jac_rozmac(p, R, C)
    @vars ε
    @unpack r, k, a, h, e, m = p
    dfdr = ((r * (1 - R / k) + r * R * -(1 / k)) - ((a * C) * (1 + a * h * R) - (a * R * C) * (a * h)) / (1 + a * h * R) ^ 2)
    dfdc = -(((a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2)
    dgdr = (ε * (((e * a * C) * (1 + a * h * R) - (e * a * R * C) * (a * h)) / (1 + a * h * R) ^ 2))
    dgdc = (ε * (((e * a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2 - m))

    return [dfdr dfdc; dgdr dgdc]
end

jac_rozmac(par_rozmac, eq_II(par_rozmac)[1], eq_II(par_rozmac)[2])

jac_rozmac(par_rozmac, 0, 0)

function jac_rozmac_woep(p, R, C)
    @vars ε
    @unpack r, k, a, h, e, m = p
    dfdr = (r * (1 - R / k) + r * R * -(1 / k)) - ((a * C) * (1 + a * h * R) - (a * R * C) * (a * h)) / (1 + a * h * R) ^ 2
    dfdc = -(((a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2)
    dgdr = ((e * a * C) * (1 + a * h * R) - (e * a * R * C) * (a * h)) / (1 + a * h * R) ^ 2
    dgdc = ((e * a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2 - m

    return [dfdr dfdc; dgdr dgdc]
end

jac_rozmac_woep(par_rozmac, eq_II(par_rozmac)[1], eq_II(par_rozmac)[2])

jac_rozmac_woep(par_rozmac, 0, 0)

function jac_rozmac_woepe(p)
    @vars ε R C
    @unpack r, k, a, h, m = p
    dfdr = (r * (1 - R / k) + r * R * -(1 / k)) - ((a * C) * (1 + a * h * R) - (a * R * C) * (a * h)) / (1 + a * h * R) ^ 2
    dfdc = -(((a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2)
    dgdr = ((a * C) * (1 + a * h * R) - (a * R * C) * (a * h)) / (1 + a * h * R) ^ 2
    dgdc = ((a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2 - m

    return [dfdr dfdc; dgdr dgdc]
end

jac_rozmac_woepe(par_rozmac)

jac_rozmac_woepe(par_rozmac, 0, 0)

#Type I
# with e and ep
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C", :R)))
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C", :C)))
Calculus.simplify(percapC(differentiate("ε * ( e * a * R * C - m * C )", :R)))
Calculus.simplify(percapC(differentiate("ε * ( e * a * R * C - m * C )", :C)))

#with e, not ep,
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C", :R)))
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C", :C)))
Calculus.simplify(percapC(differentiate("e * a * R * C - m * C", :R)))
Calculus.simplify(percapC(differentiate("e * a * R * C - m * C", :C)))

#witout e and ep
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C", :R)))
Calculus.simplify(percapR(differentiate("r * R * (1 - R / k) - a * R * C", :C)))
Calculus.simplify(percapC(differentiate("a * R * C - m * C", :R)))
Calculus.simplify(percapC(differentiate("a * R * C - m * C", :C)))
