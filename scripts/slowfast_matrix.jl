include("packages.jl")
include("slowfast_commoncode.jl")

@vars R C r a h e m k ε

function percapR(equ)
    return SymPy.simplify(eval(Expr(:call, :(/), equ, :R)))
end

function percapC(equ)
    return SymPy.simplify(eval(Expr(:call, :(/), equ, :C)))
end

function trace_fun(eq11, eq22)
    return SymPy.simplify(eval(Expr(:call, :(+), eq11, eq22)))
end

function det_fun(eq11, eq22, eq12, eq21)
    diag = Expr(:call, :(*), eq11, eq22)
    offdiag = Expr(:call, :(*), eq12, eq21)
    return SymPy.simplify(eval(Expr(:call, :(-), eq11, eq22)))
end


function eig_ineq(tra, dete)
    tra2 = Expr(:call, :(*), tra, tra)
    dete4 = Expr(:call, :(*), 4, dete)
    return SymPy.simplify(eval(Expr(:call, :(-), tra2, dete4)))
end

# Type II

#With ε and e
jac11_II = differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :R)
jac12_II = differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :C)
jac21_II = differentiate("ε * ( e * a * R * C / (1 + a * h * R) - m * C )", :R)
jac22_II = differentiate("ε * ( e * a * R * C / (1 + a * h * R) - m * C )", :C)

percapR(jac11_II)
percapR(jac12_II)
percapC(jac21_II)
percapC(jac22_II)

#trace
trace_fun(jac11_II, jac22_II)

#determinant
det_fun(jac11_II, jac22_II, jac12_II, jac21_II)

#Without ε
jac11_IIep = differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :R)
jac12_IIep = differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :C)
jac21_IIep = differentiate("e * a * R * C / (1 + a * h * R) - m * C", :R)
jac22_IIep = differentiate("e * a * R * C / (1 + a * h * R) - m * C", :C)

percapR(jac11_IIep)
percapR(jac12_IIep)
percapC(jac21_IIep)
percapC(jac22_IIep)

#trace
trace_fun(jac11_IIep, jac22_IIep)

#determinant
det_fun(jac11_IIep, jac22_IIep, jac12_IIep, jac21_IIep)


#Without ep and e
jac11_IIepe = differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :R)
jac12_IIepe = differentiate("r * R * (1 - R / k) - a * R * C / (1 + a * h * R)", :C)
jac21_IIepe = differentiate("a * R * C / (1 + a * h * R) - m * C", :R)
jac22_IIepe = differentiate("a * R * C / (1 + a * h * R) - m * C", :C)

percapR(jac11_IIepe)
percapR(jac12_IIepe)
percapC(jac21_IIepe)
percapC(jac22_IIepe)

#trace
trace_fun(jac11_IIepe, jac22_IIepe)

#determinant
det_fun(jac11_IIepe, jac22_IIepe, jac12_IIepe, jac21_IIepe)


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
#
# function jac_rozmac_woep(p, R, C)
#     @vars ε
#     @unpack r, k, a, h, e, m = p
#     dfdr = (r * (1 - R / k) + r * R * -(1 / k)) - ((a * C) * (1 + a * h * R) - (a * R * C) * (a * h)) / (1 + a * h * R) ^ 2
#     dfdc = -(((a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2)
#     dgdr = ((e * a * C) * (1 + a * h * R) - (e * a * R * C) * (a * h)) / (1 + a * h * R) ^ 2
#     dgdc = ((e * a * R) * (1 + a * h * R)) / (1 + a * h * R) ^ 2 - m
#
#     return [dfdr dfdc; dgdr dgdc]
# end
#
# jac_rozmac_woep(par_rozmac, eq_II(par_rozmac)[1], eq_II(par_rozmac)[2])
#
# jac_rozmac_woep(par_rozmac, 0, 0)
#
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

# Type I
# with e and ep
jac11_I = differentiate("r * R * (1 - R / k) - a * R * C", :R)
jac12_I = differentiate("r * R * (1 - R / k) - a * R * C", :C)
jac21_I = differentiate("ε * ( e * a * R * C - m * C )", :R)
jac22_I = differentiate("ε * ( e * a * R * C - m * C )", :C)

percapR(jac11_I)
percapR(jac12_I)
percapC(jac21_I)
percapC(jac22_I)

#trace
trace_fun(jac11_I, jac22_I)

#determinant
det_fun(jac11_I, jac22_I, jac12_I, jac21_I)

#with e, not ep
jac11_Iep = differentiate("r * R * (1 - R / k) - a * R * C", :R)
jac12_Iep = differentiate("r * R * (1 - R / k) - a * R * C", :C)
jac21_Iep = differentiate("e * a * R * C - m * C", :R)
jac22_Iep = differentiate("e * a * R * C - m * C", :C)

percapR(jac11_Iep)
percapR(jac12_Iep)
percapC(jac21_Iep)
percapC(jac22_Iep)

#trace
trace_fun(jac11_Iep, jac22_Iep)

#determinant
det_fun(jac11_Iep, jac22_Iep, jac12_Iep, jac21_Iep)


#with ep, not e
jac11_Ie = differentiate("r * R * (1 - R / k) - a * R * C", :R)
jac12_Ie = differentiate("r * R * (1 - R / k) - a * R * C", :C)
jac21_Ie = differentiate("ε * ( a * R * C - m * C )", :R)
jac22_Ie = differentiate("ε * ( a * R * C - m * C )", :C)

percapR(jac11_Ie)
percapR(jac12_Ie)
percapC(jac21_Ie)
percapC(jac22_Ie)

#trace
trace_fun(jac11_Iep, jac22_Iep)

#determinant
det_fun(jac11_Iep, jac22_Iep, jac12_Iep, jac21_Iep)


#witout e and ep
jac11_Iepe = Calculus.simplify(differentiate("r * R * (1 - R / k) - a * R * C", :R))
jac12_Iepe = Calculus.simplify(differentiate("r * R * (1 - R / k) - a * R * C", :C))
jac21_Iepe = Calculus.simplify(differentiate("a * R * C - m * C", :R))
jac22_Iepe = Calculus.simplify(differentiate("a * R * C - m * C", :C))

percapR(jac11_Iepe)
percapR(jac12_Iepe)
percapC(jac21_Iepe)
percapC(jac22_Iepe)

#trace
trace_fun(jac11_Iepe, jac22_Iepe)

#determinant
det_fun(jac11_Iepe, jac22_Iepe, jac12_Iepe, jac21_Iepe)

eig_ineq(trace_fun(jac11_Iepe, jac22_Iepe), det_fun(jac11_Iepe, jac22_Iepe, jac12_Iepe, jac21_Iepe))
