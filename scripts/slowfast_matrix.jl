include("packages.jl")
include("slowfast_commoncode.jl")
using Reduce
SymPy.@vars R C r a h e m k ε y

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

SymPy.simplify(r * (a * h * k * ( m / (a * (e - h * m))) - a * h * ( m / (a * (e - h * m)))^2 + k - ( m / (a * (e - h * m)))) / (a * k))

SymPy.simplify(subs(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), R), (R, ( m / (a * (e - h * m)))), (C, (e*r*(a*e*k - a*h*k*m - m)/(a^2*k*(e^2 - 2*e*h*m + h^2*m^2))))))

# Type II

#With ε and e
jac11_II = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), R))
jac12_II = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), C))
jac21_II = SymPy.simplify(diff(ε * ( e * a * R * C / (1 + a * h * R) - m * C ), R))
jac22_II = SymPy.simplify(diff(ε * ( e * a * R * C / (1 + a * h * R) - m * C), C))

percapR(jac11_II)
percapR(jac12_II)
percapC(jac21_II)
percapC(jac22_II)

#trace
trace_fun(jac11_II, jac22_II)

#determinant
det_fun(jac11_II, jac22_II, jac12_II, jac21_II)

#Without ε
jac11_IIep = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), R))
jac12_IIep = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), C))
jac21_IIep = SymPy.simplify(diff(e * a * R * C / (1 + a * h * R) - m * C, R))
jac22_IIep = SymPy.simplify(diff(e * a * R * C / (1 + a * h * R) - m * C, C))

percapR(jac11_IIep)
percapR(jac12_IIep)
percapC(jac21_IIep)
percapC(jac22_IIep)

#trace
trace_fun(jac11_IIep, jac22_IIep)

#determinant
det_fun(jac11_IIep, jac22_IIep, jac12_IIep, jac21_IIep)


#Without ep and e
jac11_IIepe = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), R))
jac12_IIepe = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C / (1 + a * h * R), C))
jac21_IIepe = SymPy.simplify(diff(a * R * C / (1 + a * h * R) - m * C, R))
jac22_IIepe = SymPy.simplify(diff(a * R * C / (1 + a * h * R) - m * C, C))

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
jac11_I = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, R))
jac12_I = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, C))
jac21_I = SymPy.simplify(diff(ε * ( e * a * R * C - m * C ), R))
jac22_I = SymPy.simplify(diff(ε * ( e * a * R * C - m * C ), C))

percapR(jac11_I)
percapR(jac12_I)
percapC(jac21_I)
percapC(jac22_I)

#trace
trace_fun(jac11_I, jac22_I)

#determinant
det_fun(jac11_I, jac22_I, jac12_I, jac21_I)

#with e, not ep
jac11_Iep = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, R))
jac12_Iep = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, C))
jac21_Iep = SymPy.simplify(diff(e * a * R * C - m * C, R))
jac22_Iep = SymPy.simplify(diff(e * a * R * C - m * C, C))

percapR(jac11_Iep)
percapR(jac12_Iep)
percapC(jac21_Iep)
percapC(jac22_Iep)

#trace
trace_fun(jac11_Iep, jac22_Iep)

#determinant
det_fun(jac11_Iep, jac22_Iep, jac12_Iep, jac21_Iep)


#with ep, not e
jac11_Ie = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, R))
jac12_Ie = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, C))
jac21_Ie = SymPy.simplify(diff(ε * ( a * R * C - m * C ), R))
jac22_Ie = SymPy.simplify(diff(ε * ( a * R * C - m * C ), C))

percapR(jac11_Ie)
percapR(jac12_Ie)
percapC(jac21_Ie)
percapC(jac22_Ie)

#trace
trace_fun(jac11_Iep, jac22_Iep)

#determinant
det_fun(jac11_Iep, jac22_Iep, jac12_Iep, jac21_Iep)


#witout e and ep
jac11_Iepe = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, R))
jac12_Iepe = SymPy.simplify(diff(r * R * (1 - R / k) - a * R * C, C))
jac21_Iepe = SymPy.simplify(diff(a * R * C - m * C, R))
jac22_Iepe = SymPy.simplify(diff(a * R * C - m * C, C))

percapR(jac11_Iepe)
percapR(jac12_Iepe)
percapC(jac21_Iepe)
percapC(jac22_Iepe)

#trace
trace_fun(jac11_Iepe, jac22_Iepe)

#determinant
det_fun(jac11_Iepe, jac22_Iepe, jac12_Iepe, jac21_Iepe)

eig_ineq(trace_fun(jac11_Iepe, jac22_Iepe), det_fun(jac11_Iepe, jac22_Iepe, jac12_Iepe, jac21_Iepe))
