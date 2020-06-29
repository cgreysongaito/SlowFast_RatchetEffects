include("packages.jl")
include("slowfast_commoncode.jl")

# Singular perturbation of epsilon Roz-Mac model

# Keep Consumer as a parameter and vary parameter to identify bifurcations

@with_kw mutable struct RozMacParC
    r = 2.0
    k = 3.0
    a = 1.1
    h = 0.8
    c = 0.1
end

par_rozmacC = RozMacParC()

function roz_mac_C_II!(du, u, p, t,)
    @unpack r, k, a, h, c = p
    R = u
    du = r * R * (1 - R / k) - a * R * c / (1 + a * h * R)
    return
end

function roz_mac_C_II(u, par)
    du = similar(u)
    roz_mac_C_II!(du, u, par, 0.0)
    return du
end

# Symbolic solving for equilibria
@vars u1 r k a h c

f(u1) = r * u1 * (1 - u1 / k) - a * u1 * c / (1 + a * h * u1)
#-
SymPy.simplify.(SymPy.solve(f(u1), u1))

function eq_C_II(p)
    @unpack r, a, k, h, c = p
    eq_II_R1 = (r*(a*h*k - 1) - sqrt(r*(-4*a^2*c*h*k + a^2*h^2*k^2*r + 2*a*h*k*r + r)))/(2*a*h*r)
    eq_II_R2 = (r*(a*h*k - 1) + sqrt(r*(-4*a^2*c*h*k + a^2*h^2*k^2*r + 2*a*h*k*r + r)))/(2*a*h*r)
    return vcat(0.0, eq_II_R1, eq_II_R2)
end


function eq_C_II_maxC(p)
    @unpack r, a, k, h = p
    eq_II_R1 = (r*(a*h*k - 1) - sqrt(r*(-4*a^2*c*h*k + a^2*h^2*k^2*r + 2*a*h*k*r + r)))/(2*a*h*r)
    return eq_II_R1
end

eq_C_II_maxC(par_rozmacC)

#Jacobian analysis as change Consumer parameter
#One dimensional model so bifurcations possible are saddle-node, transcritical, pitchfork

Calculus.simplify(differentiate("r * u1 * (1 - u1 / k) - a * u1 * c / (1 + a * h * u1)", :u1))

function roz_mac_C_eigen(p)
    @unpack r, a, k, h, c = p
    eq = eq_C_II(p)
    eig = [0.0, 0.0, 0.0]
    for i in 1:3
        R = eq[i]
        eig[i] = ((r * (1 - R / k) + r * R * -(1 / k)) - ((a * c) * (1 + a * h * R) - (a * R * c) * (a * h)) / (1 + a * h * R) ^ 2)
    end
    return hcat(eq,eig)
end

roz_mac_C_eigen(RozMacParC(c = 2.289))

function roz_mac_C_bifurc(par, cvals)
    eq = Any[]
    eig = Any[]

    for (ci, cval) in enumerate(cvals)
        par.c = cval
        equ_eig = roz_mac_C_eigen(par)
        push!(eq, equ_eig[:, 1])
        push!(eig, equ_eig[:, 2])
    end
    return [transpose(convert(Array, VectorOfArray(eq))), transpose(convert(Array, VectorOfArray(eig)))]
end

hcat(-0.1:0.01:2.28, roz_mac_C_bifurc(par_rozmacC, -0.1:0.01:2.28)[2])

function find_trans(cvals)
    data = hcat(cvals, roz_mac_C_bifurc(par_rozmacC, cvals)[2][:,1])

    for i in 1:length(data)
        if data[i, 2] < 0
            return data[i-1, 1]
            break
         #println(epval)
        end
    end
end

find_trans(0.0:0.01:2.28) #1.81 where transcritical happens

let
    fullcvals = 0.0:0.001:2.281
    before_trans_cvals = 0.0:0.001:1.81
    after_trans_cvals = 1.82:0.001:2.281
    sing_pert = figure()
    plot(before_trans_cvals, roz_mac_C_bifurc(par_rozmacC, before_trans_cvals)[1][:,1], color = "black", linestyle = "dashed")
    plot(after_trans_cvals, roz_mac_C_bifurc(par_rozmacC, after_trans_cvals)[1][:,1], color = "black", linestyle = "solid")
    plot(before_trans_cvals, roz_mac_C_bifurc(par_rozmacC, before_trans_cvals)[1][:,2], color = "black", linestyle = "solid")
    plot(after_trans_cvals, roz_mac_C_bifurc(par_rozmacC, after_trans_cvals)[1][:,2], color = "black", linestyle = "dashed")
    plot(fullcvals, roz_mac_C_bifurc(par_rozmacC, fullcvals)[1][:,3], color  = "black")
    plot(fullcvals, fill(con_iso(RozMacPar(e = 0.46)), length(fullcvals)))
    ylabel("Resource")
    xlabel("Consumer")
    return sing_pert
end

#Prediction transcritical when C = 0 Incorrect (that is for changing e)
#Should be bifurcation at C = max Resource isocline  - maybe CORRECT

#Transcritical at C = 1.81
#Saddle node (fold) bifurcation at C = 2.281


#Identify where consumer increasing and decreasing
#impacted by e because consumer isocline dictates whether dc/dt + or -   (m / (a * (e - h * m)))
