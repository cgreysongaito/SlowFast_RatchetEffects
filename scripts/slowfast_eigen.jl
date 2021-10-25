#Script for slow fast examination of time delays

# Find the real complex divide for different values of efficiency
function findRCdivide_eff(eff)
    par = RozMacPar()
    par.e = eff
    epvals = 0.001:0.0005:1.0

    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        equ = eq_II(par)
        eig1 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        if eig1 < 0 || eig1 > 0
            return epval
            # break
        end
    end
end

# Convert efficiency value to effective proportion real value
function converteff_prop_RCdivide(eff)
    hopf = hopf_rozmac(RozMacPar(e = eff))
    transcrit = transcrit_rozmac(RozMacPar(e = eff))
    hopf_minus_eff = hopf - eff
    propeff = hopf_minus_eff / (hopf-transcrit)
    finalprop = 1 .- propeff
    return finalprop
end

# # 
function findRCdivide_epx(ep, hopf, transcrit)
    par = RozMacPar()
    par.ε = ep
    evals = transcrit:0.00005:hopf

    for (ei, eval) in enumerate(evals)
        par.e = eval
        equ = eq_II(par)
        eig1 = imag.(eigvals(jacmat(roz_mac_II, equ, par))[1])
        if eig1 < 0 || eig1 > 0
            return eval
            break
        end
    end
end

function findRCdivide_epx_data()
    hopf = hopf_rozmac(RozMacPar())
    transcrit = transcrit_rozmac(RozMacPar())
    epvals = 0.001:0.001:100.0
    effRC = zeros(length(epvals))
    @threads for i in eachindex(epvals)
        @inbounds effRC[i] =  findRCdivide_epx(epvals[i], hopf, transcrit)
    end
    effRC_minus_hopf = hopf .- effRC
    effRC_propC = effRC_minus_hopf ./ (hopf-transcrit)
    effRC_propR = 1 .- effRC_propC
    return hcat(collect(1 ./epvals), effRC_propR)
end

function eff_maxeigen_data(ep)
    par = RozMacPar()
    par.ε = ep
    evals = 0.4:0.0001:0.8
    max_eig = fill(0.0, length(evals))

    for (ei, eval) in enumerate(evals)
        par.e = eval
        equ = eq_II(par)
        max_eig[ei] = λ_stability(jacmat(roz_mac_II, equ, par))
    end
    return hcat(evals, max_eig)
end

function imag_epx_data(eff)
    par = RozMacPar()
    par.e = eff
    epvals = 0.01:0.0005:1.0
    imagdata = zeros(length(epvals))
    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        equ = eq_II(par)
        imagdata[epi] = abs(imag.(eigvals(jacmat(roz_mac_II, equ, par))[1]))
    end
    return hcat(collect(1 ./epvals), imagdata)
end

function real_epx_data(eff)
    par = RozMacPar()
    par.e = eff
    epvals = 0.01:0.0005:1.0
    realdata = zeros(length(epvals))
    for (epi, epval) in enumerate(epvals)
        par.ε = epval
        equ = eq_II(par)
        realdata[epi] = λ_stability(jacmat(roz_mac_II, equ, par))
    end
    return hcat(collect(1 ./epvals), realdata)
end