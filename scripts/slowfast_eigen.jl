#Script for slow fast examination of time delays

# When epp is 1, find the excitable and non excitable phases for e - deterministic
function findRCdivide(ep)
    par = RozMacPar()
    par.ε = ep
    evals = 0.441:0.00005:0.9

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


# # create graph of epsilon on x and efficiency value where RC divide by transposing graph of efficiency on x and epsilon value that creates RC divide
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
    effRC
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
