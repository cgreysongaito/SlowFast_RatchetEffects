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
function findRCdivide_epx(ep)
    par = RozMacPar()
    par.ε = ep
    evals = 0.441:0.00005:0.71

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
    epvals = 0.01:0.005:2.0
    effRC = fill(0.0, length(epvals))
    for (epi, epval) in enumerate(epvals)
        effRC[epi] = findRCdivide_epx(epval)
    end
    effRC_minus_hopf = 0.71 .- effRC
    effRC_propC = effRC_minus_hopf ./ (0.71-0.441)
    effRC_propR = 1 .- effRC_propC
    return hcat(collect(epvals), effRC_propR)
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
