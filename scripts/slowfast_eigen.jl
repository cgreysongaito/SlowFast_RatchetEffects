#Script for slow fast examination of time delays
#### "Slow organisms exhibit sudden population disappearances in a reddened world" by Greyson-Gaito, Gellner, & McCann.

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