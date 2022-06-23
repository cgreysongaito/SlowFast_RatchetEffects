include("packages.jl")
include("slowfast_commoncode.jl")

function calc_min(ep, eff, r, reps, tsend)
    minvaldata = zeros(reps)
    @threads for i in 1:reps
        sol = RozMac_pert(ep, eff, 1.0, r, i, tsend, 9000.0:1.0:tsend)
        @inbounds minvaldata[i] = minimum(sol[1, 1:end])
    end
    return mean(minvaldata)
end

function mindata_ep(eprange, eff, r, reps, tsend)
    avmindata = zeros(length(eprange))
    for (epi, epval) in enumerate(eprange)
        avmindata[epi] = calc_min(epval, eff, r, reps, tsend)
    end
    return reverse(avmindata)
end

let 
    eprangeslow = 0.001:0.01:0.10
    eprangefast = 0.10:0.05:1.0
    eff07whitedata = vcat(mindata_ep(eprangefast, 0.6, 0.0, 100, 10000), mindata_ep(eprangeslow, 0.6, 0.0, 100, 10000))
    test = figure()
    plot( log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangeslow))) , eff07whitedata)
    return test
end


let 
    eprangefast = 0.4:0.01:1.0
    eff07whitedata = mindata_ep(eprangefast, 0.6, 0.0, 100, 10000)
    test = figure()
    plot( log10.(1 ./ reverse(eprangefast)) , eff07whitedata)
    xlabel("Log10 of 1/ϵ")
    ylabel("Average minimum of consumer")
    # return test
    savefig(joinpath(abpath(), "figs/greysongaitoparameters_minimumanalysis_R.pdf"))
end
    