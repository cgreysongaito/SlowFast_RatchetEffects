include("packages.jl")
include("slowfast_commoncode.jl")

# CV
function CVcalc(conres, ep, eff, r, rep, tsend)
    sol = RozMac_pert(ep, eff, 1.0, r, rep, tsend, 2000.0:1.0:tsend)
    if isapprox(sol[1,end], 3.0; atol = 1e-2)  && isapprox(sol[2,end], 0.00; atol = 1e-2)
        return missing
    elseif conres == "Consumer"
        return std(sol[2,1:end])/mean(sol[2,1:end])
    else
        return std(sol[1,1:end])/mean(sol[1,1:end])
    end
end

function CVdata(conres, ep, eff, r, reps, tsend)
    cvdata = Vector{Any}(undef, reps)
    @threads for i in 1:reps
        @inbounds cvdata[i] = CVcalc(conres, ep, eff, r, i, tsend)
    end
    return mean(skipmissing(cvdata))
end

function epcvdata(conres, eprange, eff, r, reps, tsend)
    epdata = zeros(length(eprange))
    for (epi, epval) in enumerate(eprange)
        epdata[epi] = CVdata(conres, epval, eff, r, reps, tsend)
    end
    return reverse(epdata)
end

let 
    const eprangeslow = 0.001:0.001:0.05
    const eprangemed = 0.05:0.01:0.2
    cont eprangefast = 0.2:0.05:1.0
    eff = 0.71
    reps = 8
    data2CV = vcat(epcvdata("Consumer", eprangefast, eff, 0.0, reps, 24000), epcvdata("Consumer", eprangemed, eff, 0.0, reps, 24000), epcvdata("Consumer", eprangeslow, eff, 0.0, reps, 24000))
    test = figure()
    plot( log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), data2CV)
    xlabel("Log10 of 1/ϵ")
    ylabel("CV")
    return test
    # savefig(joinpath(abpath(), "figs/greysongaitoparameters_CV_C.pdf"))
end

let 
    eprangeslow = 0.001:0.001:0.05
    eprangemed = 0.05:0.01:0.2
    eprangefast = 0.2:0.05:1.0
    eff = 0.71
    reps = 5
    @time vcat(epcvdata("Consumer", eprangefast, eff, 0.0, reps, 24000), epcvdata("Consumer", eprangemed, eff, 0.0, reps, 24000), epcvdata("Consumer", eprangeslow, eff, 0.0, reps, 24000))
end
