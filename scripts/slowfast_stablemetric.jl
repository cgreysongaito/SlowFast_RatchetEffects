include("packages.jl")
include("slowfast_commoncode.jl")


#min
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

#integral of "residuals"

function residualcalc(datapoint, equil)
    resid = datapoint - equil
    return abs(resid)
end

function integralcalc(conres, ep, eff, r, rep, tsend)
    equil = eq_II(RozMacPar(e=eff))
    sol = RozMac_pert(ep, eff, 1.0, r, rep, tsend, 9000.0:1.0:tsend)
    int = 0
    for i in 1:tsend-9000+1
        if conres == "Consumer"
            int += residualcalc(sol[2, i], equil[2])
        else
            int += residualcalc(sol[1, i], equil[1])
        end
    end
    return int
end

integralcalc("Resource", 1.0, 0.6, 0.0, 1, 10000)

function integraldata(conres, ep, eff, r, reps, tsend)
    intdata = zeros(reps)
    @threads for i in 1:reps
        @inbounds intdata[i] = integralcalc(conres, ep, eff, r, i, tsend)
    end
    return mean(intdata)
end

function epintegraldata(conres, eprange, eff, r, reps, tsend)
    epdata = zeros(length(eprange))
    for (epi, epval) in enumerate(eprange)
        epdata[epi] = integraldata(conres, epval, eff, r, reps, tsend)
    end
    return epdata
end

let
    data = epintegraldata("Resource", 0.8:0.01:1.0, 0.6, 0.0, 5, 10000)
    test = figure()
    plot(0.8:0.01:1.0, data)
    return test
end


# CV
function CVcalc(conres, standard, ep, eff, r, rep, tsend)
    sol = RozMac_pert(ep, eff, 1.0, r, rep, tsend, 2000.0:1.0:tsend)
    if conres == "Consumer"
        cv = std(sol[2,1:end])/mean(sol[2,1:end])
    else
        cv = std(sol[1,1:end])/mean(sol[1,1:end])
    end
    if standard == "Standard"
        return cv/(1/ep)
    else
        return cv
    end
end

function CVdata(conres, standard, ep, eff, r, reps, tsend)
    cvdata = zeros(reps)
    @threads for i in 1:reps
        @inbounds cvdata[i] = CVcalc(conres, standard, ep, eff, r, i, tsend)
    end
    return mean(cvdata)
end

function epcvdata(conres, standard, eprange, eff, r, reps, tsend)
    epdata = zeros(length(eprange))
    for (epi, epval) in enumerate(eprange)
        epdata[epi] = CVdata(conres, standard, epval, eff, r, reps, tsend)
    end
    return reverse(epdata)
end

let 
    eprangeslow = 0.001:0.001:0.05
    eprangemed = 0.05:0.01:0.2
    eprangefast = 0.2:0.05:1.0
    data1CV = vcat(epcvdata("Consumer", "CV", eprangefast, 0.6, 0.0, 5, 10000), epcvdata("Consumer", "CV", eprangemed, 0.6, 0.0, 5, 10000), epcvdata("Consumer", "CV", eprangeslow, 0.6, 0.0, 5, 10000))
    data2CV = vcat(epcvdata("Consumer", "CV", eprangefast, 0.6, 0.0, 5, 35000), epcvdata("Consumer", "CV", eprangemed, 0.6, 0.0, 5, 35000), epcvdata("Consumer", "CV", eprangeslow, 0.6, 0.0, 5, 35000))
    data1standard = vcat(epcvdata("Consumer", "Standard", eprangefast, 0.6, 0.0, 5, 10000), epcvdata("Consumer", "Standard", eprangemed, 0.6, 0.0, 5, 10000), epcvdata("Consumer", "Standard", eprangeslow, 0.6, 0.0, 5, 10000))
    data2standard = vcat(epcvdata("Consumer", "Standard", eprangefast, 0.6, 0.0, 5, 35000), epcvdata("Consumer", "Standard", eprangemed, 0.6, 0.0, 5, 35000), epcvdata("Consumer", "Standard", eprangeslow, 0.6, 0.0, 5, 35000))
    test = figure()
    subplot(2,1,1)
    plot( log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), data1CV)
    plot( log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), data2CV)
    subplot(2,1,2)
    plot( log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), data1standard)
    plot( log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), data2standard)
    xlabel("Log10 of 1/ϵ")
    ylabel("CV/(1/ϵ)")
    tight_layout()
    # return test
    savefig(joinpath(abpath(), "figs/greysongaitoparameters_CV_C.pdf"))
end