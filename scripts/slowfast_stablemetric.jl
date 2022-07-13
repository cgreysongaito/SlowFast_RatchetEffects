include("packages.jl")
include("slowfast_commoncode.jl")

# CV
function CVcalc(ep, eff, r, rep, tsend)
    sol = RozMac_pert(ep, eff, 1.0, r, rep, tsend, 2000.0:10.0:tsend)
    if isapprox(sol[1,end], 3.0; atol = 1e-2)  && isapprox(sol[2,end], 0.00; atol = 1e-2)
        return missing
    else
        return std(sol[2,1:end])/mean(sol[2,1:end])
    end
end

function CVdata(ep, eff, r, reps, tsend)
    cvdata = Union{Missing, Float64}[]
    for i in 1:reps
        push!(cvdata, CVcalc(ep, eff, r, i, tsend))
    end
    return mean(skipmissing(cvdata))
end

function epcvdata(eprange, eff, r, reps, tsend)
    epdata = zeros(length(eprange))
    @threads for i in eachindex(eprange)
        @inbounds epdata[i] = CVdata(eprange[i], eff, r, reps, tsend)
    end
    return reverse(epdata)
end

eprangeslow = 0.001:0.001:0.05;
eprangemed = 0.05:0.01:0.2;
eprangefast = 0.2:0.05:1.0;
eff_real = 0.5;
eff_complex = 0.71;
reps = 5;

let 
    CV_real = vcat(epcvdata(eprangefast, eff_real, 0.0, reps, 24000), epcvdata(eprangemed, eff_real, 0.0, reps, 24000), epcvdata(eprangeslow, eff_real, 0.0, reps, 24000))
    dataCV_real = DataFrame(eprange = log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), CV = CV_real)
    CSV.write(joinpath(abpath(), "data/CVresult_eff05.csv"), dataCV_real)
end

let
    CV_complex = vcat(epcvdata(eprangefast, eff_complex, 0.0, reps, 24000), epcvdata(eprangemed, eff_complex, 0.0, reps, 24000), epcvdata(eprangeslow, eff_complex, 0.0, reps, 24000))
    dataCV_complex = DataFrame(eprange = log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), CV = CV_complex)
    CSV.write(joinpath(abpath(), "data/CVresult_eff071.csv"), dataCV_complex)
end

