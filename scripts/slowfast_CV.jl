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
end #This function calculates the CV of a time series of the Rosenzweig-MacArthur model with stochasticity.

function CVdata(ep, eff, r, reps, tsend)
    cvdata = Union{Missing, Float64}[]
    for i in 1:reps
        push!(cvdata, CVcalc(ep, eff, r, i, tsend))
    end
    return mean(skipmissing(cvdata))
end #This function calculates the mean of CV from many simulations of the Rosenzweig-MacArthur model with stochasticity.

function epcvdata(eprange, eff, r, reps, tsend)
    epdata = zeros(length(eprange))
    @threads for i in eachindex(eprange)
        @inbounds epdata[i] = CVdata(eprange[i], eff, r, reps, tsend)
    end
    return reverse(epdata)
end #This function varies ϵ and calculates the mean CV for each value of ϵ.

eprangeslow = 0.001:0.001:0.05;
eprangemed = 0.05:0.01:0.2;
eprangefast = 0.2:0.05:1.0;
eff_real = 0.5;
eff_complex = 0.71;
reps = 5;

#Non-excitable
let 
    CV_real = vcat(epcvdata(eprangefast, eff_real, 0.0, reps, 24000), epcvdata(eprangemed, eff_real, 0.0, reps, 24000), epcvdata(eprangeslow, eff_real, 0.0, reps, 24000))
    dataCV_real = DataFrame(eprange = log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), CV = CV_real)
    CSV.write(joinpath(abpath(), "data/CVresult_eff05.csv"), dataCV_real)
end

#Highly excitable
let
    CV_complex = vcat(epcvdata(eprangefast, eff_complex, 0.0, reps, 24000), epcvdata(eprangemed, eff_complex, 0.0, reps, 24000), epcvdata(eprangeslow, eff_complex, 0.0, reps, 24000))
    dataCV_complex = DataFrame(eprange = log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), CV = CV_complex)
    CSV.write(joinpath(abpath(), "data/CVresult_eff071.csv"), dataCV_complex)
end

#Moderately excitable (in Supporting Information)
eff_smallcomplex = 0.6;
let
    CV_smallcomplex = vcat(epcvdata(eprangefast, eff_smallcomplex, 0.0, reps, 24000), epcvdata(eprangemed, eff_smallcomplex, 0.0, reps, 24000), epcvdata(eprangeslow, eff_smallcomplex, 0.0, reps, 24000))
    dataCV_smallcomplex = DataFrame(eprange = log10.(vcat(1 ./ reverse(eprangefast), 1 ./ reverse(eprangemed), 1 ./ reverse(eprangeslow))), CV = CV_smallcomplex)
    CSV.write(joinpath(abpath(), "data/CVresult_eff06.csv"), dataCV_smallcomplex)
end

