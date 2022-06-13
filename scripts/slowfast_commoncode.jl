#### Common code and functions
#### "Slow organisms exhibit sudden population disappearances in a reddened world" by Greyson-Gaito, Gellner, & McCann.

function abpath()
    replace(@__DIR__, "scripts" => "")
end #function to create a absolute path as a character string

@with_kw mutable struct RozMacPar
    r = 2.0
    k = 3.0
    a = 1.1
    h = 0.8
    e = 0.7
    m = 0.4
    ε = 1.0
end

par_rozmac = RozMacPar()

function roz_mac_II!(du, u, p, t,)
    @unpack r, k, a, h, e, m, ε = p
    R, C = u
    du[1] = r * R * (1 - R / k) - a * R * C / (1 + a * h * R)
    du[2] = ε * ( e * a * R * C / (1 + a * h * R) - m * C )
    return
end #setting up the Rosenzweig-MacArthur model

function roz_mac_II(u, par)
    du = similar(u)
    roz_mac_II!(du, u, par, 0.0)
    return du
end #setting up the Rosenzweig-MacArthur model 

function eq_II(p)
    @unpack r, a, k, h, m, e = p
    eq_II_R = m / (a * (e - h * m))
    eq_II_C = r * (a * h * k * eq_II_R - a * h * eq_II_R^2 + k - eq_II_R) / (a * k)
    return vcat(eq_II_R, eq_II_C)
end #calculating the interior equilibrium of the Rosenzweig-MacArthur model 


randeq(x) = x * ( 1 + rand(Uniform(1e-7, 1e-6)))

jacmat(model, eq, par) = ForwardDiff.jacobian(eq -> model(eq, par), eq) #calculate the jacobian matrix

λ_stability(M) = maximum(real.(eigvals(M))) #calculate the maximum real eigenvalue

function roz_mac_res(R, C, p)
    @unpack r, k, h, a, m = p
    return r * R * (1 - R / k) - (a * R * C / (1 + a * h * R) )
end #returns the rate of change in resource

function roz_mac_con(R, C, eff, ep, p)
    @unpack h, a, m = p
    return ep * ( ( eff * a * R * C ) / (1 + a * h * R) - m * C )
end #returns the rate of change in consumer

function con_iso_roz_mac(p)
    @unpack m, a, h, e = p
    m / (a * (e - h * m))
end #calculate the consumer isocline

function res_iso_roz_mac(R, p)
    @unpack a, k, r, h = p
    r * (a * h * k * R - a * h * R^2 + k - R) / (a * k)
end #function to help calculate the resource isocline

function iso_plot(resrange, par)
    data = [res_iso_roz_mac(R, par) for R in resrange]
    plot(collect(resrange), data)
    return vlines(con_iso_roz_mac(par), 0, maximum(data)+(0.1* maximum(data)), colors="orange")
end #plots the resource and consumer isoclines for the Rosenzweig-MacArthur consumer-resource model

function roz_mac_plot(ep, eff, width, dens)
    resconrange = range(0, stop = 1.7, length=100)
    U = [roz_mac_res(R, C, par_rozmac) for C in resconrange, R in resconrange]
    V = [roz_mac_con(R, C, eff, ep, par_rozmac) for C in resconrange, R in resconrange]
    speed = sqrt.(U.^2 .+ V.^2)
    lw = width .* speed ./ maximum(speed) # Line Widths
    streamplot(collect(resconrange), collect(resconrange), U, V, density = dens, color = "k", linewidth = lw)
    return iso_plot(resconrange, RozMacPar(e = eff))
end #plots the streamplot for the Rosenzweig-MacArthur consumer-resource model with the isoclines

function noise_creation(r, len)
    white = rand(Normal(0.0, 0.01), Int64(len))
    intnoise = [white[1]]
    for i in 2:Int64(len)
        intnoise = append!(intnoise, r * intnoise[i-1] + white[i] )
    end
    c = std(white)/std(intnoise)
    meanintnoise = mean(intnoise)
    scalednoise = zeros(Int64(len))
    for i in 1:Int64(len)
        scalednoise[i] = c * (intnoise[i] - meanintnoise)
    end
    return scalednoise
end #produces noise with a certain autocorrelation. variance of the noise is scaled using method in Wichmann et al. 2005

function RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = RozMacPar()
    par.ε = ep
    par.e = eff
    noise = noise_creation(r, tsend / freq)
    count = 1
    u0 = [eq_II(par)[1], eq_II(par)[2] + noise[1]]
    tspan = (0.0, tsend)

    function pert_cb2(integrator)
        count += 1
        if isapprox(integrator.u[2], 0.00000000; atol = 1e-8)
            integrator.u[2] = 0.00000000
        else
            integrator.u[2] = maximum([integrator.u[2] + noise[count], 0.0])
        end
    end

    cb = PeriodicCallback(pert_cb2, freq, initial_affect = false) #as of may 29th - initial_affect does not actually do the affect on the first time point
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return solend = sol(tvals)
end #function to numerically solve the Rosenzweig-MacArthur consumer-resource model with added noise to the consumer


function pert_timeseries_plot(ep, eff, freq, r, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    return plot(sol.t, sol.u)
end #plots the time series of the stochastically perturbed Rosenzweig-MacArthur consumer-resource model

function pert_consumer_timeseries_plot(ep, eff, freq, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    plot(sol.t, sol[2, :])
    return ylabel("Consumer biomass")
end #plots the consumer time series of the stochastically perturbed Rosenzweig-MacArthur consumer-resource model

function pert_phase_plot(ep, eff, freq, r, seed, tsend, tvals)
    sol = RozMac_pert(ep, eff, freq, r, seed, tsend, tvals)
    plot(sol[1, :], sol[2, :])
    xlabel("Resource")
    return ylabel("Consumer")
end #plots the time series of the stochastically perturbed Rosenzweig-MacArthur consumer-resource model

function prep_data(data, range)
    canard_prop = zeros(length(data))
    canard_plus_axial_prop = zeros(length(data))
    for i in 1:length(data)
        canard_prop[i] = data[i][1]
        canard_plus_axial_prop[i] = data[i][1] + data[i][2]
    end
    return DataFrame(xrange = range, canard = canard_prop, canard_plus_axial = canard_plus_axial_prop)
end #helper function to take data from canard finder simulations and convert to DataFrame

function transcrit_rozmac(p)
    @unpack a, k, h, m = p
    return h * m + m / (a * k)
end #calculate the efficiency parameter value where the transcritical bifurcation occurs

function hopf_rozmac(p)
    @unpack k, a, h, r, m= p
    R = k / 2 - 1 / (2 * a * h)
    return m / (R * a) + ( h * m)
end #calculate the efficiency parameter value where the Hopf bifurcation occurs

function converteff_prop_RCdivide(eff)
    hopf = hopf_rozmac(RozMacPar(e = eff))
    transcrit = transcrit_rozmac(RozMacPar(e = eff))
    hopf_minus_eff = hopf - eff
    propeff = hopf_minus_eff / (hopf-transcrit)
    finalprop = 1 .- propeff
    return finalprop
end #Function to convert efficiency value to effective proportion real value

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
end #Function to calculate the proportion of real in efficiency "parameter space" (real/complex divide between the transcritical and Hopf bifurcation) for a single ε value. 

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
end # Function to calculate the proportion of real in efficiency "parameter space" for multiple values of ε.


function quasicycle_data(ep, reps)
    lrange = 0:1:40
    ACFdata = Vector{Vector{Float64}}(undef,reps)
    avprep = zeros(reps)
    avfinal = zeros(length(lrange))
    @threads for i in 1:reps
        sol = RozMac_pert(ep, 0.6, 1.0, 0.0, i, 10000.0, 9000.0:1.0:10000.0)
        @inbounds ACFdata[i] = autocor(sol[2, 1:end], collect(lrange))
    end
    for j in 1:length(lrange)
        for l in 1:reps
        avprep[l] = ACFdata[l][j]
        end
        avfinal[j] = mean(avprep)
    end
    return avfinal
end #calculates the average ACF value for each lag point for a set of simulations of the Rosenzweig-MacArthur model with added stochasticity (where ε is 1.0)