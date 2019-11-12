using StatsBase
using IterTools

function find_delay(eq, par)
    u1 = eq .+ 0.1 # this is just an assummed perturbation, can make this a parameter
    afact = eigen(cmat(eq, par))
    if imag(afact.values[1]) == 0.0
        # not a damped oscillator
        return 0.0
    end
    evecs = afact.vectors
    # We only need one eigenvector because of conjugacy
    a = real.(evecs[:, 2])
    b = imag.(evecs[:, 2])
    c = hcat(a, b) \ u1
    ϕR = atan(c[1] * a[1] + c[2] * b[1], c[2] * a[1] - c[1] * b[1])
    ϕC = atan(c[1] * a[2] + c[2] * b[2], c[2] * a[2] - c[1] * b[2])
    # we just look at fractions of 2π since we are not looking at the phase
    # in the time domain, since abs(phase) ∈ (0, π) than we are looking at 0, 1/2 for the
    # range of this relative lag
    period = 2 * π
    return abs(ϕR - ϕC) / period
end

#
# #Tools for calculating Lags and Phase difference between timeseries
#
function find_peaks(data, tol)
    imax = Vector{Int}()
    maxes = Vector{Float64}()
    for (i, win) in enumerate(partition(data, 3, 1))
        if win[1] < win[2] && win[2] > win[3]
            if abs(win[1] - win[2]) > tol
                push!(maxes, win[2])
                push!(imax, i + 1) # Save the 2nd index (middle point)
            end
        end
    end
    return (imax, maxes)
end

function find_period(data)
    alags = -(length(data) - 1):(length(data) - 1)
    ac = crosscor(data, data, alags)
    locs, peaks = find_peaks(ac)
    return mean(diff(locs) * 0.1)
end

#TODO: what does the `x` mean?, I guess it must be the times
function find_phase(x, y1, y2)
    lags = -(length(x) - 1):(length(x) - 1)
    xcs = crosscor(y1, y2, lags)
    # assume the stepsize is fixed
    step_size = diff(x[1:2])[1]
    return step_size * lags[argmax(xcs)]
end

find_delay(x, y1, y2) = abs(find_phase(x, y1, y2)) / find_period(y1)
