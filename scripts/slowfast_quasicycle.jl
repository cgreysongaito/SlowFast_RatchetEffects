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
end


