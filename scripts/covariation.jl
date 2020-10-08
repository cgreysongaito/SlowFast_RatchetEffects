include("packages.jl")

function covartest(sd1, sd2)
    x = rand(Normal(0.0, sd1), 100)
    y = rand(Normal(0.0, sd2), 100)
    xme = mean(x)
    yme = mean(y)
    sumvec = fill(0.0, 100, 1)
    for i in 1:100
        sumvec[i] = (x[i] - xme) * (y[i] - yme)
    end

    return sum(sumvec)/99
end

covartest(0.01, 0.04)
covartest(0.01, 0.08)
covartest(0.01, 0.1)
covartest(0.01, 0.2)
covartest(0.01, 0.3)

#but above is not a wiener process so need to program wiener process
# https://www.youtube.com/watch?v=ld0rxwAJpkM
