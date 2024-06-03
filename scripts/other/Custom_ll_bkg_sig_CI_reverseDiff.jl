using DrWatson
@quickactivate

push!(LOAD_PATH, srcdir())
using SensitivityModule,Turing, Distributions, StatsPlots, ReverseDiff

mu, s, l, n1, n2 = 5, 1, 0.5, 3, 100

sig = rand(Normal(mu,s), n1);
bkg = rand(Exponential(l), n2);

d= vcat(sig, bkg);

@model function bkg_sig(d)
      
    lambda ~ Uniform(10^-2, 10.)
    mu ~ truncated(Normal(5, 1), 0, Inf)
    sigma1 ~ truncated(Normal(1, 0.01), 0, Inf)
    # Prior distributions for the proportions
    nSignal ~ Uniform(0.0, 100)  # Prior for proportion of Normal data
    nBkg ~ Uniform(0.0, 100)  # Prior for proportion of Normal data
    
    # Likelihood of data
    for x in d
        likelihood_normal = pdf(Normal(mu, sigma1), x)
        likelihood_exponential = pdf(Exponential(lambda), x)
        likelihood_mixture = (nSignal * likelihood_normal + nBkg * likelihood_exponential) 
        likelihood_extended = pdf(Poisson(nSignal + nBkg), length(d))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end
ch = sample(bkg_sig(d), NUTS(0.65), MCMCSerial(), 1_000, 1)

plot(ch)
histogram(d)


# quantiles?
fit_f(x,params, bw) = bw*(params[4]*pdf(Normal(params[2], params[3]),x) + params[5]*pdf(Exponential(params[1]),x))

bw=0.1
bins = 0:bw:10

ff(x,params) = fit_f(x,params, bw)

m, l, u = get_mean_conf(bins, ch, ff, 0.1)

histogram(d, bins=bins, label="data")
plot!(bins, m, label="mean fit", lw=3, c=5, ls = :dash)
plot!(bins, u+m, fill=m-l, label ="90% credible interval", fa = 0.7, lw=1, c=5)

