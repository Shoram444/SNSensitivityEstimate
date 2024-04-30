using Turing, Distributions, StatsPlots, ReverseDiff

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
ch = sample(bkg_sig(d), NUTS(0.65; adtype=AutoReverseDiff()), MCMCSerial(), 1_000, 1)

plot(ch)
histogram(d)


# quantiles?
samples = Array(ch[:,:,1])


lower_qs = mapslices(x -> quantile(x, 0.025), samples; dims=1)
upper_qs = mapslices(x -> quantile(x, 0.975), samples; dims=1)
mid_qs = mapslices(x -> quantile(x, 0.5), samples; dims=1)

fit_f(x,qs, bw) = bw*(qs[4]*pdf(Normal(qs[2], qs[3]),x) + qs[5]*pdf(Exponential(qs[1]),x))

bw=0.1
bins = 0:bw:10

mids = [fit_f(x,mid_qs, bw) for x in collect(bins)];
lows = [fit_f(x,lower_qs, bw) for x in collect(bins)];
ups = [fit_f(x,upper_qs, bw) for x in collect(bins)];

histogram(d, bins=bins, label="data")
plot!(bins, mids, label="mid", lw=1, ls = :dash)
plot!(bins, lows, fill = ups, label ="95% credible interval", fa = 0.7, lw=0)


