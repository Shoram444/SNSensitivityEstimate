using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using Turing, Distributions, StatsPlots, ReverseDiff, SensitivityModule

mu, s, n = 5, 1, 1000

d = rand(Normal(mu,s), n);

@model function normal_fit(d)
    mu ~ truncated(Normal(5, 1), 0, Inf)
    sigma ~ truncated(Normal(1, 0.01), 0, Inf)
    for i in eachindex(d)
        d[i] ~ Normal(mu,sigma)
    end
    (; mu, sigma)
end
ch = sample(normal_fit(d), NUTS(1000, 0.65), MCMCSerial(), 1000, 1)

plot(ch)
histogram(d)


fit_f(x, params, bw,n) = bw*n*pdf(Normal(params[1], params[2]),x)

bw=0.1
bins = 0:bw:10

fit_f1(x, params) = fit_f(x, params, bw,n)
means, lowers, uppers = get_mean_conf(bins, ch, fit_f1, 0.1)

histogram(d, bins=bins, label="data")
plot!(bins, means, label = "mean", lw =1, ls=:dash, c=5)
plot!(bins, means, ribbon=(lowers,uppers),c=5, label = "90% credible interval", lw =0, fa = 0.7)
ylims!(0, 60)
savefig("plots/test_turing_CI.png")

############################################
using BAT, DensityInterface

prior = distprod(
    mu = truncated(Normal(5, 1), 0, Inf),
    sigma = truncated(Normal(1, 0.01), 0, Inf)
)

l( p::NamedTuple{(:mu, :sigma)}, x::Real) = pdf(Normal(p.mu, p.sigma),x)

struct BatModel
    data
end

@inline DensityInterface.DensityKind(::BatModel) = IsDensity()

function DensityInterface.logdensityof(likelihood::BatModel, p)
    ll = [log(l(p, x)) for x in likelihood.data]
    return sum(ll)
end

posterior = PosteriorMeasure(BatModel(d), prior)
samples = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result

fit_f(p::@NamedTuple{mu::Float64, sigma::Float64}, x; bw=bw, n=n) = bw*n*pdf(Normal(p.mu, p.sigma),x)
histogram(d, bins=bins, label="data")
plot!(bins, fit_f, samples, thickness_scaling =1, intervals=[0.9],lw=1, mean=false, size=(800,600))
savefig("plots/test_BAT_CI.png")

