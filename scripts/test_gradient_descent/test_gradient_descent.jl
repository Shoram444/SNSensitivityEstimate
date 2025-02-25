import Pkg
Pkg.activate(".")
using BlackBoxOptim, Distributions, Plots

#mock signal, background data
N = Int(1e7)
Es, As = rand(Normal(3, 1), N), rand(Normal(7, 3), N)
Eb, Ab = rand(Exponential(7), N), rand(Exponential(0.5), N)

histogram(Es, label="Signal", alpha=0.5)
histogram!(Eb, label="Background", alpha=0.5)

histogram!(As, label="Signal", alpha=0.5)
histogram!(Ab, label="Background", alpha=0.5)

signal = (
    angle = As,
    energy = Es
)

background = (
    angle = Ab,
    energy = Eb
)

n = Int(1e8)

Emin = 0:0.1:20
Emax = 0:0.1:20
Amin = 0:0.01:20
Amax = 0:0.01:20

function get_efficiency(x, Amin, Amax, Emin, Emax, n)
    np = 0
    for i in 1:length(x)
        if Amin <= x.angle[i] < Amax && Emin <= x.energy[i] < Emax
            np += 1
        end
    end
    
    return np/n
end

function get_s_to_b(signal, background, n, rois)
    Amin, Amax, Emin, Emax = rois
    s = get_efficiency(signal, Amin, Amax, Emin, Emax, n)
    b = get_efficiency(background, Amin, Amax, Emin, Emax, n)
    s == 0 && return 0.0 
    return b != 0 ? s/sqrt(b) : 0.0
end

f(x) = -get_s_to_b(signal, background, n, x)

res = bboptimize(
    f; 
    SearchRange = [(Amin[1], Amax[end]), (Emin[1], Emax[end]), (Amin[1], Amax[end]), (Emin[1], Emax[end])], 
    NumDimensions = 4,
    Method=:adaptive_de_rand_1_bin, 
    DX_max_change=0.1,
    MaxSteps = 1000,
    TraceMode = :full,
    TraceInterval = 1,
    SaveTrace = true,
    MaxTime = 0.1
    )

histogram(As, label="Angle Signal", alpha=0.5, xlims=(0.28, 15.2))
histogram!(Ab, label="Angle Bkg", alpha=0.5, xlims=(0.28, 15.2))


histogram(Es, label="Energy Signal", alpha=0.5, xlims=(1.13, 7.25))
histogram!(Eb, label="Energy Bkg", alpha=0.5, xlims=(1.13, 7.25))