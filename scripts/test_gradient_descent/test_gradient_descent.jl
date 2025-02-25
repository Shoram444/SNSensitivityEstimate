import Pkg
Pkg.activate(".")
using BlackBoxOptim, Distributions, Plots

#mock signal, background data
Es, As = rand(Normal(3, 1), 100), rand(Normal(4, 1), 100)
Eb, Ab = rand(Exponential(), 100), rand(Exponential(0.5), 100)

histogram(Es, label="Signal", alpha=0.5)
histogram!(Eb, label="Background", alpha=0.5)

histogram(As, label="Signal", alpha=0.5)
histogram!(Ab, label="Background", alpha=0.5)

signal = (
    angle = As,
    energy = Es
)

background = (
    angle = Ab,
    energy = Eb
)

n = 100

Emin = 0:1:6
Emax = 0:1:6
Amin = 0:1:10
Amax = 0:1:10

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
    return b != 0 ? s/sqrt(b) : 0.0
end

roi1 = [0,1, 0, 1]

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