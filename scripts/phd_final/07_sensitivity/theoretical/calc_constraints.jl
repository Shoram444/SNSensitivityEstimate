using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using ColorSchemes, SNSensitivityEstimate, CairoMakie, DataFramesMeta, LaTeXStrings,StatsBase, FHist, Distributions, PrettyTables

#### MAJORON
proc = [
    "M1 - freq",
    "M2 - freq",
    "M1 - Bay",
    "M2 - Bay"
]

t = [ 1.7e23,2.71e22, 2.34e23, 9.96e22]

G = [3.85e-16, 3.824e-18, 3.85e-16, 3.824e-18]

NME_low = [ 2.18, 3.233e-3, 2.18, 3.233e-3 ]
NME_high = [5.40, 3.233e-3, 5.40, 3.233e-3]
m = [1,2,1,2]
n = [1,3,1,3]

get_gee(t, G, NME,m) = inv(t * G * NME^2)^inv(2*m)

gee_low =  [get_gee(t[i], G[i], NME_low[i], m[i]) for i in 1:length(proc)]
gee_high = [get_gee(t[i], G[i], NME_high[i], m[i]) for i in 1:length(proc)]

d = DataFrame(proc = proc, m = m, n = n, t = t, G = G, NME_low = NME_low, NME_high = NME_high,  gee_low = gee_low, gee_high = gee_high)

pretty_table(d, header = ["Process - method", "m", "n", "sensitivity (yr)", "G (yr-1)", "NME (low)", "NME (high)", "gee (low)", "gee (high)"])

# Sterile

p = ["0.5MeV", "1.5MeV"]
Gsm = [ 1.723e-18, 1.723e-18 ]
GN = [ 1.14e-18, 9.3e-20 ]
tSM = [8.69e19, 8.69e19]
tN = [1.5e22, 5.6e21]

function get_sin2(tSM, tN, Gsm, GN ) 

    gamma_SM = log(2) / tSM
    gamma_N = log(2) / tN

    g = (Gsm/(2*GN)) 
    t = gamma_N / gamma_SM
    return t*g 
end

sin2 = [get_sin2(tSM[i], tN[i], Gsm[i], GN[i]) for i in 1:length(p)]

d = DataFrame(p = p, Gsm = Gsm, GN = GN, tSM = tSM, tN = tN, sin2 = sin2)
pretty_table(d, header = ["Sterile neutrino mass", "Gsm (yr-1)", "GN (yr-1)", "tSM (yr)", "tN (yr)", "sin2"])


