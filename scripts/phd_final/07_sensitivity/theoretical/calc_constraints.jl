using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using ColorSchemes, SNSensitivityEstimate, CairoMakie, DataFramesMeta, LaTeXStrings, StatsBase, FHist, Distributions, PrettyTables
using PrettyTables

#### MAJORON CALCULATIONS
proc = [
    "M1 - freq",
    "M2 - freq",
    "M1 - Bay",
    "M2 - Bay"
]

t = [ 1.7e23, 2.71e22, 2.34e23, 9.96e22]
G = [3.85e-16, 3.824e-18, 3.85e-16, 3.824e-18]
NME_low = [ 2.18, 3.233e-3, 2.18, 3.233e-3 ]
NME_high = [5.40, 3.233e-3, 5.40, 3.233e-3]
m = [1, 2, 1, 2]
n = [1, 3, 1, 3]

get_gee(t, G, NME, m) = inv(t * G * NME^2)^inv(2*m)

gee_low =  [get_gee(t[i], G[i], NME_low[i], m[i]) for i in 1:length(proc)]
gee_high = [get_gee(t[i], G[i], NME_high[i], m[i]) for i in 1:length(proc)]

df_majoron = DataFrame(
    process = proc,
    m = m,
    n = n,
    sensitivity_yr = t,
    G = G,
    NME_low = NME_low,
    NME_high = NME_high,
    gee_low = gee_low,
    gee_high = gee_high
)

pretty_table(df_majoron)

#### STERILE NEUTRINO MIXING CALCULATIONS

# Helper functions for sin²θ calculations
function get_sin2_cupid(tSM, tN, Gsm, GN)
    """CUPID-0 PSF convention: sin²θ = (Gsm/2GN) * (γ_N/γ_SM)"""
    gamma_SM = log(2) / tSM
    gamma_N = log(2) / tN
    g = Gsm / (2 * GN)
    t = gamma_N / gamma_SM
    return t * g
end

function get_sin2_ovidiu(tSM, tN, Gsm, GN)
    """Ovidiu PSF convention: sin²θ = (Gsm/GN) * (γ_N/γ_SM)"""
    gamma_SM = log(2) / tSM
    gamma_N = log(2) / tN
    g = Gsm / (2*GN)
    t = gamma_N / gamma_SM
    return t * g
end

# Mass range
masses = collect(0.5:0.1:1.5)
p_new = [string(m, "MeV") for m in masses]

# SuperNEMO with Ovidiu PSFs
Gsm_ovidiu = 1.723e-18 .* ones(length(p_new))
GN_ovidiu = [1.140e-18, 9.694e-19, 8.061e-19, 6.561e-19, 5.228e-19, 4.078e-19, 3.113e-19, 2.323e-19, 1.692e-19, 1.202e-19, 8.308e-20]
tSM_ovidiu = 8.69e19 .* ones(length(p_new))
tN_ovidiu_freq = [
    1.40980292616615e22, 1.3173436106221405e22, 1.235071515554051e22,
    1.132669423516164e22, 1.0318495539088396e22, 9.224944005087773e21,
    8.179269472427468e21, 7.045765410456978e21, 5.969848879270155e21,
    4.963923382759394e21, 3.980077731841228e21
]
tN_ovidiu_bay =[
1.764978742839722e21,
1.1755980273668492e21,
1.685610614314832e21,
8.252354335739507e20,
1.4652394394211852e21,
1.1515019518930116e21,
6.638618511168462e20,
8.809480555806187e20,
7.734810575038939e20,
6.574709137527557e20,
6.211342182093535e20,
]

sin2_ovidiu_freq= [get_sin2_ovidiu(tSM_ovidiu[i], tN_ovidiu_freq[i], Gsm_ovidiu[i], GN_ovidiu[i]) for i in 1:length(p_new)]
sin2_ovidiu_bay = [get_sin2_ovidiu(tSM_ovidiu[i], tN_ovidiu_bay[i], Gsm_ovidiu[i], GN_ovidiu[i]) for i in 1:length(p_new)]

# SuperNEMO with CUPID-0 PSFs
GN_cupid_psf = [10.5, 8.92, 7.41, 6.02, 4.79, 3.73, 2.84, 2.11, 1.54, 1.09, 0.75] .* 1e-19
Gsm_cupid_psf = 15.9e-19 .* ones(length(p_new))
tSM_cupid_psf = tSM_ovidiu  # same SM half-life
tN_cupid_psf_freq = tN_ovidiu_freq     # same N half-life
tN_cupid_psf_bay = tN_ovidiu_bay     # same N half-life
sin2_cupid_psf_freq = [get_sin2_cupid(tSM_cupid_psf[i], tN_cupid_psf_freq[i], Gsm_cupid_psf[i], GN_cupid_psf[i]) for i in 1:length(p_new)]
sin2_cupid_psf_bay = [get_sin2_cupid(tSM_cupid_psf[i], tN_cupid_psf_bay[i], Gsm_cupid_psf[i], GN_cupid_psf[i]) for i in 1:length(p_new)]

# CUPID-0 published values
sin2_cupid_published = [0.0088, 0.0083, 0.0080, 0.0084, 0.0088, 0.0100, 0.0120, 0.0159, 0.0250, 0.0497, 0.1755]

# CUPID-Mo published values
sin2_cupid_mo_published = [0.033, 0.033, 0.035, 0.039, 0.045, 0.047, 0.049, 0.051, 0.053, 0.063, 0.074] 

# Combine all methodologies into single DataFrame
df_sterile = DataFrame(
    mass_MeV = masses,
    mass_label = p_new,
    sin2_SuperNEMO_Ovidiu_freq = sin2_ovidiu_freq,
    sin2_SuperNEMO_Ovidiu_bay = sin2_ovidiu_bay,
    sin2_SuperNEMO_CUPID0PSF_freq = sin2_cupid_psf_freq,
    sin2_SuperNEMO_CUPID0PSF_bay = sin2_cupid_psf_bay,
    sin2_CUPID0_published = sin2_cupid_published,
    sin2_CUPIDMo_published = sin2_cupid_mo_published
)

pretty_table(df_sterile)

#### PLOTTING
let
    colors = ["#003865", "#FFB948", "#52473B", "#9A3A06", "#5B1A18", "#3B0A45", "#1B0033"]

    with_theme(theme_latexfonts()) do
        fig = Figure(size = (900, 600), fontsize = 18)
        ax = Axis(
            fig[1, 1],
            xlabel = L"m_N \ \mathrm{(MeV)}",
            ylabel = L"\sin^2\theta",
            title  = L"\textrm{Sensitivity to sterile neutrino mixing}",
            xticks = (masses, [latexstring(m) for m in masses]),
            yscale = log10,
            yminorticksvisible = true,
            yminorgridvisible  = true,
            yminorticks = IntervalsBetween(9),
            limits = ((0.4, 1.6), (3e-3, 1e1))
        )

        lines!(ax, masses, df_sterile.sin2_SuperNEMO_Ovidiu_freq;
            label     = L"\textrm{SuperNEMO (Ovidiu PSFs)}",
            linewidth = 3.5,
            color     = colors[1],
        )
        CairoMakie.scatter!(ax, masses, df_sterile.sin2_SuperNEMO_Ovidiu_freq;
            label     = L"\textrm{SuperNEMO (Ovidiu PSFs)}",
            color      = colors[1],
            markersize = 15,
        )

        CairoMakie.lines!(ax, masses, df_sterile.sin2_SuperNEMO_Ovidiu_bay;
            label     = L"\textrm{SuperNEMO (Ovidiu PSFs, Bayes)}",
            linewidth = 3.5,
            linestyle = :dash,
            color     = colors[2],
        )
        CairoMakie.scatter!(ax, masses, df_sterile.sin2_SuperNEMO_Ovidiu_bay;
            label     = L"\textrm{SuperNEMO (Ovidiu PSFs, Bayes)}",
            color      = colors[2],
            markersize = 15,
            # marker     = :cross,
        )

        lines!(ax, masses, df_sterile.sin2_SuperNEMO_CUPID0PSF_freq;
            label     = L"\textrm{SuperNEMO (CUPID-0 PSFs)}",
            linewidth = 3.5,
            # linestyle = :dash,
            color     = colors[3],
        )
        CairoMakie.scatter!(ax, masses, df_sterile.sin2_SuperNEMO_CUPID0PSF_freq;
            label     = L"\textrm{SuperNEMO (CUPID-0 PSFs)}",
            color      = colors[3],
            markersize = 15,
            marker     = :rect,
        )

        CairoMakie.lines!(ax, masses, df_sterile.sin2_SuperNEMO_CUPID0PSF_bay;
            label     = L"\textrm{SuperNEMO (CUPID-0 PSFs, Bayes)}",
            linewidth = 3.5,
            linestyle = :dash,
            color     = colors[4],
        )
        CairoMakie.scatter!(ax, masses, df_sterile.sin2_SuperNEMO_CUPID0PSF_bay;
            label     = L"\textrm{SuperNEMO (CUPID-0 PSFs, Bayes)}",
            color      = colors[4],
            markersize = 15,
            marker     = :rect,
        )

        lines!(ax, masses, df_sterile.sin2_CUPID0_published;
            label     = L"\textrm{CUPID-0 (published)}",
            linewidth = 3.5,
            linestyle = :dot,
            color     = colors[5],
        )
        CairoMakie.scatter!(ax, masses, df_sterile.sin2_CUPID0_published;
            label     = L"\textrm{CUPID-0 (published)}",
            color      = colors[5],
            markersize = 15,
            marker     = :diamond,
        )

        CairoMakie.lines!(ax, masses, df_sterile.sin2_CUPIDMo_published;
            label     = L"\textrm{CUPID-Mo (published)}",
            linewidth = 3.5,
            linestyle = :dashdot,
            color     = colors[6],
        )
        CairoMakie.scatter!(ax, masses, df_sterile.sin2_CUPIDMo_published;
            label     = L"\textrm{CUPID-Mo (published)}",
            color      = colors[6],
            markersize = 15,
            marker     = :utriangle,
        )

        axislegend(
            ax, 
            position = :lt, 
            framevisible = true, 
            labelsize = 18,
            merge = true, 
            patchsize = (45, 30), 
            patchlabelgap = 15,
            nbanks =2
            )

        save(scriptsdir("phd_final/07_sensitivity/theoretical/sterile_sensitivity.png"), fig, px_per_unit = 5)
        fig

    end
end


df_psfs = DataFrame(
    mass_MeV = masses,
    mass_label = p_new,
    GN_Ovidiu = GN_ovidiu,
    Gsm_Ovidiu = Gsm_ovidiu,
    GN_CUPID0PSF = GN_cupid_psf,
    Gsm_CUPID0PSF = Gsm_cupid_psf
)

pretty_table(df_psfs)

df_tNs = DataFrame(
    mass_MeV = masses,
    mass_label = p_new,
    tN_Ovidiu_freq = tN_ovidiu_freq,
    tN_Ovidiu_bay = tN_ovidiu_bay,
)
pretty_table(df_tNs)