using DrWatson
@quickactivate "SNSensitivityEstimate"

using CairoMakie

# --- Oscillation parameters (approx PDG values) ---
Δm21 = 7.4e-5
Δm31_NO = 2.5e-3
Δm31_IO = -2.5e-3

s12 = 0.307
s13 = 0.022

c12 = 1 - s12
c13 = 1 - s13

# --- Lightest mass grid ---
m_light = 10 .^ range(-4, 0, length=400)

# --- Mass definitions ---
masses_NO(m1) = (m1,
                 sqrt(m1^2 + Δm21),
                 sqrt(m1^2 + Δm31_NO))

masses_IO(m3) = (sqrt(m3^2 + abs(Δm31_IO)),
                 sqrt(m3^2 + abs(Δm31_IO) + Δm21),
                 m3)

# --- Effective mass range ---
function mee_range(m1, m2, m3)
    a = c12*c13*m1
    b = s12*c13*m2
    c = s13*m3

    mmax = a + b + c
    mmin = max(a - b - c, 0.0)

    return mmin, mmax
end

# --- Compute bands ---
mee_min_NO = Float64[]
mee_max_NO = Float64[]
mee_min_IO = Float64[]
mee_max_IO = Float64[]

for m in m_light
    m1, m2, m3 = masses_NO(m)
    mn, mx = mee_range(m1, m2, m3)
    push!(mee_min_NO, mn)
    push!(mee_max_NO, mx)

    m1, m2, m3 = masses_IO(m)
    mn, mx = mee_range(m1, m2, m3)
    push!(mee_min_IO, mn)
    push!(mee_max_IO, mx)
end

# --- Plot ---
fig = Figure(size = (700, 500))
ax = Axis(fig[1, 1],
    xlabel = "Lightest neutrino mass [eV]",
    ylabel = "⟨mββ⟩ [eV]",
    xscale = log10,
    yscale = log10
)

# --- Bands ---
band!(ax, m_light, mee_min_NO, mee_max_NO, transparency=true, label="Normal ordering")
band!(ax, m_light, mee_min_IO, mee_max_IO, transparency=true, label="Inverted ordering")

# --- Experimental limits ---
lines!(ax, m_light, fill(0.2, length(m_light)), linestyle = :dash, label = "KATRIN (approx)")
lines!(ax, m_light, fill(0.05, length(m_light)), linestyle = :dot, label = "Cosmology (approx)")

axislegend(ax, position = :rb)

fig