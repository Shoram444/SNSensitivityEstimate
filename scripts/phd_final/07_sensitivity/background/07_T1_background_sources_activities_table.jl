using DrWatson
@quickactivate "SNSensitivityEstimate"

using SNSensitivityEstimate, DataFramesMeta, LaTeXStrings, PrettyTables

SNparams = SNSensitivityEstimate.SNparams
Na = SNparams["Nₐ"]
W = SNparams["W"]
tHalf = SNparams["SeThalf2nu"] # seconds
tHalf / 3600/24/365 # years
a_2nu = halfLife_to_activity(Na, W, tHalf) # Bq/kg

h = [
    L"2\nu\beta\beta",
    L"{}^{214}\textrm{Bi}", 
    L"{}^{208}\textrm{Tl}", 
    L"{}^{40}\textrm{K}",
    L"{}^{234m}\textrm{Pa}", 
    L"\textrm{Radon gas}",
    L"\textrm{external neutrons}",
    L"\textrm{external \gamma\textrm{-rays}}"
]

p = [
    "bb2nu", 
    "Bi214", 
    "Tl208", 
    "K40", 
    "Pa234m", 
    "Radon", 
    "Neutrons", 
    "Gamma"
]

a = [
    L"%$(a_2nu*1000) \textrm{mBq/kg}",
    L"10 \textrm{\mu Bq/kg}",
    L"2 \textrm{\mu Bq/kg}",
    L"58.7 \pm 0.1 \textrm{mBq/kg}",
    L"17.3 \pm 0.1 \textrm{mBq/kg}",
    L"150 \textrm{\mu Bq/m}^3",
    L"--",
    L"--"
]

df = DataFrame(Symbol = h, Activity = a)
pretty_table(df, backend = Val(:latex), alignment = :l)