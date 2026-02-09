# nu0bb 150uBq/m3 radon tag1

tHalf: 4.53e24
gamma = 1.53e-25
signalEff: 0.19
bkg count: 1.65
best ROI:
  phi : (15, 180.0)
  sumE : (2700, 3000)
  r : (0.0, 200.0)
  singleE : (0.0, 2650)
  dy : (0.0, 100)
  dz : (0.0, 100)
  lPint : (0.0, 4.0) -> 0.001
  lPext : (2.0, 50.0) -> 0.01
  trackLength1 : (0.0, 3000.0)
  trackLength2 : (0.0, 3000.0)


tHalf: 4.5913875567922416e24
signalEff: 0.1780742
bkg count: 1.238946815496673
best ROI:
  phi : (15.931636543956474, 178.37044028733078)
  sumE : (2718.5178609744485, 3011.755549622048)
  dy : (0.0010815545018436755, 136.0952588036475)
  dz : (0.16439178774161617, 141.38737589743906)
  lPint : (0.0, 4.45829380685108)
  lPext : (2.461977743955876, 50.0)


# nu0bb radon tag2 2mBq

tHalf: 3.507943685110317e24
signalEff: 0.170782
bkg count: 2.800524051716381
best ROI:
  phi : (34.20293856455696, 177.96183459178883)
  sumE : (2717.3224078520398, 3006.808420273045)
  dy : (0.030408720971108093, 62.92205336545711)
  dz : (0.027592791283349842, 92.84048088366544)
  lPint : (0.0, 3.9387080203444462)
  lPext : (2.1460116284084854, 50.0)


# compare nu0

 Row │ radon_level  method  sens_3yr  eff      bkg      sens_1yr   
     │ Float64      String  Float64   Float64  Float64  Float64    
─────┼─────────────────────────────────────────────────────────────
   1 │        0.15  1D       4.03e24    0.163    1.46   1.76951e24
   2 │        2.0   1D       2.88e24    0.163    4.56   1.37759e24
   3 │        0.15  ND       4.59e24    0.178    1.238  1.97516e24
   4 │        2.0   ND       3.25e24    0.16     2.9    1.51607e24


# nu0bbM1 tag1 radon

CUPID-0
https://journals.aps.org/prd/pdf/10.1103/PhysRevD.107.032006 
n=1 1.2e23 

tHalf: 1.70e23
gamma = 4.10e-24
signalEff: 0.022
bkg count: 35.67
best ROI:
  phi : (20.0, 180.0)
  sumE : (2500, 3000)
  r : (0.0, 200.0)
  singleE : (0.0, 2650.0)
  dy : (0.0, 120)
  dz : (0.0, 120)
  lPint : (0.0, 4.0) -> 0.001
  lPext : (2.0, 50.0) -> 0.01
  trackLength1 : (0.0, 2200)
  trackLength2 : (0.0, 2200)

# nu0bbM2 tag1 radon

CUPID-0
https://journals.aps.org/prd/pdf/10.1103/PhysRevD.107.032006 
n=3 1.4e22 


tHalf: 2.71e22
gamma = 2.56e-23
signalEff: 0.10
bkg count: 43282.28
best ROI:
  phi : (15.0, 175.0)
  sumE : (1400.0, 3000.0)
  r : (0.0, 200.0)
  singleE : (0.0, 2500.0)
  dy : (0.0, 100.0)
  dz : (0.0, 130.0)
  lPint : (0.0, 4.0)
  lPext : (1, 50.0)
  trackLength1 : (0.0, 2500.0)
  trackLength2 : (0.0, 2500.0)

# 2nu RH 050 tag1 radon

tHalf: 1.5238547925234648e22
signalEff: 0.1228176
bkg count: 190356.07687586214
best ROI:
  phi : (10, 180)
  sumE : (400, 2700)
  dy : (0, 135)
  dz : (0, 145)
  lPint : (0.0, 4)
  lPext : (1.3, 50.0)

# RH020

tHalf: 1.5269358234557567e22
signalEff: 0.1240719
bkg count: 193480.85620783508
best ROI:
  phi : (15, 180)
  sumE : (300, 2700)
  dy : (0, 135)
  dz : (0, 145)
  lPint : (0.0, 3)
  lPext : (1.0, 50.0)

# RH040
tHalf: 1.525970035947891e22
signalEff: 0.1246952
bkg count: 195677.17222946716
best ROI:
  phi : (10, 180)
  sumE : (300, 2700)
  dy : (0, 135)
  dz : (0, 145)
  lPint : (0.0, 3)
  lPext : (1.0, 50.0)

# bayes
| **Process**                         | **Radon level**    | **Median sensitivity (yr)** |
|------------------------------------:|-------------------:|----------------------------:|
| $0\nu\beta\beta$                    | 150\,$\mu$Bq/m$^3$ | 4.67198e24                  |
| $0\nu\beta\beta$                    | 2\,$m$Bq/m$^3$     | 3.16956e24                  |
| $0\nu\beta\beta M1$                 | 150\,$\mu$Bq/m$^3$ | 2.33663e23                  |
| $0\nu\beta\beta M2$                 | 150\,$\mu$Bq/m$^3$ | 9.95743e22                  |
| $\nu\_R\nu\_L\beta\beta \, (K=0.5)$ | 150\,$\mu$Bq/m$^3$ | 7.87119e20                  |
