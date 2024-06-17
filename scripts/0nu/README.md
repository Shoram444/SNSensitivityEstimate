## Data cuts used:
- exactly 2 foil vertices less than 5cm apart
- exactly 2 trigered OM hits
- exactly 2 associated OM hits
- exactly 2 tracks (negative charged in case of B-on; charged undefined in case of B-off)
- ToF cut: Pint > 4% and Pext < 1%
- Sum energy in 0-3500keV


## Used backgrounds/acitivities:
|process           |activitiy [Bq/kg or Bq/m3]|amount [kg of m3]|
|------------------|--------------------------|-----------------|
|2nubb_foil_bulk   | 0.001721                 | 6.07 kg         |
|Bi214_foil_bulk   | 1e-5                     | 6.25 kg         |
|Bi214_wire_furface| 0.00015                  | 15 m3           |
|Tl208_foil_bulk   | 2e-6                     | 6.25 kg         |


## Sensitivity table:
|FWHM |B field| sig_eff |ROI      | bkg | Nbexcl | T12 freq  | mu_U | T12 Bay|
|-----|-------|---------|---------|-----|--------|-----------|------|--------|
| 8%  | on    | 11.73%  |2700-3100| 0.61| 2.92   |3.58e24    |2.31  | 4.41e24|
|12%  | on    | 11.58%  |2700-3200| 1.30| 3.40   |3.04e24    |2.47  | 4.36e24|
|12%  | off   | 14.32%  |2700-3200| 1.65| 3.62   |3.53e24    |2.34  | 5.39e24|
