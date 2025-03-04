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
|12%  | off   | 14.18%  |2700-3200| 1.64| 3.49   |3.53e24    |2.34  | 5.39e24|


## Bkg contributions by source
#### 12% + Boff + (2700, 3200)keV 
| **process** <br> | **counts** <br> | **activity\_used** <br> (Bq/"amount") |
|------------------------:|------------------------:|----------------------------:|
| bb\_foil\_bulk          | 1.44227                 | 0.001721±3.1e-5             |
| Bi214\_foil\_bulk       | 0.0845939               | 1.0e-5                      |
| Bi214\_radon    | 0.0710698               | 0.00015                     |
| Tl208\_foil\_bulk       | 0.0383909               | 2.0e-6                      |
| total                   | 1.63632                 | --                          |

#### 12% + Boff + (0, 3500)keV 
| **process** <br> | **counts** <br> | **activity\_used** <br> (Bq/"amount") |
|------------------------:|------------------------:|----------------------------:|
| bb\_foil\_bulk          | 85474.5                 | 0.001721±3.1e-5             |
| Bi214\_foil\_bulk       | 8.55966                 | 1.0e-5                      |
| Bi214\_radon            | 10.3197                 | 0.00015                     |
| Tl208\_foil\_bulk       | 1.05707                 | 2.0e-6                      |
| total                   | 85494.4                 | --                          |


#### 12% + Boff + (2300, 3200)keV 
| **process** <br> | **counts** <br> | **activity\_used** <br> (Bq/"amount") |
|------------------------:|------------------------:|----------------------------:|
| bb\_foil\_bulk          | 300.82                  | 0.001721±3.1e-5             |
| Bi214\_foil\_bulk       | 0.760054                | 1.0e-5                      |
| Bi214\_radon            | 0.621646                | 0.00015                     |
| Tl208\_foil\_bulk       | 0.0572866               | 2.0e-6                      |
| total                   | 302.259                 | --                          |

#### 12% + Boff
| **left\_bin\_edge**<br>`Float64` | **bb\_foil\_bulk**<br>`Float64` | **Bi214\_foil\_bulk**<br>`Float64` | **Bi214\_radon**<br>`Float64` | **Tl208\_foil\_bulk**<br>`Float64` | **total**<br>`Float64` |
|---------------------------------:|--------------------------------:|-----------------------------------:|------------------------------:|-----------------------------------:|-----------------------:|
| 0.0                              | 0.0                             | 0.0                                | 0.0                           | 0.0                                | 0.0                    |
| 100.0                            | 0.0                             | 0.0                                | 0.0                           | 0.0                                | 0.0                    |
| 200.0                            | 85.8422                         | 0.00812344                         | 0.0282567                     | 0.00176386                         | 85.8804                |
| 300.0                            | 506.785                         | 0.040845                           | 0.126727                      | 0.0101731                          | 506.963                |
| 400.0                            | 1251.85                         | 0.0940839                          | 0.267154                      | 0.0261484                          | 1252.24                |
| 500.0                            | 2338.79                         | 0.168618                           | 0.362199                      | 0.050502                           | 2339.37                |
| 600.0                            | 3649.92                         | 0.257786                           | 0.460669                      | 0.0762945                          | 3650.71                |
| 700.0                            | 5001.0                          | 0.329664                           | 0.572839                      | 0.0940413                          | 5001.99                |
| 800.0                            | 6248.66                         | 0.369844                           | 0.566846                      | 0.103077                           | 6249.7                 |
| 900.0                            | 7222.89                         | 0.398409                           | 0.581402                      | 0.106152                           | 7223.98                |
| 1000.0                           | 7848.78                         | 0.416516                           | 0.62764                       | 0.102729                           | 7849.92                |
| 1100.0                           | 8075.45                         | 0.426974                           | 0.596815                      | 0.0938673                          | 8076.57                |
| 1200.0                           | 7917.9                          | 0.44362                            | 0.635347                      | 0.0798299                          | 7919.06                |
| 1300.0                           | 7422.64                         | 0.489039                           | 0.601952                      | 0.0661407                          | 7423.8                 |
| 1400.0                           | 6658.47                         | 0.536944                           | 0.637915                      | 0.0537164                          | 6659.7                 |
| 1500.0                           | 5711.03                         | 0.563023                           | 0.562564                      | 0.0408162                          | 5712.2                 |
| 1600.0                           | 4678.7                          | 0.561656                           | 0.597671                      | 0.0291771                          | 4679.89                |
| 1700.0                           | 3643.3                          | 0.551749                           | 0.553145                      | 0.0198859                          | 3644.42                |
| 1800.0                           | 2695.76                         | 0.532769                           | 0.492351                      | 0.0126255                          | 2696.8                 |
| 1900.0                           | 1878.55                         | 0.485774                           | 0.442688                      | 0.00749252                         | 1879.48                |
| 2000.0                           | 1214.63                         | 0.430941                           | 0.378468                      | 0.00434388                         | 1215.44                |
| 2100.0                           | 728.556                         | 0.374134                           | 0.328805                      | 0.00321053                         | 729.262                |
| 2200.0                           | 394.132                         | 0.319092                           | 0.276573                      | 0.00322987                         | 394.731                |
| 2300.0                           | 188.941                         | 0.256287                           | 0.190946                      | 0.00339233                         | 189.391                |
| 2400.0                           | 77.9406                         | 0.195608                           | 0.165259                      | 0.00392999                         | 78.3054                |
| 2500.0                           | 25.8555                         | 0.135783                           | 0.119877                      | 0.00512911                         | 26.1163                |
| 2600.0                           | 6.64097                         | 0.0877825                          | 0.0744948                     | 0.00644426                         | 6.8097                 |
| 2700.0                           | 1.26062                         | 0.0492341                          | 0.0428131                     | 0.00760469                         | 1.36027                |
| 2800.0                           | 0.17438                         | 0.0249397                          | 0.0128439                     | 0.00837831                         | 0.220542               |
| 2900.0                           | 0.00726584                      | 0.00791466                         | 0.0111314                     | 0.00810368                         | 0.0344156              |
| 3000.0                           | 0.0                             | 0.00233454                         | 0.00428131                    | 0.0076395                          | 0.0142554              |
| 3100.0                           | 0.0                             | 0.00017082                         | 0.0                           | 0.00666474                         | 0.00683556             |
| 3200.0                           | 0.0                             | 0.0                                | 0.0                           | 0.00574413                         | 0.00574413             |
| 3300.0                           | 0.0                             | 0.0                                | 0.0                           | 0.00489702                         | 0.00489702             |
| 3400.0                           | 0.0                             | 0.0                                | 0.0                           | 0.00392226                         | 0.00392226             |


## Boff 12% + TKReconstruct + gamma tracking + 2 distinct assoc caloHits + ROI (2700, 3200)keV
| **process**<br>`String` | **counts**<br>`Float64` | **activity\_used**<br>`Any` |
|------------------------:|------------------------:|----------------------------:|
| bb\_foil\_bulk          | 1.17369                 | 0.001721±3.1e-5             |
| Bi214\_foil\_bulk       | 0.0897984               | 1.0e-5                      |
| Bi214\_radon            | 0.0791443               | 0.00015                     |
| Tl208\_foil\_bulk       | 0.0376894               | 2.0e-6                      |
| total                   | 1.38032                 | --                          |

## Boff 12% + TKReconstruct + gamma tracking + 2 distinct assoc caloHits + ROI (0, 3500)keV
| **process**<br>`String` | **counts**<br>`Float64` | **activity\_used**<br>`Any` |
|------------------------:|------------------------:|----------------------------:|
| bb\_foil\_bulk          | 75735.1                 | 0.001721±3.1e-5             |
| Bi214\_foil\_bulk       | 7.91996                 | 1.0e-5                      |
| Bi214\_radon            | 14.0674                 | 0.00015                     |
| Tl208\_foil\_bulk       | 0.913531                | 2.0e-6                      |
| total                   | 75758.0                 | --                          |



## Boff 12% + TKReco + gamma tracking + 2 distinct assoc caloHits + ROI (2700, 3100) + neutron fullshielding no florr
| **process**<br>`String`                             | **counts**<br>`Float64` |
|----------------------------------------------------:|------------------------:|
| bb\_foil\_bulk                                      | 1.29397                 |
| Bi214\_foil\_bulk                                   | 0.0849167               |
| Bi214\_radon                                        | 0.0963836               |
| Tl208\_foil\_bulk                                   | 0.0295838               |
| K40\_foil\_bulk                                     | 0.0                     |
| Pa234m\_foil\_bulk                                  | 0.0                     |
| neutron\_external\nfull\_shielding\_no\_floor\_flux | 0.136586                |
| total                                               | 1.64144                 |

## Boff 12% + TKReco + gamma tracking + 2 distinct assoc caloHits + ROI (2700, 3000) + neutron current setup
| **process**<br>`String`               | **counts**<br>`Float64` |
|--------------------------------------:|------------------------:|
| bb\_foil\_bulk                        | 1.29397                 |
| Bi214\_foil\_bulk                     | 0.0824589               |
| Bi214\_radon                          | 0.0943502               |
| Tl208\_foil\_bulk                     | 0.0224918               |
| K40\_foil\_bulk                       | 0.0                     |
| Pa234m\_foil\_bulk                    | 0.0                     |
| neutron\_external\ncurrent\_shielding | 0.356625                |
| total                                 | 1.84989                 |
