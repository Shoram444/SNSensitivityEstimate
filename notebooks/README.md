This folder contains example notebooks for the sensitivity estimate analysis. There are always (unless I forget..) two pairs of files. One for the actual `Pluto` notebook with `*.jl` and another a `PDF` file for when you just want to see the notebook without running it. 

EXAMPLES
=======
1. The first example in this is the `Example_Xi31` file. The file contains a simple analysis - search for sensitivity toward the *refined* spectrum parametrized by the $\xi_{31}, \xi_{51}$ parameters. Here I generated $10^8$ events of the standard internal background processes and same for the *refined* spectrum, which I consider signal. The example shows how sensitivity is calculated using functionality of `SensitivityModule.jl` for the $E_{sum}$ channel only. 