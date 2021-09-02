# Sensitivity of discrete symmetry metrics: implications for metric choice

This repository contains the source code, Jupyter notebooks, and data used to produce the results for the above titled project.

## Instructions

To run this analysis on your computer, both Julia and Jupyter Notebook must be available. A version of Julia appropriate for your OS can be downloaded from the [Julia language website](https://julialang.org/downloads/), and Jupyter can be installed from within Julia (in the REPL) with

```julia
] add IJulia
```

Alternate instructions for installing Jupyter can be found on the [IJulia github](https://github.com/JuliaLang/IJulia.jl) or the [Jupyter homepage](https://jupyter.org/install) (advanced).

From within the main repository directory, start Julia and then start Jupyter in the Julia REPL

```julia
using IJulia
notebook(;dir=pwd())
```

or if using a system Jupyter installation, start Jupyter from your favorite available shell (e.g. Powershell on Windows, bash on any *nix variant, etc.).

The [`Swing asymmetry`](/notebooks/Swing%20asymmetry.ipynb) notebook is written such that running the entire notebook will reproduce the results reported in the paper. The [`Monte Carlo power simulation`](/notebooks/Monte%20Carlo%20power%20simulation.ipynb) notebook enables reproduction of the simulation (~14 hours using 20 threads with the given parameters) and/or analysing the simulation results and generating figures.
Simulation results are stored in `notebooks/powersimulation-results.h5`. See the power simulation notebook for more details about the storage of the simulation results.
