# CH Competition Mouse Models

Julia codebase for simulating clonal hematopoiesis (CH) competition in mice. Models variant allele frequency evolution under selection using stochastic (Moran process) and PDE-based approaches, with inference via growth model fitting and ABC rejection sampling.

## Quick start

```julia
include("src/competitiveSelection.jl")
using .CompetitiveSelection
```

## Dependencies

Install packages manually via Julia Pkg:

```julia
using Pkg
Pkg.add(["DifferentialEquations", "DataFrames", "CSV", "LsqFit", "StatsBase",
         "Distributions", "Optim", "Optimization", "Ipopt", "JLD2",
         "ElasticArrays", "Dierckx", "QuadGK", "Interpolations", "Parameters",
         "Distances", "OptimizationBBO", "OptimizationMOI", "ProgressMeter"])
```

## Module structure

| Module | File | Purpose |
|---|---|---|
| `CompetitiveSelection` | `competitiveSelection.jl` | Entrypoint, growth/selection model types |
| `DataStructuring` | `dataStructuring.jl` | Data loading, binning, MLE fitting |
| `AnalysisTools` | `analysisTools.jl` | Logistic/fitness growth functions |
| (included) | `simEvolver.jl` | Stochastic Moran / birth-death simulation |
| (included) | `pdeEvolver.jl` | PDE-based VAF spectrum evolution |
| (included) | `scientist.jl` | Distribution binning, ABC distance metrics |
| (included) | `fitModels.jl` | Model fit structs (Logistic, Constant) |

## Attribution

Forked from [competitive-selection](https://github.com/natevmp/competitive-selection).
