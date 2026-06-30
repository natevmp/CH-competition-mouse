# AGENTS.md — CH Competition Mouse Models

Julia research codebase for simulating clonal hematopoiesis (CH) competition in mice.

## Entry point

- module `CompetitiveSelection` in `src/competitiveSelection.jl`
- Load: `include("src/competitiveSelection.jl")` then `using .CompetitiveSelection`
- Submodules are auto-included; no separate `using` needed.

## Module layout

| File | Module | Purpose |
|---|---|---|
| `competitiveSelection.jl` | `CompetitiveSelection` | Main entrypoint, growth/selection model types |
| `dataStructuring.jl` | `DataStructuring` | Data loading, binning, ML fitting |
| `analysisTools.jl` | `AnalysisTools` | Logistic/fitness growth functions |
| `simEvolver.jl` | (included) | Stochastic Moran/birth-death simulation |
| `pdeEvolver.jl` | (included) | PDE-based VAF spectrum evolution |
| `scientist.jl` | (included) | Analysis: dist binning, ABC distance |
| `fitModels.jl` | (included) | Model fit structs (Logistic, Constant) |

## Setup

No `Project.toml`. Install dependencies via Julia Pkg:

```julia
using Pkg
Pkg.add(["DifferentialEquations", "DataFrames", "CSV", "LsqFit", "StatsBase",
         "Distributions", "Optim", "Optimization", "Ipopt", "JLD2",
         "ElasticArrays", "Dierckx", "QuadGK", "Interpolations", "Parameters",
         "Distances", "OptimizationBBO", "OptimizationMOI", "ProgressMeter",
         "CairoMakie"])
```

## Testing / verification

There are **no tests or CI**. Verify changes by:
1. Loading the module: `include("src/competitiveSelection.jl")`
2. Running a small simulation or calling exported functions interactively.

## Important notes

- All source files are in `src/` — no package structure, no `Project.toml`.
- Exports are declared at module tops — check `export` lines before adding new public symbols.
- Follow existing Julia style: 4-space indent, docstrings above functions, spaces around operators.
