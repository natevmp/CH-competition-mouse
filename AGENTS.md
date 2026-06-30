# AGENTS.md — CH Competition Mouse Models

Julia research codebase for simulating clonal hematopoiesis (CH) competition in mice.

## Entry point

- module `CompetitiveSelection` in `src/competitiveSelection.jl`
- Load: `include("src/competitiveSelection.jl")` then `using .CompetitiveSelection`

## Module layout

| File | Module/scope | Purpose |
|---|---|---|
| `competitiveSelection.jl` | `CompetitiveSelection` | Main entrypoint, growth/selection model types |
| `simEvolver.jl` | (included) | Stochastic Moran/birth-death simulation |
| `scientist.jl` | (included) | Analysis: binning, sampling, ABC distance |

All functions in `simEvolver.jl` and `scientist.jl` land directly in `CompetitiveSelection` (they are `include()`d, not separate modules).

## Setup

No `Project.toml`. Install dependencies via Julia Pkg:

```julia
using Pkg
Pkg.add(["DifferentialEquations", "DataFrames", "StatsBase", "Distributions",
         "Distances", "Interpolations", "Parameters", "JLD2", "LinearAlgebra"])
```

## Testing / verification

Run the smoke test:
```bash
julia test/runtests.jl
```
This loads the module, runs a short `evolvePopSim` simulation, and asserts basic type/shape invariants.

## Important notes

- All source files are in `src/` — no package structure, no `Project.toml`.
- Single exported function: `evolvePopSim` (declared at `competitiveSelection.jl:15`).
- Follow existing Julia style: 4-space indent, docstrings above functions, spaces around operators.
- Array indexing follows the underscore protocol (see global `AGENTS.md` for the naming convention).
