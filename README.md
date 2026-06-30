# CH Competition Mouse Models

Julia codebase for simulating clonal hematopoiesis (CH) competition in mice. Models variant allele frequency evolution under selection using a stochastic (Moran process) approach, with inference via growth model fitting and ABC rejection sampling.

## Quick start

```julia
include("src/competitiveSelection.jl")
using .CompetitiveSelection
```

## Dependencies

Dependencies are managed via `Project.toml`. After cloning:

```bash
julia -e 'import Pkg; Pkg.instantiate()'
```

## Module structure

| File | Scope | Purpose |
|---|---|---|
| `competitiveSelection.jl` | `CompetitiveSelection` | Entrypoint, growth/selection model types, `SimArgs` |
| `simEvolver.jl` | (included) | Stochastic Moran / birth-death simulation |
| `scientist.jl` | (included) | Analysis: binning, sampling, ABC distance metrics |

## Attribution

Forked from [competitive-selection](https://github.com/natevmp/competitive-selection).
