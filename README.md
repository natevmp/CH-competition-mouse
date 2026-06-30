# CH Competition Mouse Models

Julia codebase for simulating clonal hematopoiesis (CH) competition in mice. Models variant allele frequency evolution under selection using a stochastic (Moran process) approach, with inference via growth model fitting and ABC rejection sampling.

## Quick start

```julia
include("src/competitiveSelection.jl")
using .CompetitiveSelection
```

For simulation-only work (minimal dependencies):

```julia
include("src/competitionSDE.jl")
using .CompetitionSDE
```

## Dependencies

Dependencies are managed via `Project.toml`. After cloning:

```bash
julia -e 'import Pkg; Pkg.instantiate()'
```

## Module structure

Three modules in a two-level hierarchy:

| Module | File | Purpose |
|---|---|---|
| `CompetitiveSelection` | `competitiveSelection.jl` | Top-level convenience — re-exports both submodules |
| ↳ `CompetitionSDE` | `competitionSDE.jl` | Simulation core: types, SDE engine, `evolvePopSim` |
| ↳ `Scientist` | `Scientist.jl` | Analysis: binning, trajectory sampling, ABC distance metrics |

`CompetitionSDE` is self-contained and can be used standalone. `Scientist` depends on `CompetitionSDE` and is loaded through the parent module.

Simulation code lives in `simEvolver.jl` (included by `CompetitionSDE`).

## Attribution

Forked from [competitive-selection](https://github.com/natevmp/competitive-selection).
