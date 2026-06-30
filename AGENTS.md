# AGENTS.md — CH Competition Mouse Models

Julia research codebase for simulating clonal hematopoiesis (CH) competition in mice.

## Entry point

Load everything:
```julia
include("src/competitiveSelection.jl")
using .CompetitiveSelection
```

Simulation-only (standalone):
```julia
include("src/competitionSDE.jl")
using .CompetitionSDE
```

## Module layout

Three modules in a two-level hierarchy:

| Module | File | Purpose |
|---|---|---|
| `CompetitiveSelection` | `competitiveSelection.jl` | Top-level convenience, re-exports both submodules |
| ↳ `CompetitionSDE` | `competitionSDE.jl` | Simulation core: types, SDE engine, `evolvePopSim` |
| ↳ `Scientist` | `Scientist.jl` | Analysis: binning, trajectory sampling, ABC |

`CompetitionSDE` is self-contained (standalone). `Scientist` depends on `CompetitionSDE` via `using ..CompetitionSDE`.

## Source files

| File | Included by | Purpose |
|---|---|---|
| `simEvolver.jl` | `competitionSDE.jl` | Stochastic Moran/birth-death simulation functions |
| `competitionSDE.jl` | `competitiveSelection.jl` | CompetitionSDE module definition + types |
| `Scientist.jl` | `competitiveSelection.jl` | Scientist module definition with analysis code |

## Setup

Project environment managed via `Project.toml`. After cloning, instantiate:
```bash
julia -e 'import Pkg; Pkg.instantiate()'
```

## Testing / verification

Run the smoke test:
```bash
julia test/runtests.jl
```
This loads the module, runs a short `evolvePopSim` simulation, and asserts basic type/shape invariants.

## Important notes

- All source files are in `src/` — not a formal Julia package (no UUID in `Project.toml`), but project dependencies are managed via `Project.toml`.
- Key exported function: `evolvePopSim` (from `CompetitionSDE`).
- Follow existing Julia style: 4-space indent, docstrings above functions, spaces around operators.
- Array indexing follows the underscore protocol (see global `AGENTS.md` for the naming convention).
