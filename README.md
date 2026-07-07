# CH Competition Mouse Models

Julia codebase for simulating clonal hematopoiesis (CH) competition in mice. Models variant allele frequency evolution under selection using a stochastic (Moran process) approach.

## Quick start

```julia
include("src/competitiveSelection.jl")
using .CompetitiveSelection

params = (;
    sType="gamma", η=0.1, σ=0.02, q=0.5,
    N=1000, T=100, μ=2, τ=1.0,
    growthModel="fixed size",
)
solEns, simArgs = evolvePopSim(params; runs=10)
```

For simulation-only work (minimal dependency load):

```julia
include("src/competitionSDE.jl")
using .CompetitionSDE
```

## Setup

Dependencies are managed via `Project.toml`. After cloning:

```bash
julia -e 'import Pkg; Pkg.instantiate()'
```

Run the smoke test to verify:

```bash
julia test/runtests.jl
```

## Parameters

Simulation parameters are passed as a `NamedTuple` to `evolvePopSim`.
Required fields are marked with **req.** below; all others have defaults
set by [`complete`](@ref).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sType` | `String` | **req.** | Fitness distribution type (see below) |
| `growthModel` | `String` | `"fixed size"` | `"fixed size"` or `"unconstrained"` population |
| `T` | `Real` | **req.** | Simulation end time |
| `η` | `Float64` | **req.** | Mean of fitness distribution |
| `σ` | `Float64` | `0.0` | Std dev of fitness distribution (for gamma/gaussian) |
| `N` | `Int` | **req.** | Population size |
| `μ` | `Real` | **req.** | Mutation rate (expected arrivals per generation) |
| `τ` | `Real` | **req.** | Wild-type generation time |
| `α` | `Float64` | `1/τ` | Wild-type division rate (inferred from τ if not given) |
| `q` | `Float64` | `0.0` | Double-hit mutation probability |
| `tMature` | `Real` | `0.0` | Time to maturity (for prenatal growth phase) |
| `sMax` | `Real` | `Inf` | Ceiling on per-variant selection coefficient draws |

## Selection models

The `sType` string selects how per-variant selection coefficients are
drawn from a distribution with mean `η`.

| `sType` value | Struct | Distribution |
|---|---|---|
| `"fixed"` | `FixedSelectionModel(η, q)` | All variants get coefficient `η` |
| `"free"` | `FreeFixedModel(η, q)` | Identical to Fixed; used in fitting contexts |
| `"gaussian"` | `GaussianSelectionModel(η, σ, q)` | `Normal(η, σ)` |
| `"exponential"` | `ExponentialSelectionModel(η, q)` | `Exponential(η)` |
| `"gamma"` | `GammaSelectionModel(η, σ, q)` | `Gamma(η²/σ², σ²/η)` |

Each struct stores a `q` field — the probability that a new variant
acquires a second hit on its parent's fitness (`0` disables).

## Growth models

| Growth model | Struct | Behaviour |
|---|---|---|
| `"fixed size"` | `FixedSizeGrowthModel` | Constant population size. Drift includes competition term `sᵢ - s̄`. Diffusion scales with `√(α(2 + sᵢ + s̄)x/N)`. |
| `"unconstrained"` | `UnconstrainedGrowthModel` | Exponential growth. Drift is `α·sᵢ·x` (no competition). Diffusion scales with `√(α·x·(2 + sᵢ))`. |

## Usage patterns

### Basic simulation

```julia
params = (sType="gamma", η=0.1, σ=0.02, q=0.5, N=1000, T=100, μ=2, τ=1.0)
solEns, simArgs = evolvePopSim(params; runs=50)
```

### With prenatal growth phase

Set `tMature` to the age at maturity and pass `growthPhase=true`.
Variants arising before `tMature` evolve under the unconstrained
(exponential) model; after `tMature` the specified growth model applies.

```julia
params = (..., tMature=0.5)
solEns, simArgs = evolvePopSim(params; runs=20, growthPhase=true)
```

### Tracked variants

Pass birth times of specific variants to track their trajectories.
Provide as 1-element tuples `(t0,)` (just birth time, draws a random
selection coefficient) or 2-element tuples `(t0, s)` (birth time with
a user-specified selection coefficient).

Multiple tracked variants can share the same birth time — each is
assigned a unique index. The `_trackerID` column in the returned
DataFrame gives the position of each tracked variant in the `t₀_vid`,
`s_vid`, `x₀_vid`, etc. columns.

```julia
# Single variant with a specified fitness
solEns, simArgs = evolvePopSim(params; runs=10, _trackerVariant=[(5.0, 0.2)])
# Multiple variants at the same time
solEns, simArgs = evolvePopSim(params; runs=10, _trackerVariant=fill((0.001, 0.4), 5))
# Birth time only (fitness drawn from distribution)
solEns, simArgs = evolvePopSim(params; runs=10, _trackerVariant=[(0.5,), (10.0,)])

# Access tracked variant trajectories
ids = simArgs._trackerID[1]
solEns.u[1](t)[ids]  # frequencies at time t
```

### Changing the solver

Available algorithms (see [`ALGS`](@ref)): `:LambaEM` (default),
`:SOSRI2`, `:SOSRA2`, `:SRA3`, `:ImplicitEulerHeun`, `:SKenCarp`.

```julia
solEns, simArgs = evolvePopSim(params; runs=10, algorithm=:SOSRI2)
```

### Disabling diffusion

```julia
solEns, simArgs = evolvePopSim(params; runs=10, noDiffusion=true)
```

## Output

`evolvePopSim` returns a tuple `(solEns, simArgs)`:

- **`solEns::EnsembleSolution`** — one trajectory per run. Access the
  `i`-th run with `solEns.u[i]`, or interpolate at time `t` with
  `solEns[i](t)`.

- **`simArgs::DataFrame`** — one row per run with columns:
  `k` (variant count), `t₀_vid` (arrival times), `_trackerID`,
  `s_vid` (selection coefficients), `init_vid`, `x₀_vid`, `parentId_vid`.

## Module reference

Three modules in a two-level hierarchy:

| Module | File | Purpose |
|---|---|---|
| `CompetitiveSelection` | `competitiveSelection.jl` | Top-level convenience, re-exports both submodules |
| ↳ `CompetitionSDE` | `competitionSDE.jl` | **Simulation core**: types, SDE engine, `evolvePopSim` |
| ↳ `Scientist` | `Scientist.jl` | Analysis (binning, trajectory sampling, ABC) |

`CompetitionSDE` is self-contained and can be loaded standalone.
`Scientist` depends on `CompetitionSDE` via `using ..CompetitionSDE`.

### Source files

| File | Included by | Contents |
|---|---|---|
| `competitionSDE.jl` | — | `CompetitionSDE` module: types, `complete()`, `ALGS` |
| `simEvolver.jl` | `competitionSDE.jl` | `evolvePopSim`, `prepareSims`, drift/diffusion, selection dispatch |
| `Scientist.jl` | `competitiveSelection.jl` | `Scientist` module: binning, sampling, ABC |

### Single exported function

| Export | Description |
|---|---|
| `evolvePopSim` | Main simulation entry point (see [Usage patterns](#usage-patterns)) |

All types and internal functions are accessed through the submodule—e.g.
`CompetitionSDE.SimArgs`, `CompetitionSDE.prepareSims`—for advanced use.

## Testing

```bash
julia test/runtests.jl
```

Runs smoke tests covering basic simulations, the prenatal growth phase,
tracked variants (1- and 2-element tuples, same-time and distinct-time),
and `_trackerID` correctness.

## Attribution

Forked from [competitive-selection](https://github.com/natevmp/competitive-selection).
