module CompetitionSDE

using DataFrames
using DifferentialEquations
using StochasticDiffEq
using Distributions
using Random

"""
    GrowthModel

Abstract type for population growth models. Subtypes determine how
the drift and diffusion terms are computed in the SDE.
"""
abstract type GrowthModel end

"""
    SelectionModel

Abstract type for selection coefficient distributions.
Subtypes define how per-variant fitness values are drawn.
"""
abstract type SelectionModel end

"""
    GaussianSelectionModel(η, σ, q)

Selection coefficients are drawn from `Normal(η, σ)`.
"""
struct GaussianSelectionModel <: SelectionModel
    η::Float64
    σ::Float64
    q::Float64
end

"""
    ExponentialSelectionModel(η, q)

Selection coefficients are drawn from `Exponential(η)`.
"""
struct ExponentialSelectionModel <: SelectionModel
    η::Float64
    q::Float64
end

"""
    GammaSelectionModel(η, σ, q)

Selection coefficients are drawn from `Gamma(η²/σ², σ²/η)`.
"""
struct GammaSelectionModel <: SelectionModel
    η::Float64
    σ::Float64
    q::Float64
end

"""
    FixedSelectionModel(η, q)

All variants receive the same selection coefficient `η`.
"""
struct FixedSelectionModel <: SelectionModel
    η::Float64
    q::Float64
end

"""
    FreeFixedModel(η, q)

Identical to [`FixedSelectionModel`](@ref). Used in ABC fitting
contexts where `η` is a free parameter.
"""
struct FreeFixedModel <: SelectionModel
    η::Float64
    q::Float64
end

"""
    FixedSizeGrowthModel(selection)

Population size is held constant (Moran-like dynamics).
Drift includes competition term `s_vid[i] - s̄`.
"""
struct FixedSizeGrowthModel <: GrowthModel
    selection::SelectionModel
end

"""
    UnconstrainedGrowthModel(selection)

Population grows exponentially. No competition term in drift.
"""
struct UnconstrainedGrowthModel <: GrowthModel
    selection::SelectionModel
end

"""
    SimArgs

Mutable struct holding per-simulation parameter arrays.
Each field is a column vector with one entry per simulation.

Fields (indexed by simulation `sid`):
- `k_sid::Vector{Int64}` — number of variants
- `t₀_vid_Sid::Vector{Vector{Float64}}` — variant arrival times
- `_trackerID_Sid::Vector{Vector{Int}}` — indices of tracked variants
- `s_vid_Sid::Vector{Vector{Float64}}` — per-variant selection coefficients
- `init_vid_Sid::Vector{BitVector}` — whether each variant has been initiated
- `x₀_vid_Sid::Vector{Vector{Float64}}` — initial variant sizes
- `parentId_vid_Sid::Vector{Vector{Int64}}` — parent clone index per variant

Convert to a `DataFrame` via `DataFrame(simargs)`.
"""
mutable struct SimArgs
    k_sid::Vector{Int64}
    t₀_vid_Sid::Vector{Vector{Float64}}
    _trackerID_Sid::Vector{Vector{Int}}
    s_vid_Sid::Vector{Vector{Float64}}
    init_vid_Sid::Vector{BitVector}
    x₀_vid_Sid::Vector{Vector{Float64}}
    parentId_vid_Sid::Vector{Vector{Int64}}
end

function DataFrames.DataFrame(sa::SimArgs)
    return DataFrame(
        k=sa.k_sid, t₀_vid=sa.t₀_vid_Sid,
        _trackerID=sa._trackerID_Sid,
        s_vid=sa.s_vid_Sid, init_vid=sa.init_vid_Sid,
        x₀_vid=sa.x₀_vid_Sid, parentId_vid=sa.parentId_vid_Sid,
    )
end

"""
    complete(params::NamedTuple) -> NamedTuple

Fill in default values for optional simulation parameters.
User-provided keys override defaults.

Defaults set:
- `σ = 0.0`, `q = 0.0`
- `growthModel = "fixed size"`, `tMature = 0.0`
- `α = 1/τ` (inferred if not provided), `sMax = Inf`, `sortBias = 0.0`
"""
function complete(params::NamedTuple)
    α = haskey(params, :α) ? params.α : 1.0 / params.τ
    return merge(
        (;
            σ=0.0, q=0.0,
            growthModel="fixed size", tMature=0.0,
            α=α, sMax=Inf, sortBias=0.0,
        ),
        params,
    )
end

"""
    ALGS

Named tuple mapping solver name symbols to SDE solver constructors.
Available solvers: `:LambaEM`, `:SOSRI2`, `:SOSRA2`, `:SRA3`,
`:ImplicitEulerHeun`, `:SKenCarp`.
"""
const ALGS = (LambaEM=LambaEM, SOSRI2=SOSRI2, SOSRA2=SOSRA2, SRA3=SRA3, ImplicitEulerHeun=ImplicitEulerHeun, SKenCarp=SKenCarp)

include("simEvolver.jl")

export evolvePopSim, GrowthModel, SelectionModel, SimArgs,
       FixedSizeGrowthModel, UnconstrainedGrowthModel,
       GaussianSelectionModel, ExponentialSelectionModel, GammaSelectionModel,
       FixedSelectionModel, FreeFixedModel,
       prepareSims, complete, selectModel

end
