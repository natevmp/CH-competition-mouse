module CompetitiveSelection

using DataFrames
using Distances
using Statistics, StatsBase, LinearAlgebra
using DifferentialEquations
using StochasticDiffEq
using Interpolations
using Distributions
using Random
using JLD2

# ---- API -----
# simulation exports
export evolvePopSim
# analysis exports


abstract type GrowthModel end

abstract type SelectionModel end

struct GaussianSelectionModel <: SelectionModel
    s::Float64
    σ::Float64
    q::Float64
end

struct ExponentialSelectionModel <: SelectionModel
    s::Float64
    q::Float64
end

struct GammaSelectionModel <: SelectionModel
    s::Float64
    σ::Float64
    q::Float64
end

struct FixedSelectionModel <: SelectionModel
    s::Float64
    q::Float64
end

struct FreeFixedModel <: SelectionModel
    s::Float64
    q::Float64
end

struct FixedSizeGrowthModel <: GrowthModel 
    selection::SelectionModel
end

struct UnconstrainedGrowthModel <: GrowthModel
    selection::SelectionModel
end

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
    complete(params::NamedTuple)

Fill in default values for optional simulation parameters, so downstream code
can access every field.
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

include("simEvolver.jl")
include("scientist.jl")


end

