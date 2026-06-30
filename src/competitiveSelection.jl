module CompetitiveSelection

using DataFrames
using Distances
using Statistics, StatsBase, LinearAlgebra
using DifferentialEquations
using StochasticDiffEq
using Interpolations
using Distributions
using Random
using Parameters
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

include("simEvolver.jl")
include("scientist.jl")


end

