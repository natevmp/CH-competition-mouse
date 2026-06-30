module Scientist

using ..CompetitionSDE
using DifferentialEquations
using Distributions
using StatsBase
using LinearAlgebra
using DataFrames
using JLD2
using Random

function distbin(x_vid; bins::Int=25, xMin::Real=0, xMax::Union{Nothing,Real}=nothing, normalized::Bool=true)
    isnothing(xMax) && (xMax = maximum(x_vid)+0.01)
    _xEdges = range(xMin, xMax, bins+1)
    dx = _xEdges[2]-_xEdges[1]
    _x = _xEdges[1:end-1] .+ dx/2
    inRange_vid = (x -> x>=xMin && x<xMax).(x_vid)   #apply thresholds
    n_x = zeros(Float64, bins)
    for x in x_vid[inRange_vid]
        k = ((x-xMin)/dx |> floor |> Int) + 1
        n_x[k] += 1
    end
    normalized && (n_x .*= dx/sum(inRange_vid))
    return _x, n_x
end

"""
Get sizes distribution of variants at multiple time points determined by `bins`, from ensemble of simulations.
"""
function variantsSizeDistributionTimeBinned(solEns::EnsembleSolution; bins=25, xMin::Real=0.01, xMax::Union{Nothing,Real}=nothing, normalized::Bool=true)
    tLen = length(solEns[1])
    nSims = length(solEns)
    nV_t_f = Array{Float64,2}(undef, tLen, bins)
    _f = Vector{Float64}(undef, bins)
    for tInd in 1:tLen
        x_Sim_Vid = Array{Array{Float64}}(undef, nSims)
        for simId in 1:nSims
            x_vid = solEns[simId][tInd]
            x_Sim_Vid[simId] = x_vid
        end
        xGen_vid = vcat(x_Sim_Vid...)
        _f, nV_t_f[tInd,:] = distbin(xGen_vid; bins, xMin, xMax, normalized)

    end
    return _f, nV_t_f
end

"""
Get size distribution of variants at multiple time points determined by `bins`, from ensemble of simulations.
This method uses the StatsBase Histogram type, but is much slower than the main implementation.
"""
function variantsSizeDistributionAlt(solEns::EnsembleSolution; bins=25, xMin::Real=0.01, xMax::Real=1, normalized::Bool=true)
    normMode = normalized ? :pdf : :none
    tLen = length(solEns[1])
    nSims = length(solEns)
    nV_t_f = Array{Float64,2}(undef, tLen, bins)
    _edges = range(xMin, xMax, length=bins+1)
    _f = _edges[1:end-1] .+ (_edges[2]-_edges[1])/2
    for tInd in 1:tLen
        x_Sim_Vid = Array{Array{Float64}}(undef, nSims)
        for simId in 1:nSims
            x_vid = solEns[simId][tInd]
            x_Sim_Vid[simId] = x_vid
        end
        xGen_vid = vcat(x_Sim_Vid...)
        histData = LinearAlgebra.normalize(
            fit(Histogram, xGen_vid, _edges),
            mode=normMode
        )
        nV_t_f[tInd, :] = histData.weights
    end
    return _f, nV_t_f
end

function sampleIntFromFreq(vaf, nSamples)
    vaf<0 && (vaf=0)
    vaf>1 && (vaf=1)
    rand(Binomial(nSamples, vaf)) / nSamples
end

function variantsAboveThreshold(solEns::EnsembleSolution, f0::Real; nSamples::Union{Nothing,Real}=nothing)
    T = length(solEns[1])
    runs = length(solEns)
    nVarsAv_t = Vector{Float64}(undef, T)
    for tInd in eachindex(solEns[1])
        nVarsTot = 0
        for simId in eachindex(solEns)
            if !isnothing(nSamples)
                vafSampled_vid = (vaf->sampleIntFromFreq(vaf, nSamples)).(solEns[simId][tInd])
                detectableVariants_vid = vafSampled_vid .> f0
            else
                detectableVariants_vid = solEns[simId][tInd] .> f0
            end
            nVarsTot += sum(detectableVariants_vid)
        end
        nVarsAv_t[tInd] = nVarsTot / runs
    end
    return nVarsAv_t
end

function variantsAboveThreshold(solEns::EnsembleSolution, t::Real, f0::Real; nSamples::Union{Nothing,Real}=nothing)
    nVarsDetectAv = 0
    fTotDetectAv = 0
    for simId in eachindex(solEns)
        if !isnothing(nSamples)
            vafSampled_vid = (vaf->sampleIntFromFreq(vaf, nSamples)).(solEns[simId](t))
            detectableVariants_vid = vafSampled_vid .> f0
        else
            detectableVariants_vid = solEns[simId](t) .> f0
        end
        nVarsDetectAv += sum(detectableVariants_vid) / length(solEns)
        fTotDetectAv += sum(vafSampled_vid[detectableVariants_vid]) / length(solEns)
    end
    return nVarsDetectAv, fTotDetectAv
end

"""
    Sample variant frequencies of the simulation result at timepoints `_t` with sample sizes `S_t`. Discard variants under frequency `f0`. Return a list of the total variant fractions `fTotDetectAv_t` for each timepoint and the number of variants detected `nVarsDetectAv_t`.
    `sol`     = simulation result
    `S_t`     = sample size for each timepoint
    `_t`      = list of timepoints
    `f0`      = minimum sampled frequency to accept
    returns `nVarsDetectAv_t`, `fTotDetectAv_t`
"""
function variantsAboveThreshold(solEns::EnsembleSolution, _t::AbstractVector{<:Real}, f0::Real, S_t::AbstractVector{Int})
    if length(S_t) !== length(_t)
        println("error: `S_t` should be same length as `_t`")
    end
    nVarsDetectAv_t = Vector{Float64}(undef, length(_t))
    fTotDetectAv_t = Vector{Float64}(undef, length(_t))
    for (tid, t) in enumerate(_t)
        nVarsDetectAv_t[tid], fTotDetectAv_t[tid] = variantsAboveThreshold(solEns, t, f0; nSamples=S_t[tid])
    end
    return nVarsDetectAv_t, fTotDetectAv_t
end

function sampleSimTimepoint(sim::RODESolution, t::Real, nSamples::Int, f0::Float64)
    _cid = findall(sim(t).>=f0) |> shuffle!
    # make sure number of samples to draw does not exceed number of available variants
    length(_cid) > nSamples && (nSamples=length(_cid))
    vaf_cid = sim(t)[_cid[1:nSamples]]
    return vaf_cid
end

function drawVariantsRandomly(nVariants::Int, nVarsMax::Union{Nothing,Int})
    S = isnothing(nVarsMax) ? nVariants : minimum((nVarsMax, nVariants)) # if the number of variants is less than `nVarsMax`, use that
    vid_vidS = randperm(nVariants)[1:S]
    return vid_vidS
end

function chooseExclusiveRand(v_i, compare_j)
    vSelect = v_i[rand(1:end)]
    if !in(vSelect, compare_j)
        return vSelect
    else
        return chooseExclusiveRand(v_i, compare_j)
    end
end

function chooseExclusiveNext(i, v_i, compare_j)
    vSelect = v_i[i]
    if !in(vSelect, compare_j)
        return vSelect
    else
        if i+1>length(v_i) return vSelect end
        return chooseExclusiveNext(i+1, v_i, compare_j)
    end
end

function drawVariantsWithPriority(
        s_vid::AbstractVector{Float64},
        nVarsMax::Union{Nothing,Int}=nothing,
        sortBias::Real=1.,
    )
    vid_vFit = (1:length(s_vid))[sortperm(s_vid, rev=true)]
    S =
        if isnothing(nVarsMax)
            length(s_vid)
        else
            minimum((nVarsMax, length(s_vid))) # if the number of existing variants is less than `nVarsMax`, use that
        end
    vid_vidS = Vector{Int}(undef, S)
    for i in 1:S
        if rand()<=sortBias
            # select next fittest variant
            vid_vidS[i] = chooseExclusiveNext(i, vid_vFit, vid_vidS[1:i-1])
        else
            # select next random variant
            vid_vidS[i] = chooseExclusiveRand(vid_vFit, vid_vidS[1:i-1])
        end
    end

    return vid_vidS # variant ID's of sampled variants
end

"""
    Sample trajectories from a single simulated solution.
    Return _t, x_t_vidS
    t           = first timepoint to sample at
    tStep       = length of timesteps at which successive measurements are made
    nTimeSteps  = total number of measurements per trajectory
    freqCutoff  = ignore variants below this value
    s_vid       = if passed, variant innate fitness will be used to determine sampling priority
    sortBias    = determines the tendency of the sampler to bias towards the fittest mutants
"""
function sampleSimTrajectories(sol::RODESolution, t::Real;
    tStep=1.,
    nTimeSteps::Int=4,
    freqCutoff=0.01,
    nVarsMax::Union{Nothing,Int}=nothing,
    s_vid::Union{Vector{Float64},Nothing}=nothing, # if passed, variant innate fitness will be used to determine sampling priority
    sortBias::Real=0, # determines the tendency of the sampler to bias towards the fittest mutants,
    verbose::Bool=false
    )
    obs_vid = sol(t) .>= freqCutoff
    if verbose
        println("total variants: ", length(obs_vid))
        println("observable variants:", sum(obs_vid))
    end
    vidO_vidS =
        if isnothing(s_vid)
            drawVariantsRandomly(sum(obs_vid), nVarsMax)
        else
            verbose && println("drawing variants with priority bias=$(sortBias)...")
            drawVariantsWithPriority(s_vid[obs_vid], nVarsMax, sortBias)
        end
    _t = range(t; length=nTimeSteps, step=tStep)
    x_t_vidS = Array{Float64,2}(undef, nTimeSteps, length(vidO_vidS))
    for (i,tt) in enumerate(_t)
        for (j,x) in enumerate(sol(tt)[obs_vid][vidO_vidS])
            x_t_vidS[i,j] = x
        end
    end
    return _t, x_t_vidS
end

"""
    Sample trajectories from a single simulated solution, using x_t_vid as input. This is still kind of janky because times need to equal indices. Use with caution.
    Return _t, x_t_vidS
    t           = first timepoint to sample at. Must be int with index t+1
    tStep       = length of timesteps at which successive measurements are made. Must also be int.
    nTimeSteps  = total number of measurements per trajectory
    freqCutoff  = ignore variants below this value
    s_vid       = if passed, variant innate fitness will be used to determine sampling priority
    sortBias    = determines the tendency of the sampler to bias towards the fittest mutants
"""
function sampleSimTrajectories(x_t_vid,
    t::Int; # timepoint of first measurement (also index of x_t_vid)
    tStep::Int=1,
    nTimeSteps::Int=4,
    freqCutoff=0.01,
    nVarsMax::Union{Nothing,Int}=nothing,
    s_vid::Union{Vector{Float64},Nothing}=nothing, # if passed, variant innate fitness will be used to determine sampling priority
    sortBias::Real=0, # determines the tendency of the sampler to bias towards the fittest mutants. Only applies when s_vid is passed.
    verbose::Bool=false
    )
    obs_vid = x_t_vid[t+1,:] .>= freqCutoff
    if verbose
        println("total variants: ", length(obs_vid))
        println("observable variants:", sum(obs_vid))
    end
    vidO_vidS =
        if isnothing(s_vid)
            drawVariantsRandomly(sum(obs_vid), nVarsMax)
        else
            verbose && println("drawing variants with priority bias=$(sortBias)...")
            drawVariantsWithPriority(s_vid[obs_vid], nVarsMax, sortBias)
        end
    _t = range(t; length=nTimeSteps, step=tStep)
    x_t_vidS = Array{Float64,2}(undef, nTimeSteps, length(vidO_vidS))

    for (i,tt) in enumerate(_t)
        for (j,x) in enumerate(x_t_vid[tt,obs_vid][vidO_vidS])
            x_t_vidS[i,j] = x
        end
    end
    return _t, x_t_vidS
end

struct StepUniform <: Sampleable{Univariate,Continuous}
    edges::Vector{Float64}
    counts::Vector{Int64}
end

function Base.rand(rng::AbstractRNG, d::StepUniform)
    bin = StatsBase.sample(rng, range(1,length(d.counts)), Weights(d.counts))
    sample = d.edges[bin] + (d.edges[bin+1]-d.edges[bin])*rand()
    return sample
end

function averageTrackedVariant(solEns, simArgs)
    sizeTrackedvariant_sid_t = Array{Float64}(undef, length(solEns), length(solEns[1].t))
    for (sid, sol) in enumerate(solEns)
        trackID = simArgs[sid,:trackerID]
        for t in eachindex(sol.t)
            sizeTrackedvariant_sid_t[sid,t] = sol.u[t][trackID]
        end
    end
    return vec(mean(sizeTrackedvariant_sid_t, dims=1))
end

## ================================= ABC functions ===========================

"""
sample a variant allele frequency from a true frequency `x` with coverage `coverage`.
"""
function sampleFreqFromFreq(x, coverage)
    x<0 && (return 0.)
    x>1 && (return 1.)
    rand(Binomial(coverage, x))/coverage
end

"""
    function sampleSimTimepoints(sim::RODESolution, _t::AbstractVector, nSamples::Int, coverage::Int; 
        fMin::Float64=0.,
        fMax::Float64=1.,
        s_vid::Union{Nothing, AbstractVector{Float64}}=nothing,
        sortBias::Real=0.,
    )

Sample variant trajectories from a single simulation with binomial sampling on the VAFs at times `_t`. Returns the size trajectories of the observed variants `x_t_VidSampled`, with a maximal length of `nSamples`.
"""
function sampleSimTimepoints(sim::RODESolution, _t::AbstractVector, nSamples::Int, coverage::Int;
        fMin::Float64=0.,
        fMax::Float64=1.,
        s_vid::Union{Nothing, AbstractVector{Float64}}=nothing,
        sortBias::Real=0.,
        verbose::Bool=false,
    )
    # `x`∈[0,1]; not to be confused with `vaf`∈[0,0.5]
    x_vid_T = sim(_t)
    # select only observable variants
    obs_vid = fMin .<= x_vid_T[1] .<= fMax
    vidO_vidS =
        if !isnothing(s_vid)
            # we take 4*nSamples to allow for the possibility of variants being lost later during binomial sampling
            drawVariantsWithPriority(s_vid[obs_vid], 4*nSamples, sortBias)
        else
            drawVariantsRandomly(length(s_vid[obs_vid]), 4*nSamples)
        end
    j = 0   # index for accepted observed variants
    x_t_VidSampled = Vector{Vector{Float64}}(undef, nSamples)
    for vidO in vidO_vidS
        x_t = [x_vid[obs_vid][vidO] for x_vid in x_vid_T]
        xObs_t = (x -> sampleFreqFromFreq(x, coverage)).(x_t)
        #check whether variant has enough nonzero observations; if not skip to next variant (this is why _vidS has length 4*nSamples)
        sum(xObs_t.==0)>=length(_t)/2 && continue
        j+=1
        x_t_VidSampled[j] = xObs_t
        j==nSamples && break
    end
    if verbose #! debug
        println("sampleSimTimepoints::any(obs_vid): ", any(obs_vid))
        println("sampleSimTimepoints::vidO_vidS: ", vidO_vidS)
    end
    if j<nSamples
        return x_t_VidSampled[1:j]
    end
    return x_t_VidSampled
end

"""
    sampleSimTimepoint(sim::RODESolution, t::Real, nSamples::Int, coverage::Int; fMin::Float64=0., fMax::Float64=1.)

Samples variants from a single simulation with binomial sampling on the VAFs, at a single timepoint `t`.
"""
function sampleSimTimepoint(
        sim::RODESolution,
        vid_child_Vid::Vector{Vector{T}} where T<:Integer,
        t::Real,
        nSamples::Int,
        coverage::Int;
        fMin::Float64=0.,
        fMax::Float64=1.,
        s_vid::Union{Nothing, AbstractVector{Float64}}=nothing,
        sortBias::Real=0.,
    )
    obs_vid = fMin .<= sim(t) .<= fMax
    vidO_vidS =
        if isnothing(s_vid)
            drawVariantsRandomly(length(s_vid[obs_vid]), 4*nSamples)
        else
            drawVariantsWithPriority(s_vid[obs_vid], 4*nSamples, sortBias)
        end
    vafObs_j = Vector{Float64}(undef, nSamples)
    j = 0
    for vidO in vidO_vidS
        x = @view(sim(t)[obs_vid])[vidO]
        # add size of children (if they exist)
        for vid in vid_child_Vid[obs_vid][vidO]
            x += sim(t)[vid]
        end
        (x<fMin || x>fMax) && continue
        # allele frequency of a variant is f=x/2
        vafObs = rand(Binomial(coverage, x/2)) / coverage
        vafObs==0 && continue
        j+=1
        vafObs_j[j]=vafObs
        j==nSamples && break
    end
    if j<nSamples
        return vafObs_j[1:j]
    end
    return vafObs_j
end

function buildFamilyArray(parentVid_vid::Vector{T} where T<:Integer)
    childVid_child_Vid = [Int64[] for _ in eachindex(parentVid_vid)]
    for (vid, parId) in enumerate(parentVid_vid)
        if parId==0 continue end
        push!(childVid_child_Vid[parId], vid)
    end
    return childVid_child_Vid
end

"""
Compute the size distribution of `f_vid` as a density per bin.
"""
function sizeDistDens(x_vid; bins::Int=25, xMin=0., xMax=0.5)
    _xEdge = range(xMin, xMax, bins+1)
    dx = Float64(_xEdge.step)
    n_x = zeros(Float64, bins)
    for x in x_vid
        k = ((x-xMin)/dx |> floor |> Int) + 1
        if k > bins
            n_x[end] += 1.
            continue
        end
        n_x[k] += 1.
    end
    nDens_x = (n -> n*dx/sum(n_x)).(n_x)
    return _xEdge, nDens_x
end

function sizeDistSims(
    tMeasure,
    solEns,
    parentVid_vid_Sid::Vector{Vector{T}} where T<:Integer,
    s_vid_Sid::Vector{Vector{Float64}},
    ctrlParams;
    sortBias=0.,
    )
    # -------- sample sims --------
    vaf_vidS = Vector{Float64}(undef, ctrlParams[:simRuns]*ctrlParams[:nSamples])
    vCur = 1
    vNext = 1
    for (sid,sol) in enumerate(solEns)
        childVid_child_Vid = buildFamilyArray(parentVid_vid_Sid[sid])
        f_vidSimSample = sampleSimTimepoint(
                sol,
                childVid_child_Vid,
                tMeasure,
                ctrlParams[:nSamples],
                ctrlParams[:coverage];
                fMin = 2*ctrlParams[:fMin],
                fMax = 2*ctrlParams[:fMax],
                s_vid = s_vid_Sid[sid],
                sortBias = sortBias,
            )
        vNext += length(f_vidSimSample)
        vaf_vidS[vCur:vNext-1] .= f_vidSimSample
        vCur = vNext
    end
    # -------- construct size distribution --------
    _f, n_f = sizeDistDens(@view vaf_vidS[1:vCur-1]; bins=ctrlParams[:fBins], xMin=ctrlParams[:fMin], xMax=ctrlParams[:fMax])
    if ctrlParams[:cumulativeDist]
        # construct cumulative size distribution
        n_f .= hcat([sum(@view n_f[1:i]) for i in eachindex(n_f)])
    else
        n_f .= hcat(n_f)
    end
    return n_f
end



"""
    Sample variant frequencies of the simulation result at timepoints `_t` with sample sizes `S_t`. Discard variants under frequency `f0`. Return a list `s_vid_T` with for each timepoint a list of the sampled variant frequencies above f0.
    `sol`     = simulation result
    `S_t`     = sample size for each timepoint
    `_t`      = list of timepoints
    `f0`      = minimum sampled frequency to accept
    returns `s_vid_T`
"""
function samplePatientSim(sol, S_t::Vector{Int}, _t; f0=0)
    nVars = length(sol[1])
    s_vid_T = Vector{Vector{Float64}}(undef, length(_t))
    for (tInd,t) in enumerate(_t)
        s_vid_T[tInd] = Float64[]
        for vid in 1:nVars
            # draw S times with success prob fVid to get frequency of vid sVid in sample
            binomDist = Binomial(S_t[tInd], sol(t)[vid])
            sVid = rand(binomDist) / S_t[tInd]
            if sVid>f0 push!(s_vid_T[tInd], sVid) end
        end
    end
    return s_vid_T
end

function samplePatientSim(sol, S::Int, _t; f0=0)
    S_t = fill(S, length(_t))
    samplePatientSim(sol, S_t, _t; f0)
end

function runModelSimFixedFitness(paramsABC, ctrlParams)
    # fix `s` based on value of τ
    ctrlParams[:params] = merge(ctrlParams[:params], (s=ctrlParams[:sFixed] * paramsABC[:τ],))
    runModelSim(paramsABC, ctrlParams)
end

## ----------------------------------------
#region - ABC modelfitting


"""
    runModelSim(paramsABC)

Perform a run of the model simulations to obtain a single particle with parameter set `paramsABC`.
"""
function runModelSim(paramsABC, ctrlParams; debug::Bool=false)
    modelParams = ctrlParams[:params]
    for (pName, pVal) in pairs(paramsABC)
        pName ∈ ctrlParams[:fixPar] && continue
        modelParams = merge(modelParams, (; pName => pVal))
    end

    if ctrlParams[:normalizedFitnessDist]
        modelParams = merge(modelParams, (s=modelParams.s * modelParams.τ, σ=modelParams.σ * modelParams.τ))
    end

    # run model sims
    selection = GammaSelectionModel(modelParams.s, modelParams.σ, modelParams.q)
    growthModel = FixedSizeGrowthModel(selection)
    solEns, simArgs = evolvePopSim(modelParams, growthModel; runs=ctrlParams[:simRuns], noDiffusion=false)
    tMeasure = (ctrlParams[:tBounds][1]+ctrlParams[:tBounds][2])/2
    parentVid_vid_Sid = simArgs[!, :parentId_vid]
    ctrlParams[:params] = modelParams

    # ============ size distribution ============
    sortBias =
        if haskey(ctrlParams, :sortBias)
            ctrlParams[:sortBias]
        else
            modelParams.sortBias
        end
    nVars_f = sizeDistSims(
        tMeasure,
        solEns,
        parentVid_vid_Sid,
        simArgs[!,:s_vid],
        ctrlParams;
        sortBias=sortBias,
    )

    #! ---- debug ----
    if debug
        jldsave("debug.jld2"; tMeasure, solEns, simArgs, ctrlParams, sortBias)
    end

    return nVars_f, nVarsDetectAv_t, fTotDetectAv_t
end

#endregion

export distbin, variantsSizeDistributionTimeBinned, variantsSizeDistributionAlt,
       sampleIntFromFreq, variantsAboveThreshold,
       sampleSimTimepoint, drawVariantsRandomly, drawVariantsWithPriority,
       sampleSimTrajectories, averageTrackedVariant,
       sampleFreqFromFreq, sampleSimTimepoints,
       buildFamilyArray, sizeDistDens, sizeDistSims,
       samplePatientSim, runModelSimFixedFitness, runModelSim

end
