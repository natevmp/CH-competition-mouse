
"""
    prepareSims(params, selectionModel, _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}, runs) where {U,S<:Real}

Create the `DataFrame` `simArgs`. Each row contains the arguments pertaining to a single simulation. The columns denote the arguments:

`k`: total number of variants arising \\
`t₀_vid`: `vector` of clone arrival times ordered by clone id (coincindes with ascending arrival time order). \\
`_trackerID`: `vector` with elements pertaining to inserted tracker variants. Each element is a `tuple` `(t0, s)` containg that variant's arrival time `t0` and innate fitness `s`. \\
`init_vid`: `vector` of bools denoting whether a clone has arrived (relevant for internal simulation). \\
`parentId_vid`: `vector` denoting the parent of each clone, ordered by clone id. \\
`s_vid`: `vector` with elements the innate fitness of each clone, ordered by clone id. \\
"""
function prepareSims(params, selectionModel, _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}, runs::Int; growthPhase::Bool=false) where {U,S<:Real}
    (; T, μ, N) = params
    tMature = params.tMature

    k_sid = Vector{Int64}(undef, runs)
    t₀_vid_Sid = Vector{Vector{Float64}}(undef, runs)
    for sid in eachindex(k_sid)
        kGrowthPhase, t₀_vidPB = 
            if growthPhase
                growthPhaseArrivals((N=N, μ=μ, tMature=tMature))
            else
                (0, Float64[])
        end
        kConstPhase = rand( Poisson(μ*(T-tMature)) )
        t₀_vidConstPhase = rand(Uniform(tMature,T), kConstPhase) |> sort!
        k_sid[sid] = kGrowthPhase + kConstPhase
        t₀_vid_Sid[sid] = vcat(t₀_vidPB, t₀_vidConstPhase)
    end
    _trackerID_Sid = [Vector{Int}(undef, length(_trackerVariant)) for _ in eachindex(k_sid)]
    for sid in eachindex(k_sid)
        for ts in _trackerVariant
            k_sid[sid] += 1
            push!(t₀_vid_Sid[sid], ts[1])
        end
        sort!(t₀_vid_Sid[sid])
        for (i,ts) in enumerate(_trackerVariant)
            _trackerID_Sid[sid][i] = findfirst(t₀_vid_Sid[sid] .== ts[1])
        end
    end
    s_vid_Sid = (k->selectionParamModel(selectionModel, k)).(k_sid)
    if isfinite(params.sMax)
       for s_vid in s_vid_Sid
           maxS_vid = findall(s_vid.>params.sMax)
           s_vid[maxS_vid] .= params.sMax
       end
    end
    init_vid_Sid = falses.(k_sid)
    x₀_vid_Sid = zeros.(k_sid)
    parentId_vid_Sid = zeros.(Int64, k_sid)

    simArgs = SimArgs(k_sid, t₀_vid_Sid, _trackerID_Sid, s_vid_Sid, init_vid_Sid, x₀_vid_Sid, parentId_vid_Sid)

    if length(_trackerVariant)<1 return simArgs end
    if eltype(_trackerVariant)==Tuple{Real} return simArgs end
    # place tracked variants
    for sid in eachindex(simArgs.s_vid_Sid)
        simArgs.s_vid_Sid[sid][simArgs._trackerID_Sid[sid]] .= [ts[2] for ts in _trackerVariant]
    end
    return simArgs
end

function drawInhomogeneousPoisson(λ::Function, λHom::Real, tWindow::Tuple{<:Real,<:Real})
    t0 = tWindow[1]
    T = tWindow[end]
    mHom = rand(Poisson(λHom))
    t_vidHom = rand(Uniform(t0,T), mHom) # generate locations in homogeneous process
    t_vid = Float64[]   # locations in inhomogeneous process
    for t in t_vidHom
        if rand() > λ(t)/λHom
            continue # skip to thin point
        end
        push!(t_vid, t) # add point to inhomogeneous process
    end
    sort!(t_vid)
    return t_vid
end

function growthPhaseArrivals(params)
    (; μ, N, tMature) = params
    TPreBirth = 9/12
    nPop(t) = exp(log(N)/(TPreBirth+tMature) * t)
    λ(t) = μ/N * nPop(t)
    λMax = μ*(TPreBirth+tMature)
    tPos_vid = drawInhomogeneousPoisson(λ, λMax, (0,TPreBirth+tMature))
    t_vid = @. 0 - (TPreBirth - tPos_vid)
    return (k=length(tPos_vid), t_vid=t_vid)
end

function selectionParamModel(selectionModel::FixedSelectionModel, k::Int)
    fill(selectionModel.s, k)
end

function selectionParamModel(selectionModel::FreeFixedModel, k::Int)
    fill(selectionModel.s, k)
end

function selectionParamModel(selectionModel::GaussianSelectionModel, k::Int)
    rand(Normal(selectionModel.s, selectionModel.σ), k)
end

function selectionParamModel(selectionModel::ExponentialSelectionModel, k::Int)
    rand(Exponential(selectionModel.s), k)
end

function selectionParamModel(selectionModel::GammaSelectionModel, k::Int)
    rand(Gamma(selectionModel.s^2/selectionModel.σ^2 , selectionModel.σ^2/selectionModel.s), k)
end

function fitnessDoubleHit(sChild, sParent)
    sParent + sChild + sChild*sParent
end

"""
    Sets the parent clone of a new variant, and optionally performs a double hit if `q>0`.
"""
function initiateVariant!(childId::Integer, init_vid, s_vid, x_vid, parentId_vid::Vector{Int64}, q; sMax::Real=Inf)
    pRand = rand()  # random variable ∈(0,1) to select the parent
    xTot = 0.
    for (parentId,x) in enumerate(@view x_vid[init_vid])
        xTot += x
        if xTot<pRand continue end # continue if current ID is not parent
        parentId_vid[childId] = parentId
        if (q==0 || rand()>q) break end # break if double hit fails
        s_vid[childId] = min(fitnessDoubleHit(s_vid[childId], s_vid[parentId]), sMax)
        break
    end
    init_vid[childId] = true
end

function drift(model::FixedSizeGrowthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t; sMax::Real=Inf)
    sX = 0   # non-allocating sum
    for (vid, x) in enumerate(x_vid)
        t < t0_vid[vid] && break
        sX += x*s_vid[vid]
    end
    # @views sX = sum(x_vid .* s_vid) # allocating sum
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break # skip variant if not yet initiated
        !init_vid[i] && initiateVariant!(i, init_vid, s_vid, x_vid, parentId_vid, model.selection.q; sMax=sMax) # initiate variant if new
        if x<0
            x_vid[i] = 0
            dx_vid[i] = 0
            continue
        end
        if x>1
            x_vid[i] = 1
            dx_vid[i] = 0
            continue
        end
        dx_vid[i] = α*x*(s_vid[i] - sX) # evolve variant
    end
end

function diffusion(model::FixedSizeGrowthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t)
    sX = 0  # non-allocating sum
    for (vid, x) in enumerate(x_vid)
        t < t0_vid[vid] && break
        sX += x*s_vid[vid]
    end
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        dx_vid[i] = (x>0 && x<1) ? √( α*(2 + s_vid[i] + sX)*x/N ) : 0
    end
end

function drift(model::UnconstrainedGrowthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t; sMax::Real=Inf)
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break # skip variant if not yet initiated
        if !init_vid[i] # initiate variant if new
            initiateVariant!(i, init_vid, s_vid, x_vid, parentId_vid, model.selection.q; sMax=sMax)
        end
        if x<0
            x_vid[i] = 0
            dx_vid[i] = 0
            continue
        end
        dx_vid[i] = α*x*s_vid[i] # evolve variant
    end
end

function diffusion(model::UnconstrainedGrowthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s_vid, N, parentId_vid), t)
    for (i, x) in enumerate(x_vid)
        t < t0_vid[i] && break
        dx_vid[i] = (x>0) ? √( α*x*(2 + s_vid[i]) ) : 0
    end
end

function newCloneSize(model::UnconstrainedGrowthModel, N)
    1.
end

function newCloneSize(model::FixedSizeGrowthModel, N)
    1/N
end

function evolveGrowthPhase!(
    simArgs::SimArgs,
    params,
    selectionModel::SelectionModel;
    runs::Int,
    _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}=Vector{Float64}[],
    noDiffusion::Bool=false,
    algorithm::Symbol=:LambaEM,
    ) where {U,S<:Real}

    (; N, s, μ, α, tMature, sMax) = params
    T = 9/12 # full time is from conception until birth

    growthModel = UnconstrainedGrowthModel(selectionModel)
    f!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t) = 
        drift(growthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t; sMax)
    g!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t) = begin
        noDiffusion ? 0 : diffusion(growthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t)
    end

    # get allocation-free minimum
    t0Min = -T + 0.1
    for t_vid in simArgs.t₀_vid_Sid
        if length(t_vid)==0 continue end # ingnore empty arrays
        if t_vid[1]<t0Min
            t0Min = t_vid[1]
        end
    end

    condt0(u,t,integrator) = t==t0Min
    function addStops!(integrator)
        for t0 in integrator.p[1]
            if t0<tMature
                add_tstop!(integrator, t0)
            end
        end
    end
    callBackAddStops = DiscreteCallback(condt0, addStops!; save_positions=(false,false))
    function affect!(integrator)
        i = findfirst(integrator.p[1].==integrator.t)
        integrator.u[i] = newCloneSize(growthModel, integrator.p[5]) # set size of new variant
    end
    condition(u,t,integrator) = t ∈ integrator.p[1]
    callBackAddClone = DiscreteCallback(condition, affect!; save_positions=(false,false))
    callbacks = CallbackSet(callBackAddStops, callBackAddClone)

    prob = SDEProblem(f!, g!, simArgs.x₀_vid_Sid[1], (-T, tMature), (simArgs.t₀_vid_Sid[1], simArgs.init_vid_Sid[1], α, simArgs.s_vid_Sid[1], N, simArgs.parentId_vid_Sid[1]))
    solver = ALGS[algorithm]
    function probFunc(prob, ctx)
        i = ctx.sim_id
        remake(prob, u0=simArgs.x₀_vid_Sid[i], p=(simArgs.t₀_vid_Sid[i], simArgs.init_vid_Sid[i], α, simArgs.s_vid_Sid[i], N, simArgs.parentId_vid_Sid[i]))
    end
    ensembleProb = EnsembleProblem(prob, prob_func=probFunc)
    solEns = solve(ensembleProb, solver(), EnsembleThreads(); callback=callbacks, tstops=[t0Min,], saveat=[tMature,], trajectories=runs)
    x0_vid_Sid = Vector{Vector{Float64}}(undef, runs)
    for (sid, sol) in enumerate(solEns.u)
        x0_vid_Sid[sid] = sol.u[end] ./ (sum(sol.u[end]) + N)
    end
    simArgs.x₀_vid_Sid = x0_vid_Sid
    return solEns
end

function selectModel(sType::String, s, σ, q, growthModel::String)
    selection =
        if sType=="fixed"
            FixedSelectionModel(s, q)
        elseif sType=="exponential"
            ExponentialSelectionModel(s, q)
        elseif sType=="gaussian"
            GaussianSelectionModel(s, σ, q)
        elseif sType=="gamma"
            GammaSelectionModel(s, σ, q)
        elseif sType=="free"
            FreeFixedModel(s, q)
        else
            error("Error: selection model undefined.")
        end
    growth =
        if growthModel=="fixed size"
            FixedSizeGrowthModel(selection)
        elseif growthModel=="unconstrained"
            UnconstrainedGrowthModel(selection)
        else
            error("Error: growth model undefined.")
        end
    return selection, growth
end

"""
    evolvePopSim()

`_trackerVariant` is a vector containing clones to be tracked. The elements are either the birth times, or a tuple containing both birth time and fitness of the form `(t0, s)`.
"""
function evolvePopSim(
        params;
        runs::Int=1,
        _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}=Vector{Float64}[],
        noDiffusion::Bool=false,
        algorithm::Symbol=:LambaEM,
        growthPhase::Bool=false,
    ) where {U,S<:Real}
    params = complete(params)
    (; sType, s, σ, q, growthModel) = params
    selectionModel, growthModel = selectModel(sType, s, σ, q, growthModel)
    evolvePopSim(params, growthModel; runs, _trackerVariant, noDiffusion, algorithm, growthPhase)
end

function evolvePopSim(
    params,
    growthModel::GrowthModel; 
    runs::Int=1,
    _trackerVariant::Union{Vector{U},Vector{Tuple{U,S}}}=Vector{Float64}[],
    noDiffusion::Bool=false,
    algorithm::Symbol=:LambaEM,
    growthPhase::Bool=false,
    ) where {U,S<:Real}
    
    (; N, s, T, μ, α, tMature) = params

    f!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t) = 
        drift(growthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t; sMax=params.sMax)
    g!(dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t) = begin
        noDiffusion ? 0 : diffusion(growthModel, dx_vid, x_vid, (t0_vid, init_vid, α, s, N, parentId_vid), t)
    end

    simArgs = prepareSims(params, growthModel.selection, _trackerVariant, runs; growthPhase)

    if growthPhase
        evolveGrowthPhase!(simArgs, params, growthModel.selection; runs, _trackerVariant, noDiffusion, algorithm)
    end

    # Get first post-growth variant arrival time from all sims
    # (allocation-free method)
    t0Min = tMature + 1.
    for t_vid in simArgs.t₀_vid_Sid
        if length(t_vid)==0 continue end # ingnore empty arrays (no variants occur)
        ind = findfirst(t_vid.>tMature) # ignore developmental arrival times
        if t_vid[ind]<t0Min
            t0Min = t_vid[ind]
        end
    end
    condt0(u,t,integrator) = t==t0Min
    function addStops!(integrator)
        for t0 in integrator.p[1]
            if t0<tMature continue end # skip prenatal stops
            add_tstop!(integrator, t0)
        end
    end
    callBackAddStops = DiscreteCallback(condt0, addStops!; save_positions=(false,false))
    function affect!(integrator)
        i = findfirst(integrator.p[1].==integrator.t)
        integrator.u[i] = newCloneSize(growthModel, integrator.p[5]) # set size of new variant
    end
    condition(u,t,integrator) = t ∈ integrator.p[1]
    callBackAddClone = DiscreteCallback(condition, affect!; save_positions=(false,false))
    callbacks = CallbackSet(callBackAddStops, callBackAddClone)

    prob = SDEProblem(
        f!,
        g!,
        simArgs.x₀_vid_Sid[1],
        (tMature, T),
        (simArgs.t₀_vid_Sid[1], simArgs.init_vid_Sid[1], α, simArgs.s_vid_Sid[1], N, simArgs.parentId_vid_Sid[1])
    )
    solver = ALGS[algorithm]
    function probFunc(prob, ctx)
        i = ctx.sim_id
        remake(prob, u0=simArgs.x₀_vid_Sid[i], p=(simArgs.t₀_vid_Sid[i], simArgs.init_vid_Sid[i], α, simArgs.s_vid_Sid[i], N, simArgs.parentId_vid_Sid[i]))
    end
    ensembleProb = EnsembleProblem(prob, prob_func=probFunc)
    _t = range(Int(ceil(tMature)), T)
    solEns = solve(ensembleProb, solver(), EnsembleThreads(); callback=callbacks, tstops=[t0Min,], saveat=_t, trajectories=runs)
    return solEns, DataFrame(simArgs)
end
