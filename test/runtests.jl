import Pkg; Pkg.activate(".")

includet("../src/competitiveSelection.jl")
using .CompetitiveSelection
using DataFrames
import DifferentialEquations as DE
println("=== smoke test: evolvePopSim ===")

params = (;
    sType="gamma",
    η=0.1,
    σ=0.02,
    q=0.5,
    N=1000,
    T=100,
    μ=2,
    τ=1.,
    growthModel="fixed size",
)

solEns, simArgs = evolvePopSim(params; runs=10, algorithm=:LambaEM);

println("nSims: ", length(solEns.u))
println("simArgs rows: ", nrow(simArgs))
println("simArgs cols: ", names(simArgs))

@assert length(solEns.u) == 10 "expected 10 simulations"
@assert nrow(simArgs) == 10 "expected 10 sim arg rows"

println("PASSED: evolvePopSim returns correct types and shapes")

println("\n=== smoke test: evolvePopSim with growthPhase ===")

solEns2, simArgs2 = evolvePopSim(params; runs=5, growthPhase=true, algorithm=:LambaEM)

println("nSims: ", length(solEns2.u))
println("simArgs rows: ", nrow(simArgs2))

@assert length(solEns2.u) == 5
@assert nrow(simArgs2) == 5

println("PASSED: evolvePopSim with growthPhase")

println("\n=== smoke test: evolvePopSim with _trackerVariant (1-element tuples) ===")

params2 = merge(params, (T=1,))  # shorter T for tracked variant tests
_tv1 = [(0.0,), (0.5,)]
solEns3, simArgs3 = evolvePopSim(params2; runs=3, _trackerVariant=_tv1, tSaves=10)
@assert length(solEns3.u) == 3
@assert length(simArgs3._trackerID[1]) == 2 "expected 2 tracker IDs, got $(simArgs3._trackerID[1])"
@assert all(1 .<= simArgs3._trackerID[1] .<= simArgs3.k[1]) "tracker IDs out of range"
println("PASSED: 1-element tuples — no error, $(nrow(simArgs3)) sims, IDs=$(simArgs3._trackerID[1])")

println("\n=== smoke test: evolvePopSim with _trackerVariant (2-element tuples, distinct times) ===")

_tv2 = [(0.001, 0.3), (0.5, 0.15)]
solEns4, simArgs4 = evolvePopSim(params2; runs=3, _trackerVariant=_tv2, tSaves=10)
@assert length(solEns4.u) == 3
@assert length(simArgs4._trackerID[1]) == 2
s_tracked = simArgs4.s_vid[1][simArgs4._trackerID[1]]
@assert sort(s_tracked) == [0.15, 0.3] "tracked s_vid should match assigned values: got $s_tracked"
u_mat = hcat(solEns4.u[1].u...)  # [nVars × nSaves]
@assert all(i -> any(u_mat[i,:] .> 0), simArgs4._trackerID[1]) "each tracked variant should be >0 at some save point"
println("PASSED: 2-element tuples — s_vid matches, variants grow")

println("\n=== smoke test: evolvePopSim with _trackerVariant (5 variants, same time) ===")

_tv3 = fill((0.001, 0.4), 5)
solEns5, simArgs5 = evolvePopSim(params2; runs=2, _trackerVariant=_tv3, tSaves=5)
@assert length(solEns5.u) == 2
@assert length(simArgs5._trackerID[1]) == 5
s_tracked5 = simArgs5.s_vid[1][simArgs5._trackerID[1]]
@assert all(==(0.4), s_tracked5) "all tracked s_vid should be 0.4, got $s_tracked5"
u_mat5 = hcat(solEns5.u[1].u...)
@assert all(i -> any(u_mat5[i,:] .> 0), simArgs5._trackerID[1]) "each tracked variant should be >0 at some save point"
println("PASSED: 5 variants at same time — all initialized correctly")

println("\n=== _trackerID verification: distinct times ===")

_tv6 = [(0.001, 0.1), (0.5, 0.2), (0.9, 0.3)]
solEns6, simArgs6 = evolvePopSim(params2; runs=1, _trackerVariant=_tv6, tSaves=5)
ids = simArgs6._trackerID[1]
@assert length(unique(ids)) == 3 "trackerIDs must be unique for distinct times: $ids"
s_tracked6 = sort(simArgs6.s_vid[1][ids])
@assert s_tracked6 == [0.1, 0.2, 0.3] "s_vid mapping wrong: $s_tracked6"
println("PASSED: distinct times → unique IDs, correct s_vid")

println("\n=== _trackerID verification: same time ===")

_tv7 = fill((0.001, 0.4), 4)
solEns7, simArgs7 = evolvePopSim(params2; runs=1, _trackerVariant=_tv7, tSaves=5)
ids7 = simArgs7._trackerID[1]
@assert length(unique(ids7)) == 4 "trackerIDs must all be unique at same time: $ids7"
@assert all(==(0.4), simArgs7.s_vid[1][ids7]) "s_vid mapping wrong"
println("PASSED: same time → all unique IDs, correct s_vid")

println("\n=== _trackerID verification: mixed times (some same, some distinct) ===")

_tv8 = [(0.001, 0.1), (0.001, 0.2), (0.5, 0.3)]
solEns8, simArgs8 = evolvePopSim(params2; runs=1, _trackerVariant=_tv8, tSaves=5)
ids8 = simArgs8._trackerID[1]
@assert length(unique(ids8)) == 3 "trackerIDs must be unique for mixed times: $ids8"
s_sorted8 = sort(simArgs8.s_vid[1][ids8])
@assert s_sorted8 == [0.1, 0.2, 0.3] "s_vid mapping wrong: $s_sorted8"
println("PASSED: mixed times → all unique IDs, correct s_vid")

println("\nAll tests passed.")
