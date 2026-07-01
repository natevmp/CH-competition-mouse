import Pkg; Pkg.activate(".")

includet("../src/competitiveSelection.jl")
using .CompetitiveSelection
using DataFrames
import DifferentialEquations as DE
println("=== smoke test: evolvePopSim ===")

params = (;
    sType="gamma",
    s=0.1,
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

println("\nAll tests passed.")
