module CompetitiveSelection

include("CompetitionSDE/competitionSDE.jl")
include("Scientist/scientist.jl")

using .CompetitionSDE
using .Scientist

export evolvePopSim

end
