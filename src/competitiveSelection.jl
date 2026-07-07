module CompetitiveSelection

include("CompetitionSDE/competitionSDE.jl")
include("OldScientist/oldScientist.jl")

using .CompetitionSDE
using .OldScientist

export evolvePopSim

end
