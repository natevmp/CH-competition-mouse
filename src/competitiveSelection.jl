module CompetitiveSelection

include("CompetitionSDE/competitionSDE.jl")
include("Scientist/scientist.jl")

using .CompetitionSDE
using .Scientist

export evolvePopSim, GrowthModel, SelectionModel, SimArgs, FixedSizeGrowthModel,
       UnconstrainedGrowthModel, FixedSelectionModel, ExponentialSelectionModel,
       GaussianSelectionModel, GammaSelectionModel, FreeFixedModel,
       prepareSims, complete,
       distbin, variantsSizeDistributionTimeBinned, variantsSizeDistributionAlt,
       sampleIntFromFreq, variantsAboveThreshold,
       sampleSimTimepoint, drawVariantsRandomly, drawVariantsWithPriority,
       sampleSimTrajectories, averageTrackedVariant,
       sampleFreqFromFreq, sampleSimTimepoints,
       buildFamilyArray, sizeDistDens, sizeDistSims,
       samplePatientSim, runModelSimFixedFitness, runModelSim

end
