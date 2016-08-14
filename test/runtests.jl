module CrystalTest
using Crystals
using FactCheck: @fact, facts, context, exitstatus, roughly,
                 exactly, @fact_throws, greater_than
using DataFrames: nrow, ncol, NA, DataArray, DataFrame

contains(x) = y -> x ∈ y

# include("Crystal.jl")
# include("SpaceGroup.jl")
# include("Gruber.jl")
include("SNF.jl")

exitstatus()
end
