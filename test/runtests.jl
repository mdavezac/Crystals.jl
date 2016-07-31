module CrystalTest
using Crystals
using FactCheck: @fact, facts, context, exitstatus, roughly, exactly, @fact_throws
using DataFrames: nrow, ncol, NA, DataArray, DataFrame

contains(x) = y -> x ∈ y

include("Crystal.jl")
include("SpaceGroup.jl")

exitstatus()
end
