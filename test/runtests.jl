module CrystalTest
using Crystals
using FactCheck: @fact, facts, context, exitstatus, roughly,
                 exactly, @fact_throws, greater_than, not
using DataFrames: nrow, ncol, NA, DataArray, DataFrame, deleterows!, NA

contains(x) = y -> x ∈ y
all_integers(x::Array, ε::AbstractFloat=1e-8) = all(abs(x - round(Integer, x)) .< ε)
all_integers(ε::AbstractFloat=1e-8) = y -> all_integers(y, ε)

include("Crystal.jl")
include("SpaceGroup.jl")
include("Gruber.jl")
include("SNF.jl")
include("utilities.jl")

exitstatus()
end
