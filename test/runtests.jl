module CrystalTest
using Crystals
using FactCheck: @fact, facts, context, exitstatus, roughly,
                 exactly, @fact_throws, greater_than, not
using DataFrames: nrow, ncol, NA, DataArray, DataFrame, deleterows!, NA,
                  DataFrameRow, eachrow, @data

contains(x) = y -> x ∈ y
all_integers(x::Array, ε::AbstractFloat=1e-8) = all(abs(x - round(Integer, x)) .< ε)
all_integers(ε::AbstractFloat=1e-8) = y -> all_integers(y, ε)
is_subtype(x::Type) = y -> (y <: x)

# include("Positions.jl")
# include("Crystal.jl")
include("SpaceGroup.jl")
# include("SNF.jl")
# include("utilities.jl")
# include("Gruber.jl")

exitstatus()
end
