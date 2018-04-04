module CrystalTest
using Crystals
using DataFrames: nrow, DataFrame, ncol, deleterows!, missing, ismissing, eachrow
using Base.Test
using Unitful
using Documenter

contains(x) = y -> x ∈ y
all_integers(x::Array, ε::AbstractFloat=1e-8) = all(abs(x - round(Integer, x)) .< ε)
all_integers(ε::AbstractFloat=1e-8) = y -> all_integers(y, ε)
is_subtype(x::Type) = y -> (y <: x)

@testset "Crystal" begin
    include("Crystal.jl")
end
@testset "Atoms and iteration" begin
    include("CrystalAtoms.jl")
end
@testset "SNF" begin
    include("SNF.jl")
end
@testset "Utilities" begin
    include("utilities.jl")
end
@testset "SpaceGroup" begin
    include("SpaceGroup.jl")
end
@testset "Gruber" begin
    include("Gruber.jl")
end

@testset "Docstests" begin
    mktempdir() do path
        makedocs(modules = [Crystals], clean = true,
             build = joinpath(path, "build"),
             pages = Crystals._doc_pages(),
             root=joinpath(dirname(@__DIR__), "docs"))
    end
end

end
