module Structure
export AbstractCrystal, Crystal, RealCrystal, FractionalCrystal, is_fractional
using Unitful: Quantity, Dimensions, Units, unit, ustrip

using DataFrames: DataFrame, nrow
using Crystals: Log
import Base
import Unitful

"""
Base type for all crystals

Crystal should all be structures with a cell, positions and atomic properties.
The exact implementation thereof is left to the derived classes.
"""
abstract Crystal

"""
Crystal that lies in a space with physical units
"""
abstract UnitfulCrystal{T <: Number, D <: Dimensions, U <: Units} <: Crystal

""" Crystal with positions in real coordinates """
type RealCrystal{T, D, U} <: UnitfulCrystal{T, D, U}
    """ Periodicity of the crystal structure """
    cell::Matrix{Quantity{T, D, U}}
    """ Atomic positions """
    positions::Matrix{Quantity{T, D, U}}
    """ Atomic properties """
    properties::DataFrame

    function RealCrystal{T, D, U}(cell::Matrix{Quantity{T, D, U}},
                                  positions::Matrix{Quantity{T, D, U}},
                                  args...; kwargs...)
        size(cell, 1) == size(cell, 2) ||
            Log.error("Cell matrix is not square")

        size(positions, 1)  == size(cell, 1) ||
            Log.error("Positions and cell have different sizes")

        properties = DataFrame(args...; kwargs...)

        nrow(properties) == 0 || nrow(properties) == size(positions, 2) ||
            Log.error("atomic properties and positions have different lengths")

        new(cell, positions, properties)
    end
end

""" Crystal with positions in fractional coordinates """
type FractionalCrystal{T, D, U} <: UnitfulCrystal{T, D, U}
    """ Periodicity of the crystal structure """
    cell::Matrix{Quantity{T, D, U}}
    """ Atomic positions """
    positions::Matrix{T}
    """ Atomic properties """
    properties::DataFrame

    function FractionalCrystal{T, D, U}(cell::Matrix{Quantity{T, D, U}},
                                        positions::Matrix{T},
                                        args...; kwargs...)
        size(cell, 1) == size(cell, 2) ||
            Log.error("Cell matrix is not square")

        size(positions, 1)  == size(cell, 1) ||
            Log.error("Positions and cell have different sizes")

        properties = DataFrame(args...; kwargs...)

        nrow(properties) == 0 || nrow(properties) == size(positions, 2) ||
            Log.error("atomic properties and positions have different lengths")

        new(cell, positions, properties)
    end
end

function Crystal{T, D, U, Tp, Up}(cell::Matrix{Quantity{T, D, U}},
                                  positions::Matrix{Quantity{Tp, D, Up}},
                                  args...; kwargs...)
    const TT = promote_type(T, Tp)
    RealCrystal{TT, D, U}(convert(Matrix{Quantity{TT, D, U}}, cell),
                          convert(Matrix{Quantity{TT, D, U}}, positions),
                          args...; kwargs...)
end
Crystal(cell::Matrix, positions::Vector, args...; kwargs...) =
   Crystal(cell, reshape(positions, (length(positions), 1)), args...; kwargs...)

function Crystal{T, D, U}(cell::Matrix{Quantity{T, D, U}},
                          positions::Matrix{Quantity{T, D, U}},
                          args...; kwargs...)
    RealCrystal{T, D, U}(cell, positions, args...; kwargs...)
end

function Crystal{T, D, U, Tp <: Number}(cell::Matrix{Quantity{T, D, U}},
                                        positions::Matrix{Tp},
                                        args...; kwargs...)
    const TT = promote_type(T, Tp)
    FractionalCrystal{TT, D, U}(convert(Matrix{Quantity{TT, D, U}}, cell),
                                convert(Matrix{TT}, positions),
                                args...; kwargs...)
end

function Crystal{T, D, U}(cell::Matrix{Quantity{T, D, U}};
                          tpositions::Matrix=Matrix{Quantity{T, D, U}}(),
                          kwargs...)
    if length(tpositions) ≠ 0
        return Crystal(cell, transpose(tpositions); kwargs...)
    end
    Crystal(cell, zeros(Quantity{T, D, U}, (size(cell, 1), 0)); kwargs...)
end

""" Underlying physical units of the crystal """
Unitful.unit(crystal::UnitfulCrystal) = typeof(crystal).parameters[3]()
Unitful.dimension(crystal::UnitfulCrystal) = typeof(crystal).parameters[2]()
Base.length(crystal::Crystal) = size(crystal.positions, 2)
is_fractional(::FractionalCrystal) = true
is_fractional(::RealCrystal) = false

"""
     position_for_crystal(crystal::Crystal, position::Array)

Converts position to same type as in the input crystal. If real positions are
given for a crystal with fractional positions, then the positions are converted.
And vice-versa.
"""
position_for_crystal{T, D, U}(crystal::RealCrystal,
                              position::Array{Quantity{T, D, U}}) =
    convert(Array{eltype(crystal.positions)}, position)
position_for_crystal{T <: Number}(crystal::RealCrystal, position::Array{T}) =
    convert(Array{eltype(crystal.positions)}, crystal.cell * position)
position_for_crystal{T, D, U}(crystal::FractionalCrystal,
                              position::Array{Quantity{T, D, U}}) =
    convert(Array{eltype(crystal.positions)}, inv(crystal.cell) * position)
position_for_crystal{T <: Number}(crystal::Crystal, position::Array{T}) =
    convert(Array{eltype(crystal.positions)}, position)

function Base.show(io::IO, crystal::Crystal)
    println(io, "cell(", unit(crystal), "):")
    for i in 1:size(crystal.cell, 1)
        print(io, "  ")
        join(io, ustrip(crystal.cell[i, :]), ' ')
        println(io)
    end
    with_pos = copy(crystal.properties)
    const name = is_fractional(crystal) ? :fractional : :position
    const positions = [
        tuple(ustrip(crystal.positions[:, i])...) for i in 1:length(crystal)]
    show(io, hcat(DataFrame(Any[positions], [name]), crystal.properties),
         false, :Atom, false)
end

function Base.push!(crystal::Crystal, position::Vector; kwargs...)
    length(position) ≠ size(crystal.cell, 1) &&
        Log.error("Dimensionality of input position and crystal do not match")
    row = Any[NA for u in 1:length(crystal.atoms)]
    row[index(crystal.atoms)[:position]] =
        convert(eltype(crystal.atoms[:position]), position)

    missing = Tuple{Symbol, Any}[]
    const colnames = names(crystal)
    for (name, value) in kwargs
        if name ∉ colnames
            push!(missing, (name, value))
        else
            row[index(crystal.atoms)[name]] = value
        end
    end
    push!(crystal, row)
    for (name, value) in missing
        # given twice, makes no sense
        @assert name ∉ names(crystal)
        crystal[name] = fill(value, size(crystal, 1))
        crystal[1:(size(crystal, 1) - 1), name] = NA
    end
end
#
# function Base.push!(crystal::Crystal, row::DataFrameRow;
#                     no_new_properties::Bool=false)
#     push!(crystal, [NA for u in 1:ncol(crystal)])
#     for (key, value) in row[intersect(names(row), names(crystal))]
#         crystal[end, key] = value
#     end
#     no_new_properties && return
#     for (key, value) in row[setdiff(names(row), names(crystal))]
#         crystal[:, key] = value
#         crystal[1:end - 1, key] = NA
#     end
# end
#
# function Base.push!{T <: Number}(
#             crystal::Crystal, position::Vector{T}; kwargs...)
#     pos = convert(Position{T}, position)
#     Base.push!(crystal, pos; kwargs...)
# end

#
# # Forwards all indexing to the DataFrame
# """ Equivalent to `getindex(crystal.atoms, ...)` """
# Base.getindex(crystal::Crystal, col::Any) = crystal.atoms[col]
# Base.getindex(crystal::Crystal, col::Any, row::Any) = crystal.atoms[col, row]
# Base.names(crystal::Crystal) = names(crystal.atoms)
# """ Equivalent to `endof(crystal.atoms)` """
# Base.endof(crystal::Crystal) = endof(crystal.atoms)
# """ Equivalent to `size(crystal.atoms)` """
# Base.size(crystal::Crystal) = size(crystal.atoms)
# Base.size(crystal::Crystal, i) = size(crystal.atoms, i)
# """ Equivalent to `ndims(crystal.atoms)` """
# Base.ndims(crystal::Crystal) = ndims(crystal.atoms)
#
# # Single column
# Base.setindex!(crystal::Crystal, v::Crystal, col::Any) =
#     setindex!(crystal.atoms, v.atoms, col)
# Base.setindex!(crystal::Crystal, v::Matrix, col::Any) =
#     setindex!(crystal.atoms, convert(PositionDataArray, v), col)
# Base.setindex!(crystal::Crystal, v::Any, col::Any) =
#     setindex!(crystal.atoms, v, col)
# Base.setindex!(crystal::Crystal, v::Crystal, row::Any, col::Any) =
#     setindex!(crystal.atoms, v.atoms, row, col)
# Base.setindex!(crystal::Crystal, v::Matrix, row::Any, col::Any) =
#     setindex!(crystal.atoms, convert(PositionDataArray, v), row, col)
# Base.setindex!(crystal::Crystal, v::Any, row::Any, col::Any) =
#     setindex!(crystal.atoms, v, row, col)
#
# # Special deletion assignment
# Base.setindex!(crystal::Crystal, x::Void, col::Int) = delete!(crystal, col)
# Base.empty!(crystal::Crystal) = (empty!(crystal.atoms); crystal)
# Base.insert!(crystal::Crystal, col::Int, item::AbstractVector, name::Symbol) =
#     (insert!(crystal.atoms, col, item, name); crystal)
#
# Base.merge!(crystal::Crystal, others::DataFrame...) =
#     (merge!(crystal.atoms, others...); crystal)
# Base.copy!(crystal::Crystal) =
#     Crystal{eltype(crystal.cell)}(copy(crystal.cell), copy(crystal.atoms))
# Base.deepcopy(crystal::Crystal) =
#     Crystal{eltype(crystal.cell)}(
#         deepcopy(crystal.cell), deepcopy(crystal.atoms))
#
# Base.delete!(crystal::Crystal, cols::Any) =
#     (delete!(crystal.atoms, cols); crystal)
# DataFrames.deleterows!(crystal::Crystal, cols::Any) =
#     (deleterows!(crystal.atoms, cols); crystal)
# DataFrames.hcat!(crystal::Crystal, x...) = (hcat!(crystal.atoms, x...); crystal)
# Base.hcat(crystal::Crystal, x...) = hcat!(copy(crystal), x...)
# Base.vcat(crystal::Crystal) = crystal
# """
# `vcat(crys::Crystal, dfs::AbstractDataFrame...)`
#
# Concatenates atoms into crystal
# """
# function Base.vcat(crys::Crystal, dfs::AbstractDataFrame...)
#     result = Crystal(crys.cell)
#     result.atoms = vcat(crys.atoms, dfs...)
#     result
# end
# DataFrames.nullable!(crystal::Crystal, x...) =
#     (nullable!(crystal.atoms, x...); crystal)
# DataFrames.pool!(crystal::Crystal, x::Any) = pool!(crystal.atoms, x::Any)
# Base.append!(crystal::Crystal, atoms::DataFrame) =
#     (append!(crystal.atoms, atoms); crystal)
# Base.push!(crystal::Crystal, x::Any) = push!(crystal.atoms, x)
#
#
# function Base.show(io::IO, crystal::Crystal, args...; kwargs...)
#     println(io, typeof(crystal))
#     println(io, "cell(", unit(crystal.cell[1]), "): ",
#         ustrip(convert(Matrix{typeof(crystal.cell[1])}, crystal.cell)))
#     println(io, crystal.atoms)
# end
#
#
# function Base.showcompact(io::IO, pos::Position)
#     result = string(pos)
#     print(io, result[findfirst(result, '('):end])
# end
# function Base.showcompact{T, D, U}(io::IO, pos::Position{Quantity{T, D, U}})
#     result = string(ustrip(pos))
#     print(io, result[findfirst(result, '('):end], "(", unit(pos[1]), ")")
# end
#
# function DataFrames.ourshowcompact(io::IO, pos::Position)
#     result = string(pos)
#     print(io, result[findfirst(result, '(') + 1:end - 1])
# end
# function DataFrames.ourshowcompact{T, D, U}(
#         io::IO, pos::Position{Quantity{T, D, U}})
#     result = string(ustrip(pos))
#     print(io, result[findfirst(result, '(') + 1:end - 1], " ", unit(pos[1]))
# end
#
# DataFrames.nrow(crystal::Crystal) = nrow(crystal.atoms)
# DataFrames.ncol(crystal::Crystal) = ncol(crystal.atoms)
# DataFrames.eachrow(crystal::Crystal) = eachrow(crystal.atoms)
# DataFrames.eachcol(crystal::Crystal) = eachcol(crystal.atoms)
#
# volume(cell::Matrix) = abs(det(cell))
# volume(crystal::Crystal) = volume(crystal.cell)
#
# """
#     round!(crystal::Crystal, args...)
#
# Rounds the cell and position of a crystal. See `round` for possible parameters.
# """
# function round!(crystal::Crystal, args...)
#     crystal.cell = round(crystal.cell, args...)
#     crystal[:position] = round(convert(Array, crystal[:position]), args...)
#     crystal
# end
#
# """
#     round(crystal::Crystal, args...)
#
# Rounds the cell and position of a crystal. See `round` for possible parameters.
# """
# Base.round(crystal::Crystal, args...) = round!(deepcopy(crystal), args...)

end
