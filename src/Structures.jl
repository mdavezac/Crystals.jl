module Structures
export AbstractCrystal, Crystal, is_fractional, volume
using Unitful: Quantity, Dimensions, Units, unit, ustrip

using DataFrames: DataFrame, nrow, NA, index, ncol, deleterows!
using Crystals: Log
import Base
import Unitful
import DataFrames
using Unitful: dimension, unit

""" Top type node for Crystals """
abstract AbstractCrystal

typealias RowIndices{T <: Integer} Union{T, AbstractVector{T}, Range{T}, Colon}

""" Crystal that lies in a space with physical units """
type Crystal{T <: Number, D, U, P <: Number} <: AbstractCrystal
    """ Periodicity of the crystal structure """
    cell::Matrix{Quantity{T, D, U}}
    """ Atomic positions """
    positions::Matrix{P}
    """ Atomic properties """
    properties::DataFrame
end

function Crystal{C, D, U, P <: Number}(cell::Matrix{Quantity{C, D, U}},
                                       positions::Matrix{P},
                                       args...; kwargs...)
    size(cell, 1) == size(cell, 2) || Log.error("Cell matrix is not square")

    if size(positions, 1)  != size(cell, 1)
        Log.error("Positions and cell have different sizes")
    end

    properties = DataFrame(args...; kwargs...)

    if nrow(properties) != 0 && nrow(properties) != size(positions, 2)
        Log.error("atomic properties and positions have different lengths")
    end

    if P <: Quantity
        const Q = promote_type(Quantity{C, D, U}, P)
        p = convert(Matrix{Q}, positions)
        Crystal{Q.parameters..., Q}(convert(Matrix{Q}, cell),
                                    convert(Matrix{Q}, positions),
                                    properties)
    else
        const T = promote_type(C, P)
        const Q = Quantity{T, D, U}
        Crystal{T, D, U, T}(convert(Matrix{Q}, cell),
                            convert(Matrix{T}, positions),
                            properties)
    end
end

function Crystal{C, D, U}(cell::Matrix{Quantity{C, D, U}}; kwargs...)
    const dpositions = filter(x -> x[1] ∈ (:position, :positions), kwargs)
    const tpositions = filter(x -> x[1] ∈ (:tposition, :tpositions), kwargs)
    if length(dpositions) ≠ 0 && length(tpositions) ≠ 0
        const positions = hcat((x[2] for x in dpositions)...,
                               transpose(hcat([x[2] for x in tpositions]...)))
    elseif length(dpositions) ≠ 0
        const positions = hcat([x[2] for x in dpositions]...)
    elseif length(tpositions) ≠ 0
        const positions = transpose(hcat([x[2] for x in tpositions]...))
    else
        const positions = reshape(Matrix{Quantity{C, D, U}}(),
                                  (size(cell, 1), 0))
    end

    filter!(x -> x[1] ∉ (:position, :positions, :tposition, :tpositions), kwargs)
    Crystal(cell, positions; kwargs...)
end

""" Underlying physical units of the crystal """
Unitful.unit(crystal::Crystal) = typeof(crystal).parameters[3]()
Unitful.dimension(crystal::Crystal) = typeof(crystal).parameters[2]()
Base.length(crystal::Crystal) = size(crystal.positions, 2)
is_fractional(crystal::Crystal) = is_fractional(typeof(crystal))
is_fractional{T <: Crystal}(::Type{T}) = !(T.parameters[end] <: Quantity)
fractional_trait(crystal::Crystal) = fractional_trait(typeof(crystal))
function fractional_trait{T <: Crystal}(::Type{T})
    is_fractional(T) ? Val{:fractional}(): Val{:real}()
end

volume(cell::Matrix) = abs(det(cell))
function volume{T, D, U}(cell::Matrix{Quantity{T, D, U}})
    abs(det(ustrip(cell))) * unit(eltype(cell))^3
end
volume(crystal::Crystal) = volume(crystal.cell)


"""
position_for_crystal(crystal::Crystal, position::Array)

Converts position to same type as in the input crystal. If real positions are
given for a crystal with fractional positions, then the positions are converted.
And vice-versa.
"""
function position_for_crystal(crystal::Crystal, position::Array)
    position_for_crystal(fractional_trait(crystal), crystal.cell, position)
end
function position_for_crystal{T <: Quantity}(::Val{:fractional},
                                             cell::Matrix,
                                             position::Array{T})
    inv(cell) * position
end
position_for_crystal(::Val{:fractional}, cell::Matrix, position::Array) = position
position_for_crystal(::Val{:real}, cell::Matrix, position::Array) = cell * position
position_for_crystal{T <: Quantity}(::Val{:real}, cell::Matrix, pos::Array{T}) = pos

@inline _is_not_position(s::Symbol) = :position ≠ s
@inline _is_not_position(s::AbstractVector{Symbol}) = :position ∉ s

function Base.show(io::IO, crystal::Crystal)
    println(io, "cell(", unit(crystal), "):")
    for i in 1:size(crystal.cell, 1)
        print(io, "  ")
        join(io, ustrip(crystal.cell[i, :]), ' ')
        println(io)
    end
    with_pos = copy(crystal.properties)
    const name = is_fractional(crystal) ? :fractional : :position
    const positions = [tuple(ustrip(crystal.positions[:, i])...) for i in 1:length(crystal)]
    show(io, hcat(DataFrame(Any[positions], [name]), crystal.properties),
         false, :Atom, false)
end

function Base.push!(crystal::Crystal, position::Vector; kwargs...)
    if length(position) ≠ size(crystal.cell, 1)
        Log.error("Dimensionality of input position and crystal do not match")
    end
    crystal.positions = hcat(crystal.positions,
                             position_for_crystal(crystal, position))
    row = Any[NA for u in 1:length(crystal.properties)]

    missing = Tuple{Symbol, Any}[]
    const colnames = names(crystal.properties)
    for (name, value) in kwargs
        if name ∉ colnames
            push!(missing, (name, value))
        else
            row[index(crystal.properties)[name]] = value
        end
    end
    if size(crystal.properties, 2) == 0 && length(missing) > 0
        const name, value = pop!(missing)
        crystal.properties[name] = [value]
    else
        push!(crystal.properties, row)
    end
    for (name, value) in missing
        # given twice, makes no sense
        @assert name ∉ names(crystal.properties)
        crystal.properties[name] = fill(value, size(crystal.properties, 1))
        crystal.properties[1:(size(crystal.properties, 1) - 1), name] = NA
    end
end

function Base.getindex(crystal::Crystal, index::RowIndices)
    typeof(crystal)(crystal.cell, hcat(crystal.positions[:, index]),
                    crystal.properties[index, :])
end

function Base.getindex(crystal::Crystal, symbol::Symbol)
    if symbol == :position
        return crystal.positions
    end
    crystal.properties[symbol]
end

Base.getindex(crystal::Crystal, ::Colon) = deepcopy(crystal)

function Base.getindex(crystal::Crystal, symbols::AbstractVector{Symbol})
    if :position ∉ symbols
        return crystal.properties[symbols]
    end
    args = filter(x -> x ≠ :position, symbols)
    typeof(crystal)(crystal.cell, crystal.positions, crystal.properties[args])
end

function Base.getindex(crystal::Crystal, index::RowIndices, symbol::Symbol)
    if symbol == :position
        return crystal.positions[:, index]
    end
    crystal.properties[index, symbol]
end

function Base.getindex(crystal::Crystal, ::Colon, symbol::Symbol)
    if symbol == :position
        return crystal.positions
    end
    crystal.properties[:, symbol]
end

Base.getindex(crystal::Crystal, ::Colon, ::Colon) = copy(crystal)
Base.getindex(crystal::Crystal, row::Any, col::Colon) = Base.getindex(crystal, row)

function Base.getindex(crystal::Crystal, index::RowIndices, symbols::AbstractVector{Symbol})
    if :position ∈ symbols
        args = filter(x -> x ≠ :position, symbols)
        return typeof(crystal)(crystal.cell,
                               hcat(crystal.positions[:, index]),
                               crystal.properties[index, args])
    end
    crystal.properties[index, symbols]
end

function Base.getindex(crystal::Crystal, symbol::Symbol, pos::Any)
    @assert symbol == :position
    crystal.positions[pos, :]
end

function Base.getindex(crystal::Crystal, atoms::Any, symbol::Symbol, pos::Any)
    @assert symbol == :position
    crystal.positions[pos, atoms]
end

Base.names(crystal::Crystal) = push!(names(crystal.properties), :position)
Base.size(crystal::Crystal) = (r = size(crystal.properties); (r[1], r[2] + 1))
Base.size(crystal::Crystal, i::Integer) = size(crystal)[i]
Base.ndims(crystal::Crystal) = 2
DataFrames.nrow(crystal::Crystal) = size(crystal.positions, 2)
DataFrames.ncol(crystal::Crystal) = ncol(crystal.properties) + 1
Base.endof(crystal::Crystal) = nrow(crystal)

function Base.copy(crystal::Crystal)
    typeof(crystal)(copy(crystal.cell), copy(crystal.positions), copy(crystal.properties))
end

function Base.deepcopy(crystal::Crystal)
    typeof(crystal)(deepcopy(crystal.cell),
                    deepcopy(crystal.positions),
                    deepcopy(crystal.properties))
end

function Base.setindex!(crystal::Crystal, v::DataFrame, cols::AbstractVector{Symbol})
    if :position ∈ names(v) || :position ∈ cols
        Log.error("Positions cannot be set from a dataframe")
    end
    setindex!(crystal.properties, v, cols)
end

function Base.setindex!(crystal::Crystal,
                        v::DataFrame,
                        rows::Any,
                        cols::AbstractVector{Symbol})
    if :position ∈ names(v) || :position ∈ cols
        Log.error("Positions cannot be set from a dataframe")
    end
    setindex!(crystal.properties, v, rows, cols)
end

function Base.setindex!(crystal::Crystal, v::Crystal, cols::AbstractVector{Symbol})
    if :position ∈ cols
        # check lattices are equivalent
        @assert isinteger(inv(crystal.cell) * v.cell)
        @assert volume(crystal) ≈ volume(v)
        real_position = position_for_crystal(Val{:real}(), v.cell, v.positions)
        crystal.positions = position_for_crystal(crystal, real_position)
        # forward to other properties
        setindex!(crystal.properties, v.properties, filter(x -> x ≠ :position, cols))
    else
        setindex!(crystal.properties, v.properties, cols)
    end
end

function Base.setindex!(crystal::Crystal, v::Any, row::Any, col::Symbol)
    if col == :position
        setindex!(crystal.positions, position_for_crystal(crystal, v), :, row)
    else
        setindex!(crystal.properties, v, row, col)
    end
end

function Base.setindex!(crystal::Crystal, v::Any, col::Symbol)
    if col == :position
        if size(v) ≠ size(crystal.positions)
            Log.error("Input has incorrect size")
        end
        crystal.positions = position_for_crystal(crystal, v)
    else
        setindex!(crystal.properties, v, col)
    end
end

function Base.setindex!(crystal::Crystal, v::Crystal,
                        row::Any, cols::AbstractVector{Symbol})
    if :position ∈ cols
        # check lattices are equivalent
        @assert isinteger(inv(crystal.cell) * v.cell)
        @assert volume(crystal) ≈ volume(v)
        real_position = position_for_crystal(Val{:real}(), v.cell, v.positions)
        crystal.positions[:, row] = position_for_crystal(crystal, real_position)
        # forward to other properties
        setindex!(crystal.properties, v.properties, row, filter(x -> x ≠ :position, cols))
    else
        setindex!(crystal.properties, v.properties, row, cols)
    end
end

function Base.setindex!(crystal::Crystal, v::Any, col::Symbol, xyz::RowIndices)
    if col != :position
        Log.error("crystal[:positions, integer] is only available for positions")
    end

    crystal.positions[xyz, :] = v
end

function Base.setindex!(crystal::Crystal, v::Any, rows::RowIndices,
                        col::Symbol, xyz::RowIndices)
    if col != :position
        Log.error("crystal[:positions, integer] is only available for positions")
    end
    crystal.positions[xyz, rows] = v
end

Base.setindex!(crys::Crystal, v::Any, ::Colon, col::Symbol) = Base.setindex!(crys, v, col)

function Base.delete!(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}})
    if !_is_not_position(col)
        Log.error("Cannot delete position column from a structure")
    end
    delete!(crystal.properties, col)
    crystal
end

Base.delete!(crystal::Crystal, ::Colon) = (empty!(crystal.properties); crystal)

Base.setindex!(crystal::Crystal, ::Void, col::Any) = delete!(crystal, col)

function Base.empty!(crystal::Crystal)
    empty!(crystal.properties)
    crystal.positions = crystal.positions[:, 1:0]
end

function DataFrames.deleterows!(crystal::Crystal, row::Integer)
    deleterows!(crystal.properties, row)
    rows = filter(i -> i ≠ row, 1:nrow(crystal))
    crystal.positions = crystal.positions[:, rows]
end

function DataFrames.deleterows!(crystal::Crystal, rows::RowIndices)
    deleterows!(crystal.properties, rows)
    prows = filter(i -> i ∉ rows, 1:nrow(crystal))
    crystal.positions = crystal.positions[:, prows]
end

# # Special deletion assignment
# Base.setindex!(crystal::Crystal, x::Void, col::Int) = delete!(crystal, col)
# Base.empty!(crystal::Crystal) = (empty!(crystal.atoms); crystal)
# Base.insert!(crystal::Crystal, col::Int, item::AbstractVector, name::Symbol) =
#     (insert!(crystal.atoms, col, item, name); crystal)
#
# Base.merge!(crystal::Crystal, others::DataFrame...) =
#     (merge!(crystal.atoms, others...); crystal)

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
# DataFrames.eachrow(crystal::Crystal) = eachrow(crystal.atoms)
# DataFrames.eachcol(crystal::Crystal) = eachcol(crystal.atoms)
#

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
