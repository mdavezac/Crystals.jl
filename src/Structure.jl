module Structure
export Crystal, volume

using DataFrames: AbstractDataFrame, DataFrame, NA, index, nrow, ncol, hcat!,
        nullable!, pool!, eachrow, eachcol, deleterows!, DataFrameRow
using Crystals.Positions: Position, PositionDataArray,
        MIN_POSITION_SIZE, MAX_POSITION_SIZE
using Crystals: Log
import DataFrames
import Base

""" Base type for all crystals """
abstract AbstractCrystal

"""
Crystal structure

Holds all data necessary to define a crystal structure.
The cell and scale are held as fields.
The atomic properties, including the positions, are held in a DataFrame.
"""
type Crystal{T} <: AbstractCrystal
    """ Scale/unit of the positions and cell """
    scale::T
    """ Periodicity of the crystal structure """
    cell::Matrix{T}
    """ Atoms and atomic properties """
    atoms::DataFrame

    function Crystal(cell, scale, atoms::DataFrame)
        size(cell, 1) == size(cell, 2) ||
            Log.error("Cell matrix is not square")
        MIN_POSITION_SIZE ≤ size(cell, 1) ≤ MAX_POSITION_SIZE ||
            Log.error("Incorrect column vector size")
        :position ∈ names(atoms) || nrow(atoms) == 0 ||
            Log.error("Input Dataframe has atoms without positions")
        if nrow(atoms) == 0
            atoms[:position] = Vector{typeof(Position(cell[:, 1]))}()
        else
            eltype(atoms[:position]) <: Position ||
                Log.error("Positions do not have acceptable type")
            size(cell, 1) == length(eltype(atoms[:position])) ||
                Log.error("Cell and position dimensionality do not match")
        end
        new(scale, cell, atoms)
    end
end

function Crystal(T::Type, cell::Matrix, scale=1::Real; kwargs...)
    nkwargs = Tuple{Symbol, Any}[]
    for i in 1:length(kwargs)
        if kwargs[i][1] ∉ (:position, :tposition)
            push!(nkwargs, kwargs[i])
            continue
        end
        if kwargs[i][1] == :position
            position = convert(PositionDataArray{T}, kwargs[i][2])
        else
            position = convert(PositionDataArray{T}, transpose(kwargs[i][2]))
        end
        size(cell, 1) == length(eltype(position)) ||
            Log.error("Dimensionality of cell and positions do not match")
        push!(nkwargs, (:position, position))
    end

    atoms = DataFrame(; nkwargs...)
    Crystal{T}(convert(Matrix{T}, cell), convert(T, scale), atoms)
end

Crystal(cell::Matrix, scale=1::Real; kwargs...) =
    Crystal(eltype(cell), cell, scale; kwargs...)


function Crystal(cell::Matrix, columns::Vector{Any}, names::Vector{Symbol},
                 scale=1::Real)
    for (i, (value, name)) in enumerate(zip(columns[:], names[:]))
        name == :position || continue

        columns[i] = convert(PositionDataArray, value)
        size(cell, 1) == length(eltype(columns[i])) ||
            Log.error("Dimensionality of cell and positions do not match")
    end

    atoms = DataFrame(columns, names)
    Crystal{eltype(cell)}(cell, scale, atoms)
end

# Forwards all indexing to the DataFrame
""" Equivalent to `getindex(crystal.atoms, ...)` """
Base.getindex(crystal::Crystal, col::Any) = crystal.atoms[col]
Base.getindex(crystal::Crystal, col::Any, row::Any) = crystal.atoms[col, row]
Base.names(crystal::Crystal) = names(crystal.atoms)
""" Equivalent to `endof(crystal.atoms)` """
Base.endof(crystal::Crystal) = endof(crystal.atoms)
""" Equivalent to `size(crystal.atoms)` """
Base.size(crystal::Crystal) = size(crystal.atoms)
Base.size(crystal::Crystal, i) = size(crystal.atoms, i)
""" Equivalent to `ndims(crystal.atoms)` """
Base.ndims(crystal::Crystal) = ndims(crystal.atoms)

# Single column
Base.setindex!(crystal::Crystal, v::Crystal, col::Any) =
    setindex!(crystal.atoms, v.atoms, col)
Base.setindex!(crystal::Crystal, v::Matrix, col::Any) =
    setindex!(crystal.atoms, convert(PositionDataArray, v), col)
Base.setindex!(crystal::Crystal, v::Any, col::Any) =
    setindex!(crystal.atoms, v, col)
Base.setindex!(crystal::Crystal, v::Crystal, row::Any, col::Any) =
    setindex!(crystal.atoms, v.atoms, row, col)
Base.setindex!(crystal::Crystal, v::Matrix, row::Any, col::Any) =
    setindex!(crystal.atoms, convert(PositionDataArray, v), row, col)
Base.setindex!(crystal::Crystal, v::Any, row::Any, col::Any) =
    setindex!(crystal.atoms, v, row, col)

# Special deletion assignment
Base.setindex!(crystal::Crystal, x::Void, col::Int) = delete!(crystal, col)
Base.empty!(crystal::Crystal) = (empty!(crystal.atoms); crystal)
Base.insert!(crystal::Crystal, col::Int, item::AbstractVector, name::Symbol) =
    (insert!(crystal.atoms, col, item, name); crystal)

Base.merge!(crystal::Crystal, others::DataFrame...) =
    (merge!(crystal.atoms, others...); crystal)
Base.copy!(crystal::Crystal) =
    Crystal{eltype(crystal.cell)}(
        copy(crystal.cell), crystal.scale, copy(crystal.atoms))
Base.deepcopy(crystal::Crystal) =
    Crystal{eltype(crystal.cell)}(
        deepcopy(crystal.cell), crystal.scale, deepcopy(crystal.atoms))

Base.delete!(crystal::Crystal, cols::Any) =
    (delete!(crystal.atoms, cols); crystal)
DataFrames.deleterows!(crystal::Crystal, cols::Any) =
    (deleterows!(crystal.atoms, cols); crystal)
DataFrames.hcat!(crystal::Crystal, x...) = (hcat!(crystal.atoms, x...); crystal)
Base.hcat(crystal::Crystal, x...) = hcat!(copy(crystal), x...)
Base.vcat(crystal::Crystal) = crystal
"""
`vcat(crys::Crystal, dfs::AbstractDataFrame...)`

Concatenates atoms into crystal
"""
function Base.vcat(crys::Crystal, dfs::AbstractDataFrame...)
    result = Crystal(crys.cell, crys.scale)
    result.atoms = vcat(crys.atoms, dfs...)
    result
end
DataFrames.nullable!(crystal::Crystal, x...) =
    (nullable!(crystal.atoms, x...); crystal)
DataFrames.pool!(crystal::Crystal, x::Any) = pool!(crystal.atoms, x::Any)
Base.append!(crystal::Crystal, atoms::DataFrame) =
    (append!(crystal.atoms, atoms); crystal)
Base.push!(crystal::Crystal, x::Any) = push!(crystal.atoms, x)

function Base.push!(crystal::Crystal, position::Position; kwargs...)
    length(position) ≠ size(crystal.cell, 1) &&
        Log.error("Dimensionality of input position and crystal do not match")
    row = Any[NA for u in 1:length(crystal.atoms)]
    row[index(crystal.atoms)[:position]] = position
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

function Base.push!(crystal::Crystal, row::DataFrameRow;
                    no_new_properties::Bool=false)
    push!(crystal, [NA for u in 1:ncol(crystal)])
    for (key, value) in row[intersect(names(row), names(crystal))]
        crystal[end, key] = value
    end
    no_new_properties && return
    for (key, value) in row[setdiff(names(row), names(crystal))]
        crystal[:, key] = value
        crystal[1:end - 1, key] = NA
    end
end

function Base.push!{T <: Real}(crystal::Crystal, position::Vector{T}; kwargs...)
    pos = convert(Position{T}, position)
    Base.push!(crystal, pos; kwargs...)
end

function Base.show(io::IO, crystal::Crystal, args...; kwargs...)
    println(io, typeof(crystal))
    println(io, "cell: ", crystal.cell)
    println(io, "scale: ", crystal.scale)
    println(io, crystal.atoms)
end


function Base.showcompact(io::IO, pos::Position)
    result = string(pos)
    print(io, result[findfirst(result, '('):end])
end

function DataFrames.ourshowcompact(io::IO, pos::Position)
    result = string(pos)
    print(io, result[findfirst(result, '(') + 1:end - 1])
end

DataFrames.nrow(crystal::Crystal) = nrow(crystal.atoms)
DataFrames.ncol(crystal::Crystal) = ncol(crystal.atoms)
DataFrames.eachrow(crystal::Crystal) = eachrow(crystal.atoms)
DataFrames.eachcol(crystal::Crystal) = eachcol(crystal.atoms)

volume(cell::Matrix) = abs(det(cell))
volume(crystal::Crystal) = volume(crystal.cell)

end
