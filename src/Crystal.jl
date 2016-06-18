module Crystal

using DataFrames: DataFrame, NA, DataArray,ColumnIndex, index
import DataFrames: names, names!, nrow, ncol

typealias POSITIONS Array{Float64, 2}
typealias CELL POSITIONS

"""
Specify a crystal structure

A structure is constructed empty of atoms. Structure

    structure = Structure([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0], 2e0)
    @assert abs(structure.scale - 2e0) < 1e-8

    push!(structure, [0, 0, 0], 'Al'; charge=3)
"""
type Structure
  " 3x3 matrix specifying the translational symmetries of the structure "
  cell::CELL
  " nx3 matrix specifying the atomic positions "
  positions::POSITIONS
  " other atomic properties "
  properties::DataFrame
  " Scale factor for the cell and atomic positions "
  scale::eltype(POSITIONS)

  function Structure(cell=eye(3), scale=1e0; specie_type=UTF8String)
    @assert size(cell, 1) == size(cell, 2)
    new(cell, POSITIONS(size(cell, 1), 0), DataFrame(specie=specie_type[]), scale)
  end
end

Base.length(s::Structure) = length(s.properties) + 1
nrow(s::Structure) = size(s.positions, 1)
ncol(s::Structure) = length(s)
Base.size(s::Structure) = (nrow(s), ncol(s))
function Base.size(s::Structure, i::Integer)
  if i == 1
    nrow(s)
  elseif i == 2
    ncol(s)
  else
    throw(ArgumentError("Structures have only two dimensions"))
  end
end
Base.ndims(s::Structure) = 2
Base.endof(s::Structure) = length(s.properties) + 1
Base.names(s::Structure) = [names(s.properties); :position]

"""
    Base.push!(structure::Structure, position, specie; <keyword arguments>)

Add an atom to the structure. The `position` and associated `specie` are
required arguments. Other properties can be added using keyword arguments.

# Examples

```julia
push!(structure, [0, 0, 0], NA)
push!(structure, [1, 2, 3], "Na"; shell="fruity")
```

The first line above adds an atom at the origin with a unidentified/missing
atomic specie (`NA`). The second line adds a Sodium atom at position `[1 2 3]`,
with an additional property `shell`. If other atoms in the structure have not
specified this property, then their `shell` is set to missing (`NA`).
"""
function Base.push!(s::Structure, position, specie; kwargs...)
  s.positions = hcat(s.positions, position)
  # Add row to properties
  push!(kwargs, (:specie, specie))
  push!(s.properties, fill(NA, ncol(s.properties)))
  for (name, value) in kwargs
    if name ∉ names(s.properties)
      valtype = typeof(value) == typeof(NA) ? Any: typeof(value)
      s.properties[name] = DataArray(valtype, nrow(s.properties))
    end
    s.properties[max(end, 1), name] = value
  end
  s
end

function Base.copy(s::Structure)
    result = Structure(s.cell, s.scale)
    result.positions = copy(s.positions)
    result.properties = copy(s.properties)
    result
end

function Base.getindex(s::Structure, col_ind::ColumnIndex)
  return col_ind ∈ [endof(s), :position] ? s.positions: s.properties[col_ind]
end

"""
    getindex{T <: ColumnIndex}(structure::Structure, indices::T[])

If requesting positions, then a new Structure object is returned, including
other requested atomic properties, and the cell and scale. Otherwise, returns a
dataframe with only the requested properties.
"""
function Base.getindex{T <: ColumnIndex}(
        s::Structure, col_inds::AbstractVector{T})
  positions = filter(𝒾 -> 𝒾 ∈ [endof(s) :position], col_inds)
  indices = filter(𝒾 -> 𝒾 ∉ [endof(s) :position], col_inds)
  length(positions) < 2 || throw(ArgumentError("Requested positions more than once"))
  if length(positions) == 0
    return s.properties[indices]
  else
    result = Structure(s.cell, s.scale)
    result.positions = copy(s.positions)
    result.properties = s.properties[indices]
    result
  end
end
Base.getindex(structure::Structure, ::Colon) = copy(structure)
Base.getindex(structure::Structure, ::Colon, ::Colon) = copy(structure)

function Base.getindex{T <: Real}(s::Structure,
                                  row_ind::Union{T, AbstractVector{T}},
                                  col_ind::ColumnIndex)
  if col_ind ∈ [endof(s) :position]
    return s.positions[:, row_ind]
  end
  s[col_ind][row_ind]
end
function Base.getindex{T <: ColumnIndex}(
        s::Structure, row_ind::Real, col_inds::AbstractVector{T})
  positions = filter(𝒾 -> 𝒾 ∈ [endof(s) :position], col_inds)
  indices = filter(𝒾 -> 𝒾 ∉ [endof(s) :position], col_inds)
  length(positions) < 2 || throw(ArgumentError("Requested positions more than once"))
  if length(positions) == 0
    return s.properties[row_ind, indices]
  else
    result = Structure(s.cell, s.scale)
    result.positions = s.positions[:, row_ind:row_ind]
    result.properties = s.properties[row_ind, indices]
    result
  end
end
function Base.getindex{R <: Real, T <: ColumnIndex}(
          s::Structure, row_inds::AbstractVector{R}, col_inds::AbstractVector{T})
  positions = filter(𝒾 -> 𝒾 ∈ [endof(s) :position], col_inds)
  indices = filter(𝒾 -> 𝒾 ∉ [endof(s) :position], col_inds)
  length(positions) < 2 || throw(ArgumentError("Requested positions more than once"))
  if length(positions) == 0
    return s.properties[row_inds, indices]
  else
    result = Structure(s.cell, s.scale)
    result.positions = s.positions[:, row_inds]
    result.properties = s.properties[row_inds, indices]
    result
  end
end
function Base.getindex{T<:ColumnIndex}(s::Structure,
                                       ::Colon,
                                       col_inds::Union{T, AbstractVector{T}})
  s[col_inds]
end
Base.getindex(s::Structure, row_ind::Real, ::Colon) = s[[row_ind], :]
function Base.getindex{R<:Real}(s::Structure, row_inds::AbstractVector{R}, ::Colon)
  result = Structure(s.cell, s.scale)
  result.positions = s.positions[:, row_inds]
  result.properties = s.properties[row_inds, :]
  result
end

export Structure

end # module
