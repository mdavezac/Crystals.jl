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
Base.size(s::Structure) = (size(s.positions, 1), length(s))
function Base.size(s::Structure, i::Int)
  @assert i == 1 || i == 2
  return i == 1 ? size(s.positions, 1): length(s)
end
Base.ndims(s::Structure) = size(s.cell, 1)
Base.endof(s::Structure) = length(s)
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


function Base.getindex(s::Structure, col_ind::Real)
  return col_ind == endof(s) ? s.positions: s.properties[col_ind]
end
function Base.getindex(s::Structure, col_ind::Symbol)
  return col_ind == :position ? s.positions: s.properties[col_ind]
end

"""
    getindex{T <: ColumnIndex}(structure::Structure, indices::Union{Real, Symbol}[])

New structure with only the specified properties. The cell, scale, and
positions are always copied.
"""
function Base.getindex{T <: ColumnIndex}(s::Structure, col_inds::AbstractVector{T})
    result = Structure(s.cell, s.scale)
    result.positions = s.positions
    indices = filter(𝒾 -> 𝒾 ∉ [1 :position], col_inds)
    result.properties = s.properties[indices]
    result
end


export Structure

end # module
