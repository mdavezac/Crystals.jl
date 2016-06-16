module Crystal

using DataFrames: DataFrame, NA, DataArray, nrow, ncol

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

Base.length(s::Structure) = size(s.positions, 2)
Base.ndims(s::Structure) = size(s.cell, 1)

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


export Structure

end # module
