module Crystal

import Base: length, ndims, push!
using DataFrames: DataFrame, NA, DataArray, nrow, ncol

typealias POSITIONS Array{Float64, 2}
typealias CELL POSITIONS

" All the information needed to specify a crystal structure "
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

import Base: push!, getindex, length
length(s::Structure) = size(s.positions, 2)
ndims(s::Structure) = size(s.cell, 1)
function push!(s::Structure, position, specie; kwargs...)
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
