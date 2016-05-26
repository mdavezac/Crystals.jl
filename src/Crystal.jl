module Crystal

using DataFrames: DataFrame

typealias POSITIONS Array{Float64, 2}
typealias CELL POSITIONS
typealias SPECIES Array{UTF16String, 1}

" All the information needed to specify a crystal structure "
type Structure
  " 3x3 matrix specifying the translational symmetries of the structure "
  cell::CELL
  " nx3 matrix specifying the atomic positions "
  positions::POSITIONS
  " array of atomic species "
  species::SPECIES
  " other atomic properties "
  properties::DataFrame
  " Scale factor for the cell and atomic positions "
  scale::eltype(POSITIONS)

  function Structure(cell=eye(3), positions=POSITIONS(), specie=SPECIES(),
                     attributes=DataFrame(), scale=1e0)
    new(eye(3), POSITIONS(), SPECIES(), DataFrame(), 1e0)
  end
end

import Base: push!, getindex, length
length(s::Structure) = length(positions)
getindex(s::Structure, i) = (s.positions[i], s.species[i], s.properties[i, :])
function push!(s::Structure, position, specie; kwargs...)
  if size(s.positions) == 0
    s.positions = position
    return
  end
  s.positions = vcat(s.positions, position)
  # s.species = vcat(s.species, [species])
end


export Structure

end # module
