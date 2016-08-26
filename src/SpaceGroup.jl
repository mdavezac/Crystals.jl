module SpaceGroup
export cell_invariants, inner_translations

using Crystals.Constants: default_tolerance
using Crystals.Structure: Crystal

"""
gvectors with equivalent norms

Figures out vectors in the lattice defined by the cell which have the same
lengths as the column vectors defining the cell. these new vectors are
potentially symmetrically equivalent.
"""
function potential_equivalents(cell::Matrix; tolerance::Real=default_tolerance)
  const volume = abs(det(cell))
  const a0 = cell[:, 1]
  const a1 = cell[:, 2]
  const a2 = cell[:, 3]

  lengths = reducedim(+, cell .* cell, 1)
  max_norm = mapreduce(i -> norm(cell[:, i]), max, 0, 1:size(cell, 2))

  const n0 = Int(ceil(max_norm * norm(cross(a1, a2)) / volume))
  const n1 = Int(ceil(max_norm * norm(cross(a2, a0)) / volume))
  const n2 = Int(ceil(max_norm * norm(cross(a0, a1)) / volume))

  gvectors = Any[Array{eltype(cell), 1}[] for u in 1:length(lengths)]
  for i in -n0:n0, j in -n1:n1, k in -n2:n2
    g = cell * eltype(cell)[i, j, k]
    glength = sum(g .* g)
    for (length, result) in zip(lengths, gvectors)
      if abs(length - glength) < tolerance
        push!(result, g)
      end
    end
  end

  [hcat(gs...) for gs in gvectors]
end

"""
Finds and stores point group operations for a given lattice

A lattice is defined by a 3x3 matrix or cell.  Rotations are determined from
G-vector triplets with the same norm as the unit-cell vectors.

Implementation taken from ENUM_.

.. _ENUM: http://enum.sourceforge.net/
"""
function cell_invariants(cell::Matrix; tolerance::Real=default_tolerance)
  @assert size(cell, 1) == size(cell, 2)
  const ndims = size(cell, 1)

  result = Matrix{eltype(cell)}[]

  avecs, bvecs, cvecs = potential_equivalents(cell, tolerance=tolerance)

  const identity = eye(size(cell, 1))
  const invcell = inv(cell)
  for i in 1:size(avecs, 2), j in 1:size(bvecs, 2), k in 1:size(cvecs, 2)
    # (potential) rotation in cartesian coordinates
    rotation = hcat(avecs[:, i], bvecs[:, j], cvecs[:, k])
    # check operator is invertible
    abs(det(rotation)) >= tolerance || continue

    # rotation in fractional coordinates
    rotation = rotation * invcell
    any(abs(rotation - identity) .> tolerance) || continue

    # check matrix is a rotation
    all(abs(rotation * transpose(rotation) - identity) .< tolerance) || continue

    # check rotation not in list 
    index = findfirst(x -> all(abs(x - rotation) .< tolerance), result)
    index != 0 || push!(result, rotation)
  end
  result
end

""" Looks for internal translations """
function inner_translations(crystal::Crystal; tolerance::Real=default_tolerance):
  any(isna(crystal[:position])) && error("Some positions are not available")
  any(isna(crystal[:species])) && error("Some species are not available")

  cell = gruber(structure.cell)

  # find specie with minimum number of atoms
  species_count = by(crystal[[:species]], :species, d -> nrow(d))
  front = species_count[findmin(species_count[:x1])[2], :]

  positions = convert(
    Array, crystal[crystal[:species] == front[:species], :position])

  translations = []
  for site in eachrow(crystal)
    site[:species] ≠ front[:species] && continue

    translation = into_voronoi(site.pos - center, cell)
    all(abs(translation) .< tolerance) && continue

    is_mapping = true
    for mapping in eachrow(structure)
      pos = into_cell(mapping[:position] + translation, cell)
      found = false
      for mappee in eachrow(structure)
        found = mapping[:specie] == mappee[:specie] &&
                all(abs(mappee[:position] - pos) .< tolerance)
        found && break
      end
      (!found) && (is_mapping = false; break)
    end
    is_mapping && push!(translations, into_voronoi(translation, cell))
  end

  translations
end

"""
Computes space-group operations
"""
function space_group(cell::Matrix, positions::Matrix, types::Vector;
                     tolerance=default_tolerance)
  @assert size(cell, 1) == size(cell, 2)
  @assert size(positions, 1) == size(cell, 1) || length(positions) == 0
  @assert length(types) == size(positions, 2)

#   from numpy import dot, allclose, zeros
#   from numpy.linalg import inv
#   from . import gruber, Atom, into_voronoi, into_cell, is_primitive
#   from .. import error
#   if len(lattice) == 0:
#       raise error.ValueError("Empty lattice")

#   if not is_primitive(lattice, tolerance):
#       raise error.ValueError("Input lattice is not primitive")

#   # Finds minimum translation.
#   translation = lattice[0].pos
#   cell = gruber(lattice.cell, tolerance=tolerance)
#   invcell = inv(cell)

#   point_group = cell_invariants(lattice.cell)
#   assert len(point_group) > 0

#   centered = [Atom(into_cell(u.pos - translation, cell, invcell), u.type) for u in lattice]

#   # translations limited to those from one atom type to othe atom of same type
#   translations = [u.pos for u in centered if u.type == lattice[0].type]

#   result = []
#   for pg in point_group:
#       for trial in translations:
#           # Checks that this is a mapping of the lattice upon itself.
#           for unmapped in centered:
#               transpos = into_cell(dot(pg, unmapped.pos) + trial, cell, invcell)
#               for atom in centered:
#                   if atom.type != unmapped.type:
#                       continue
#                   if allclose(atom.pos, transpos, tolerance):
#                       break
#               # else executes only no atom is mapping of unmapped
#               else:
#                   break
#           # else executes only if all positions in structures have mapping
#           else:
#               transform = zeros((len(trial) + 1, len(trial)), dtype='float64', order='F')
#               transform[:3, :3] = pg
#               transform[3, :] = into_voronoi(
#                   trial - dot(pg, translation) + translation, cell, invcell)
#               result.append(transform)
#               # only one trial translation is possible, so break out of loop early
#               break
#   return result
end
end
