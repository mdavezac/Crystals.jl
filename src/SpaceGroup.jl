module SpaceGroup
export point_group_operations, inner_translations, is_primitive, primitive

using Crystals.Constants: default_tolerance
using Crystals.Structure: Crystal, volume
using Crystals.Gruber: gruber
using Crystals.Utilities: into_voronoi, is_periodic, into_cell
using Crystals: Log
using AffineTransforms: AffineTransform
using DataFrames: isna, by, nrow, eachrow, DataFrame

"""
gvectors with equivalent norms

Figures out vectors in the lattice defined by the cell which have the same
lengths as the column vectors defining the cell. these new vectors are
potentially symmetrically equivalent.
"""
function potential_equivalents(cell::Matrix; tolerance::Real=default_tolerance)
    const V = volume(cell)
    const a0 = cell[:, 1]
    const a1 = cell[:, 2]
    const a2 = cell[:, 3]

    lengths = reducedim(+, cell .* cell, 1)
    max_norm = mapreduce(i -> norm(cell[:, i]), max, 0, 1:size(cell, 2))

    const n0 = ceil(Integer, max_norm * norm(cross(a1, a2)) / V)
    const n1 = ceil(Integer, max_norm * norm(cross(a2, a0)) / V)
    const n2 = ceil(Integer, max_norm * norm(cross(a0, a1)) / V)

    gvectors = Any[Array{eltype(cell), 1}[] for u in 1:length(lengths)]
    for i in -n0:n0, j in -n1:n1, k in -n2:n2
        g = cell * eltype(cell)[i, j, k]
        glength = sum(g .* g)
        for (length, result) in zip(lengths, gvectors)
            abs(length - glength) < tolerance && push!(result, g)
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
function point_group_operations(cell::Matrix; tolerance::Real=default_tolerance)
    @assert size(cell, 1) == size(cell, 2)
    const ndims = size(cell, 1)

    avecs, bvecs, cvecs = potential_equivalents(cell, tolerance=tolerance)
    result = Matrix{eltype(cell)}[eye(eltype(cell), size(cell, 1))]

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
        all(abs(rotation * transpose(rotation) - identity) .< tolerance) ||
            continue

        # check rotation not in list 
        index = findfirst(x -> all(abs(x - rotation) .< tolerance), result)
        index ≠ 0 || push!(result, rotation)
    end
    map(result) do x; AffineTransform(x, zeros(size(cell, 1))) end
end

""" Looks for internal translations within a crystal """
function inner_translations(crystal::Crystal; tolerance::Real=default_tolerance)
    any(isna(crystal[:position])) &&
        Log.error("Some positions are not available")
    any(isna(crystal[:species])) &&
        Log.error("Some species are not available")

    cell = gruber(crystal.cell)

    # find species with minimum number of atoms
    species_count = by(crystal.atoms, :species, d -> nrow(d))
    species = species_count[findmin(species_count[:x1])[2], :species]
    k = findfirst(crystal[:species], species)
    center = crystal[k, :position]

    translations = []
    for site in eachrow(crystal)
        site[:species] ≠ species && continue

        translation = into_voronoi(site[:position] - center, cell)
        all(abs(translation) .< tolerance) && continue

        is_mapping = true
        for mapping in eachrow(crystal)
            pos = mapping[:position] + translation
            found = false
            for mappee in eachrow(crystal)
                mapping[:species] ≠ mappee[:species] && continue
                found = is_periodic(pos, mappee[:position], crystal.cell;
                                    tolerance=tolerance) || continue
                found && break
            end
            (!found) && (is_mapping = false; break)
        end
        is_mapping && push!(translations, into_voronoi(translation, cell))
    end
    translations
end

""" True if the crystal structure cannot be reduced """
is_primitive(crystal::Crystal; tolerance=default_tolerance) =
    length(inner_translations(crystal; tolerance=tolerance)) == 0

"""
    primitive(crystal::Crystal; tolerance=default_tolerance)

    Computes the primitive cell of the input crystal.
    If the crystal is primitive, it is returned as is, otherwise a new crystal is
    returned.
"""
function primitive(crystal::Crystal; tolerance=default_tolerance)
    nrow(crystal) == 0 && return crystal

    cell = gruber(crystal.cell)
    trans = inner_translations(crystal, tolerance=tolerance)
    length(trans) == 0 && return crystal

    # adds original translations.
    push!(trans, cell[:, 1])
    push!(trans, cell[:, 2])
    push!(trans, cell[:, 3])

    # Looks for cell with smallest volume 
    new_cell = deepcopy(crystal.cell)
    V = volume(new_cell)
    for (i, first) in enumerate(trans), (j, second) in enumerate(trans)
        i == j && continue
        for (k, third) in enumerate(trans)
            (i == k ||  j == k) && continue
            trial = hcat(first, second, third)
            abs(det(trial) < 1e-12) && continue
            (abs(det(trial)) > V - 3.0 * tolerance) && continue

            if det(trial) < 0e0
                trial[:, 2], trial[:, 1] = second, third
            end
            det(trial) < 0e0 && Log.error("Negative volume")

            int_cell = inv(trial) * cell
            all(abs(int_cell - round(Integer, int_cell) .< 1e-8)) || continue

            new_cell = trial
            V = volume(trial)
        end
    end

    # Found the new cell with smallest volume (e.g. primivite)
    V < volume(crystal) - tolerance ||
        Log.error("Found translation but no primitive cell.")

    # now creates new lattice.
    result = Crystal(eltype(crystal.cell), gruber(new_cell), crystal.scale)
    columns = filter(names(crystal.atoms)) do x; x ∉ (:site_id, :cell_id) end
    result.atoms = crystal[1:0, columns]

    for site in eachrow(crystal)
        position = into_cell(site[:position], result.cell)
        k = findfirst(eachrow(result)) do atom
            site[:species] == atom[:species] &&
            all(abs(position - atom[:position]) .< tolerance)
        end
        if k == 0
            push!(result, site; no_new_properties=true)
            result[end, :position] = position
        end
    end

    nrow(crystal) % nrow(result) ≠ 0 &&
        Log.error("Nb of atoms in output not multiple of input.")

    abs(
        nrow(crystal) * volume(result) - nrow(result) * volume(crystal)
    ) < tolerance || Log.error("Size and volumes do not match.")
    result
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
