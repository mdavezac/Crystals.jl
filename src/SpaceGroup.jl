module SpaceGroup
export point_group_operations, is_primitive, primitive
export space_group

using Unitful: ustrip, Quantity
using Crystals.Constants: default_tolerance
using Crystals.Structures: Crystal, volume
using Crystals.Gruber: gruber
using Crystals.Utilities: into_voronoi, is_periodic, into_cell, is_unitful, to_fractional
using Crystals: Log
using AffineTransforms: AffineTransform
using DataFrames: isna, by, nrow, eachrow, DataFrame, AbstractDataFrame, groupby, ncol

"""
    potential_equivalents(cell::AbstractMatrix; tolerance::Real=$(default_tolerance))

Figures out vectors in the lattice defined by the cell which have the same
lengths as the column vectors defining the cell. these new vectors are
potentially symmetrically equivalent.
"""
function potential_equivalents(cell::AbstractMatrix; tolerance::Real=default_tolerance)
    const V = volume(cell)
    const a0 = cell[:, 1]
    const a1 = cell[:, 2]
    const a2 = cell[:, 3]

    lengths = reducedim(+, cell .* cell, 1)
    max_norm = mapreduce(i -> norm(cell[:, i]), max, zero(eltype(cell)), 1:size(cell, 2))

    const n0 = ceil(Integer, max_norm * norm(cross(a1, a2)) / V)
    const n1 = ceil(Integer, max_norm * norm(cross(a2, a0)) / V)
    const n2 = ceil(Integer, max_norm * norm(cross(a0, a1)) / V)

    gvectors = Any[Array{eltype(cell), 1}[] for u in 1:length(lengths)]
    for i in -n0:n0, j in -n1:n1, k in -n2:n2
        g = cell * [i, j, k]
        glength = sum(g .* g)
        for (length, result) in zip(lengths, gvectors)
            ustrip(abs(length - glength)) < tolerance && push!(result, g)
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
function point_group_operations{T <: Number}(cell::AbstractMatrix{T};
                                             tolerance::Real=default_tolerance)
    @assert size(cell, 1) == size(cell, 2)
    const ndims = size(cell, 1)

    avecs, bvecs, cvecs = potential_equivalents(cell, tolerance=tolerance)
    result = [eye(cell)]

    const identity = eye(cell)
    const invcell = inv(cell)
    for i in 1:size(avecs, 2), j in 1:size(bvecs, 2), k in 1:size(cvecs, 2)
        # (potential) rotation in cartesian coordinates
        rotation = hcat(avecs[:, i], bvecs[:, j], cvecs[:, k])
        # check operator is invertible
        ustrip(volume(rotation)) ≥ tolerance || continue

        # rotation in fractional coordinates
        rotation = rotation * invcell
        any(ustrip.(abs.(rotation .- identity)) .> tolerance) || continue

        # check matrix is a rotation
        if !all(ustrip.(abs.(rotation * transpose(rotation) .- identity)) .< tolerance)
            continue
        end

        # check rotation not in list
        index = findfirst(x -> all(ustrip.(abs.(x .- rotation)) .< tolerance), result)
        index ≠ 0 || push!(result, rotation)
    end
    map(x -> AffineTransform(x, zeros(size(cell, 1))), result)
end

function point_group_operations{T, D, U}(cell::AbstractMatrix{Quantity{T, D, U}};
                                         tolerance::Real=default_tolerance)
    point_group_operations(ustrip(cell), tolerance=tolerance)
end

# """ Index to first instance of least represented species in crystal """
# function min_species_index(crystal::Crystal)
#     species_count = by(crystal.properties, :species, d -> nrow(d))
#     species = species_count[findmin(species_count[:x1])[2], :species]
#     k = findfirst(crystal[:species], species)
#     @assert k ≠ 0 # that should not be possible
#     k
# end

function inner_translations_impl(cell::AbstractMatrix,
                                 fractional::AbstractMatrix,
                                 species::AbstractVector;
                                 tolerance::Real=default_tolerance)
    @assert is_unitful(fractional) == Val{:unitless}()
    if length(species) ≠ size(fractional, 2)
        Log.error("Size of species and fractional positions do not match")
    end
    if size(cell, 2) ≠ size(fractional, 1)
        Log.error("Size of cell and fractional positions do not match")
    end
    grubcell = gruber(cell)

    # find species with minimum number of atoms
    unique_species = unique(species)
    k = findmin([countnz(species .== u) for u in unique(species)])[2]
    atom_type = unique_species[k]
    k = findfirst(species, atom_type)
    center = fractional[:, k]

    translations = []
    for site ∈ eachindex(species)
        species[site] ≠ atom_type && continue

        translation = into_voronoi(fractional[:, site] - center, grubcell)
        all(abs(translation) .< tolerance) && continue

        is_mapping = true
        for mapping ∈ eachindex(species)
            pos = fractional[:, mapping] + translation
            found = false
            for mappee ∈ eachindex(species)
                species[mappee] ≠ species[mapping] && continue
                found = is_periodic(pos, fractional[:, mappee], grubcell;
                                    tolerance=tolerance)
                found && break
            end
            (!found) && (is_mapping = false; break)
        end
        is_mapping && push!(translations, into_voronoi(translation, grubcell))
    end
    translations
end

""" Internal translations which leave a crystal unchanged """
function inner_translations(cell::AbstractMatrix,
                            positions::AbstractMatrix,
                            species::AbstractVector;
                            kwargs...)
    inner_translations_impl(cell, to_fractional(positions, cell), species; kwargs...)
end
function inner_translations(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}};
                            kwargs...)
    @assert nrow(crystal.properties) == nrow(crystal)
    inner_translations(crystal.cell,
                       crystal[:fractional],
                       species_ids(crystal.properties, col);
                       kwargs...)
end
function inner_translations(crystal::Crystal; kwargs...)
    if nrow(crystal.properties) == 0 || ncol(crystal.properties) == 0
        inner_translations(crystal.cell, crystal[:fractional], 1:nrow(crystal); kwargs...)
    else
        inner_translations(crystal, names(crystal.properties); kwargs...)
    end
end


"""
    species_ids(properties::AbstractDataFrame, cols::Union{Symbol, AbstractVector{Symbol}})

Returns an array where each unique set of elements in the input columns corresponds to a
unique element in the output. This functions helps creates an array corresponding to atomic
species, where a species can be defined by more than one input columns.
"""
function species_ids(properties::AbstractDataFrame, cols::AbstractVector{Symbol})
    if length(cols) == 0
        return Vector{Int64}(1:nrow(properties))
    end
    props = copy(properties)
    row_ids_name = Symbol(prod((string(u) for u in names(properties))))
    props[row_ids_name] = 1:nrow(props)

    results = zeros(Int64, (nrow(props,)))
    for (i, group) in enumerate(groupby(props, cols))
        results[group[row_ids_name]] = i
    end
    results
end

species_ids(properties::AbstractDataFrame, col::Symbol) = properties[col]

"""
    is_primitive(cell::AbstractMatrix, fractional::AbstractMatrix, species::AbstractVector;
                 tolerance::Real=$(default_tolerance))
    is_primitive(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}}; kwargs...)
    is_primitive(crystal::Crystal; kwargs...)

True if the crystal structure is primitice, e.g. not a supercell, e.g. not reducible to a
an equivalent lattice with fewer atoms.
"""
is_primitive(args...; kwargs...) = length(inner_translations(args...; kwargs...)) == 0

# function fewest_atoms(properties::AbstractDataFrame, unicity::AbstractVector{Symbol})
#     properties = copy(properties)
#     row_id_name = Symbol(prod((string(u) for u in names(properties))))
#     properties[row_id_names] = 1:nrow(properties)
# end


# """
#     primitive(crystal::Crystal; tolerance::Real=$(default_tolerance))

# Computes the primitive cell of the input crystal. If the crystal is primitive,
# it is returned as is, otherwise a new crystal is returned.
# """
# function primitive(crystal::Crystal; tolerance::Real=default_tolerance)
#     nrow(crystal) == 0 && return crystal

#     cell = gruber(crystal.cell)
#     trans = inner_translations(crystal, tolerance=tolerance)
#     length(trans) == 0 && return crystal

#     # adds original translations.
#     push!(trans, cell[:, 1])
#     push!(trans, cell[:, 2])
#     push!(trans, cell[:, 3])

#     # Looks for cell with smallest volume
#     new_cell = deepcopy(crystal.cell)
#     V = volume(new_cell)
#     for (i, first) in enumerate(trans), (j, second) in enumerate(trans)
#         i == j && continue
#         for (k, third) in enumerate(trans)
#             (i == k ||  j == k) && continue
#             trial = hcat(first, second, third)
#             abs(det(trial) < 1e-12) && continue
#             (abs(det(trial)) > V - 3.0 * tolerance) && continue

#             if det(trial) < 0e0
#                 trial[:, 2], trial[:, 1] = second, third
#             end
#             det(trial) < 0e0 && Log.error("Negative volume")

#             int_cell = inv(trial) * cell
#             all(abs(int_cell - round(Integer, int_cell) .< 1e-8)) || continue

#             new_cell = trial
#             V = volume(trial)
#         end
#     end

#     # Found the new cell with smallest volume (e.g. primivite)
#     V < volume(crystal) - tolerance ||
#         Log.error("Found translation but no primitive cell.")

#     # now creates new lattice.
#     result = Crystal(eltype(crystal.cell), gruber(new_cell), crystal.scale)
#     columns = filter(names(crystal.atoms)) do x; x ∉ (:site_id, :cell_id) end
#     result.atoms = crystal[1:0, columns]
#     Log.debug("Found gruber cell $(result.cell)")

#     invcell = inv(result.cell)
#     for site in eachrow(crystal)
#         position = into_cell(
#             site[:position], result.cell, invcell; tolerance=tolerance)
#         k = findfirst(eachrow(result)) do atom
#             site[:species] == atom[:species] &&
#             all(abs(position - atom[:position]) .< tolerance)
#         end
#         if k == 0
#             push!(result, site; no_new_properties=true)
#             result[end, :position] = position
#             Log.debug("Adding one atom", position=position)
#         end
#     end

#     nrow(crystal) % nrow(result) ≠ 0 && Log.error(
#         "Nb of atoms in output not multiple of input: " *
#         "$(nrow(crystal)) % $(nrow(result)) ≠ 0"
#     )

#     abs(
#         nrow(crystal) * volume(result) - nrow(result) * volume(crystal)
#     ) < tolerance || Log.error(
#         "Size and volumes do not match: " *
#         "abs($(nrow(crystal)) * $(volume(result)) - $(nrow(result)) *" *
#         "$(volume(crystal))) ≥ $(tolerance)"
#     )
#     result
# end

# """
# Computes space-group operations
# """
# function space_group(crystal::Crystal, tolerance::Real=default_tolerance)

#     const cell = gruber(crystal.cell, tolerance=tolerance)
#     const invcell = inv(cell)
#     const site = crystal[min_species_index(crystal), [:position, :species]]
#     const point_group = point_group_operations(cell, tolerance=tolerance)

#     # translations limited to those from one atom type to othe atom of same type
#     translations = crystal[crystal[:species] .== site[:species], :position]
#     translations -= translations[1]
#     translations = into_cell(translations, cell, invcell; tolerance=tolerance)

#     result = AffineTransform{eltype(crystal.cell), size(crystal.cell, 1)}[]
#     for pg in point_group
#         for trial in translations
#             is_invariant = findfirst(eachrow(crystal)) do mapper
#                 position = pg * (mapper[:position] - site[1, :position]) + trial
#                 mappee = findfirst(eachrow(crystal)) do atom
#                     mapper[:species] == atom[:species] &&
#                     is_periodic(position, atom[:position] - site[1, :position], cell, invcell;
#                                 tolerance=tolerance)
#                 end
#                 mappee == 0
#             end
#             if is_invariant == 0
#                 translation = into_voronoi(
#                     trial - pg * site[1, :position] + site[1, :position],
#                     cell, invcell
#                 )
#                 push!(result, pg * AffineTransform(eye(ndims(pg)), translation))
#             end
#         end
#     end
#     Log.info(
#         "$(length(result)) symmetry operations found, with " *
#         "$(count(result) do op; all(abs(op.offset) .< 1e-8) end) " *
#         "pure symmetries."
#     )
#     result
# end
end
