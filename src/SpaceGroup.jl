module SpaceGroup
export point_group_operations, is_primitive, primitive
export space_group

using Crystals.Constants: default_tolerance
using Crystals.Structures: Crystal, volume
using Crystals.Gruber: gruber
using Crystals.Utilities: into_voronoi, is_periodic, into_cell, is_unitful, to_fractional
using Crystals.Utilities: to_cartesian, to_same_kind
using Crystals: Log
using Unitful: ustrip, Quantity, unit
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

Implementation taken from [ENUM](http://enum.sourceforge.net/).
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

function inner_translations_impl(fractional::AbstractMatrix,
                                 cell::AbstractMatrix,
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
    length(translations) == 0 && return similar(fractional, (size(fractional, 1), 0))
    hcat(translations...)
end

""" Internal translations which leave a crystal unchanged """
function inner_translations(positions::AbstractMatrix,
                            cell::AbstractMatrix,
                            species::AbstractVector;
                            kwargs...)
    inner_translations_impl(to_fractional(positions, cell), cell, species; kwargs...)
end
function inner_translations(crystal::Crystal, cols::Union{Symbol, AbstractVector{Symbol}};
                            kwargs...)
    @assert nrow(crystal.properties) == nrow(crystal)
    inner_translations(crystal[:fractional],
                       crystal.cell,
                       species_ids(nrow(crystal), crystal.properties, cols);
                       kwargs...)
end
function inner_translations(crystal::Crystal; kwargs...)
    inner_translations(crystal, names(crystal.properties); kwargs...)
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

function species_ids(nrows::Integer,
                     properties::AbstractDataFrame,
                     cols::AbstractVector{Symbol})
    if length(cols) == 0
        Int64[1:nrows]
    else
        species_ids(properties, cols)
    end
end
species_ids(::Integer, props::AbstractDataFrame, cols::Symbol) = species_ids(props, cols)

"""
    is_primitive(cartesian::AbstractMatrix, cell::AbstractMatrix, species::AbstractVector;
                 tolerance::Real=$(default_tolerance))
    is_primitive(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}}; kwargs...)
    is_primitive(crystal::Crystal; kwargs...)

True if the crystal structure is primitice, e.g. not a supercell, e.g. not reducible to a
an equivalent lattice with fewer atoms.
"""
is_primitive(args...; kwargs...) = size(inner_translations(args...; kwargs...), 2) == 0

# function fewest_atoms(properties::AbstractDataFrame, unicity::AbstractVector{Symbol})
#     properties = copy(properties)
#     row_id_name = Symbol(prod((string(u) for u in names(properties))))
#     properties[row_id_names] = 1:nrow(properties)
# end
# """ Index to first instance of least represented species in crystal """
# function min_species_index(crystal::Crystal)
#     species_count = by(crystal.properties, :species, d -> nrow(d))
#     species = species_count[findmin(species_count[:x1])[2], :species]
#     k = findfirst(crystal[:species], species)
#     @assert k ≠ 0 # that should not be possible
#     k
# end



"""
    primitive(crystal::Crystal; tolerance::Real=$(default_tolerance))

Computes the primitive cell of the input crystal. If the crystal is primitive,
it is returned as is, otherwise a new crystal is returned.
"""
function primitive(crystal::Crystal; kwargs...)
    primitive(crystal, names(crystal.properties); kwargs...)
end
function primitive(crystal::Crystal, cols::Union{Symbol, AbstractVector{Symbol}}; kwargs...)
    @assert nrow(crystal.properties) == nrow(crystal)
    const species = species_ids(nrow(crystal), crystal.properties, cols)
    new_cell, indices = primitive(crystal.cell, crystal[:position], species; kwargs...)
    const incell = into_cell(crystal[indices, :cartesian], new_cell)
    const positions = to_same_kind(incell, incell, new_cell)
    typeof(crystal)(new_cell, positions, crystal.properties[indices, :])
end
function primitive(cell::AbstractMatrix, positions::AbstractMatrix, species::AbstractVector;
                   kwargs...)
    primitive_impl(cell, to_cartesian(positions, cell), species; kwargs...)
end

function primitive_impl(cell::AbstractMatrix,
                        cartesian::AbstractMatrix,
                        species::AbstractVector;
                        tolerance::Real=default_tolerance)
    if size(cartesian, 2) ≠ length(species)
        Log.error("Positions and species have incompatible sizes")
    end
    size(cell, 2) == size(cartesian, 1) || Log.error("Cell and positions are incompatible")
    size(cartesian, 2) == 0 && return cell, Int64[1:length(species)...]

    grubcell = gruber(cell)
    trans = grubcell * inner_translations(inv(grubcell) * cartesian, grubcell, species;
                                          tolerance=tolerance)
    size(trans, 2) == 0 && return cell, Int64[1:length(species)...]

    # adds original translations.
    trans = hcat(trans, grubcell)

    # Looks for cell with smallest volume
    new_cell = deepcopy(grubcell)
    V = volume(new_cell)
    for i ∈ 1:size(trans, 2), j ∈ 1:size(trans, 2), k ∈ 1:size(trans, 2)
        (i == k ||  j == k || i == j) && continue
        trial = trans[:, [i, j, k]]
        volume(trial) < 1e-12 * unit(V) && continue
        (volume(trial) > V - 3.0 * tolerance * unit(V)) && continue

        if det(ustrip(trial)) < 0e0
            trial = trans[:, [k, j, i]]
        end
        det(ustrip(trial)) < 0e0 && Log.error("Negative volume")

        int_cell = inv(trial) * cell
        all(abs(int_cell - round(Integer, int_cell) .< 1e-8)) || continue

        new_cell = trial
        V = volume(trial)
    end

    # Found the new cell with smallest volume (e.g. primivite)
    if V ≥ volume(grubcell) - tolerance * unit(V)
        Log.error("Found translation but no primitive cell.")
    end

    # now creates new lattice.
    indices = Int64[]
    Log.debug("Found gruber cell $(grubcell)")

    for site in 1:size(cartesian, 2)
        position = into_cell(cartesian[:, site], new_cell; tolerance=tolerance)
        k = findfirst(indices) do index
            species[site] == species[index] &&
            is_periodic(position, cartesian[:, index], new_cell, tolerance=tolerance)
        end
        if k == 0
            push!(indices, site)
        end
    end

    length(species) % length(indices) ≠ 0 && Log.error(
        "Nb of atoms in output not multiple of input: " *
        "$(length(species)) % $(length(indices)) ≠ 0"
    )

    abs(
        length(species) * volume(new_cell) - length(indices) * volume(cell)
       ) < tolerance * unit(volume(new_cell))|| Log.error(
        "Size and volumes do not match: " *
        "abs($(length(species)) * $(volume(new_cell)) - $(length(indices)) *" *
        "$(volume(cell))) ≥ $(tolerance * unit(volume(new_cell)))"
    )
    new_cell, indices
end

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
