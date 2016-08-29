module Utilities
export hart_forcade, is_periodic, into_cell, origin_centered, into_voronoi,
        supercell

using Crystals.Constants: default_tolerance
using Crystals.Structure: Position, Crystal
using Crystals.SNF: smith_normal_form
using DataFrames: nrow
"""
Hart-Forcade transform

Computes the cyclic group of a supercell with respect to a lattice. It makes it
possible to identify the class of periodically equivalent cell that a given
position within the supercell belongs to.

Returns the transform and the quotient.
"""
function hart_forcade(lattice::Matrix, supercell::Matrix; digits=8)
    fractional = convert(Matrix{Int64}, round(inv(lattice) * supercell, digits))

    snf, left, right = smith_normal_form(fractional)

    left * inv(lattice), diag(snf)
end

"""
True if the two positions are periodic
  
Returns a boolean if the input are two positions, and an array of booleans if
the input are arrays of positions.
"""
is_periodic(a::Matrix, b::Matrix, cell::Matrix; tolerance=default_tolerance) =
  all(abs(origin_centered(a - b, cell)) .< tolerance, 1)

function is_periodic(a::Union{Position, Vector},
                     b::Union{Position, Vector},
                     cell::Matrix; tolerance=default_tolerance)
    all(abs(origin_centered(a - b, cell)) .< tolerance)
end

""" Folds periodic positions into cell """
into_cell(pos, cell::Matrix) = cell * mod(inv(cell) * pos, 1)

""" Folds vector back to origin """
origin_centered(pos, cell::Matrix) =
    cell * (mod(inv(cell) * pos .+ 0.5, -1) .+ 0.5)

"""
Folds vector into first Brillouin zone of the input cell

Returns the periodic image with the smallest possible norm.
"""
function into_voronoi(pos, cell::Matrix)
    zcentered = origin_centered(pos, cell)
    result = deepcopy(zcentered)
    norms = [norm(zcentered[:, i]) for i in 1:size(zcentered, 2)]
        for n in 1:length(norms)
            for i = -1:1, j = -1:1, k = -1:1
                translation = cell * [i, j, k]
                position = zcentered[:, n] + translation 
                d = norm(position)
                if d < norms[n]
                    result[:, n] = position
                    norms[n] = d
                end
            end
        end
        result
end

""" Creates a supercell of an input lattice """
function supercell(lattice::Crystal, supercell; site_id=true, cell_id=true)
    nrow(lattice) == 0 && error("Lattice is empty")

    transform, quotient = hart_forcade(lattice.cell, supercell)
    itransform = inv(transform)
    atoms = deepcopy(lattice.atoms)
    if site_id
        atoms[:site_id] = 1:nrow(atoms)
    end
    if cell_id
        atoms[:cell_id] = convert(Position, [0, 0, 0])
    end


    all_atoms = Any[]
    for i = 1:quotient[1], j = 1:quotient[2], k = 1:quotient[3]
        for n in 1:nrow(atoms)
            atoms[n, :position] = into_cell(
            lattice[n, :position] + itransform * [i, j, k], supercell)
        end
        if cell_id
            atoms[:, :cell_id] = convert(Position, [i, j, k])
        end
        push!(all_atoms, deepcopy(atoms))
    end
    vcat(Crystal(supercell, lattice.scale), all_atoms...)
end
end
