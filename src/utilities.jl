module Utilities
export hart_forcade, is_periodic, into_cell, origin_centered, into_voronoi,
       supercell, cell_parameters, cell_parameters°

using Crystals.Constants: default_tolerance
using Crystals.Structure: Position, Crystal
using Crystals.SNF: smith_normal_form
using Crystals: Log
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
    is_periodic(a::Matrix, b::Matrix, cell::Matrix; tolerance=default_tolerance)
    is_periodic(a::Union{Position, Vector}, b::Union{Position, Vector},
                cell::Matrix; tolerance=default_tolerance)

True if the positions are one-to-one periodic with respect to the input cell.
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
    into_voronoi(pos, cell::Matrix)

Folds vector into first Brillouin zone of the input cell.
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

"""
    supercell(lattice::Crystal, supercell::Matrix; site_id=true, cell_id=true)

Creates a supercell from an input lattice.

# Parameters
* `lattice::Crystal`: the original lattice
* `supercell::Matrix`: the cell of the supercell in cartesian coordinates
* `site_id::Bool`: Whether to add/modify an atomic property indicating the index
  of the site in the original lattice
* `cell_id::Bool`: Whether to add/modify an atomic property indicating the index
  of the cell the site belongs to
"""
function supercell(lattice::Crystal, supercell::Matrix;
                   site_id::Bool=true, cell_id::Bool=true)
    nrow(lattice) == 0 && Log.error("Lattice is empty")

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
    vcat(Crystal(eltype(lattice.cell), supercell, lattice.scale), all_atoms...)
end

"""
```julia
cell_parameters(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat,
                α::AbstractFloat=π/2, β::AbstractFloat=π/2,
                γ::AbstractFloat=π/2)
```
Computes cell from input cell parameters [a, b, c, α, β, γ].  α, β, γ are in
radian.
"""
function cell_parameters(a::AbstractFloat, b::AbstractFloat,
                         c::AbstractFloat, α::AbstractFloat=π/2,
						 β::AbstractFloat=π/2, γ::AbstractFloat=π/2)
    cx = cos(β)
    cy = (cos(α) - cos(β)cos(γ))/sin(γ)
    [a b * cos(γ) c * cx;
     0 b * sin(γ) c * cy;
     0 0          c * √(1 - cx * cx - cy * cy)]
end

"""
    cell_parameters(cell::Matrix)

Parameters [a, b, c, α, β, γ] of the input cell. α, β, γ are in radian.
"""
function cell_parameters(cell::Matrix)
    G = transpose(cell) * cell
    a, b, c = √diag(G)
    α = acos(0.5(G[2, 3] + G[3, 2])/(c * b))
    β = acos(0.5(G[3, 1] + G[1, 3])/(a * c))
    γ = acos(0.5(G[1, 2] + G[2, 1])/(a * b))
    a, b, c, α, β, γ
end

"""
    cell_parameters(cell::Crystal)

Parameters [a, b, c, α, β, γ] of the input cell. `α`, `β`, `γ` are in radian.
`a`, `b`, `c` include the scale factor `crystal.scale`.
"""
cell_parameters(crystal::Crystal) =
    cell_parameters(crystal.scale * crystal.cell)

""" Cell parameters with angles in degrees """
function cell_parameters°(x::Any)
    p = cell_parameters(x)
    p[1], p[2], p[3], rad2deg(p[4]), rad2deg(p[5]), rad2deg(p[6])
end

""" Cell from parameters with angles in degrees """
function cell_parameters°(a::AbstractFloat, b::AbstractFloat,
                          c::AbstractFloat, α::AbstractFloat=90,
                          β::AbstractFloat=90, γ::AbstractFloat=90)
    cell_parameters(a, b, c, deg2rad(α), deg2rad(β), deg2rad(γ))
end

end
