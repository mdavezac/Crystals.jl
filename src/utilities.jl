module Utilities
export hart_forcade, is_periodic, into_cell, origin_centered, into_voronoi,
       supercell, cell_parameters, cell_parameters°, underlying_dimension,
       to_fractional

using Crystals.Constants: default_tolerance
using Crystals.Structure: Crystal
using Crystals.Positions: Position, PositionArray, PositionDataArray
using Crystals.SNF: smith_normal_form
using Crystals: Log
using DataFrames: nrow, DataArray
using Unitful: Quantity
import Unitful: dimension, Dimensions

"""
    underlying_dimension(u::Any)

Dimension of the input, if it has one. Arrays have a dimension if their element
type has a uniform dimension.
"""
underlying_dimension{T <: AbstractArray}(u::Type{T}) =
    underlying_dimension(eltype(T))
underlying_dimension{T <: Position}(u::Type{T}) =
    underlying_dimension(eltype(T))
underlying_dimension{T, D, U}(::Type{Quantity{T, D, U}}) = D
underlying_dimension{T <: Number}(::Type{T}) = Dimensions{()}
underlying_dimension(u::DataType) = Dimensions{()}
underlying_dimension(u::Any) = underlying_dimension(typeof(u))

"""
    is_position_array(x::Union{Position, AbstractArray})

True if the input collapses to an array of positions, rather than just
positions. This function is necessary, since it can be difficult to distinguish
between a `Vector`, a `Position`, a `Matrix`, or an array of `Position`.
"""
is_position_array{T <: Position}(::Type{T}) = false
is_position_array{T <: AbstractArray}(::Type{T}) =
    ndims(T) == 2 || eltype(T) <: Position
is_position_array(x::Any) = is_position_array(typeof(x))

typealias PositionTypes Union{Position, Vector}

"""
    hart_forcade(lattice::Matrix, supercell::Matrix; digits::Integer=8)

Computes the cyclic group of a supercell with respect to a lattice. It makes it
possible to identify the class of periodically equivalent cell that a given
position within the supercell belongs to.

Returns the transform and the quotient.
"""
function hart_forcade(lattice::Matrix, supercell::Matrix; digits::Integer=8)
    fractional = convert(Matrix{Int64}, round(inv(lattice) * supercell, digits))

    snf, left, right = smith_normal_form(fractional)

    left * inv(lattice), diag(snf)
end

"""
     to_fractional{Q <: Quantity}(positions::Matrix, cell::Matrix{Q})

Converts positions to fractional units. If the positions have no units, they are
already fractional units, and this operation is a no-op.
"""
function to_fractional{Q <: Quantity}(
                positions::Union{Position, AbstractArray}, cell::Matrix{Q})
    if underlying_dimension(positions) !== Dimensions{()}
        const result = inv(cell) * positions
        @assert underlying_dimension(result) === Dimensions{()}
        return result
    end
    positions;
end

"""
    is_periodic(a::Matrix, b::Matrix, cell::Matrix;
                tolerance::Real=default_tolerance)
    is_periodic(a::Union{Position, Vector}, b::Union{Position, Vector},
                cell::Matrix; tolerance::Real=default_tolerance)

True if the positions are one-to-one periodic with respect to the input cell.
Returns a boolean if the input are two positions, and an array of booleans if
the input are arrays of positions.
"""
function is_periodic{Tc <: Quantity}(a::Union{Position, AbstractArray},
                                     b::Union{Position, AbstractArray},
                                     cell::Matrix{Tc};
                                     tolerance::Real=default_tolerance)
    const result = abs(
        mod(to_fractional(a, cell) .- to_fractional(b, cell) + 0.5, 1) - 0.5
    ) .< tolerance
    is_position_array(result) ? vec(all(result, 1)): all(result)
end

"""
    into_cell(pos::Union{Vector, Position}, cell::Matrix)

Folds periodic positions into cell
"""
function into_cell{Tc <: Quantity}(pos::PositionTypes, cell::Matrix{Tc};
                                   tolerance::Real=default_tolerance)
    frac = mod(to_fractional(pos, cell), 1)
    if tolerance > 0
        for i in 1:length(frac)
            if abs(frac[i] - 1e0) < tolerance
                frac[i] = 0e0
            end
        end
    end
    cell * frac
end


"""
    origin_centered(pos, cell::Matrix)

Folds column vector(s)/Position(s) back to origin
"""
origin_centered(pos::PositionTypes, cell::Matrix, invcell::Matrix) =
    cell * (mod(to_fractional(pos, cell) .+ 0.5, -1) .+ 0.5)

"""
    into_voronoi(pos, cell::Matrix)

Folds column vector(s)/Position(s) into first Brillouin zone of the input cell.
Returns the periodic image with the smallest possible norm.
"""
function into_voronoi(pos::PositionTypes, cell::Matrix, invcell::Matrix)
    zcentered = origin_centered(pos, cell, invcell)
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

for name in (:into_cell, :into_voronoi, :origin_centered)
    @eval begin
        function $name{Tc <: Quantity}(
                        positions::AbstractArray, cell::Matrix{Tc}; kwargs...)
            result = similar(positions)
            for i = 1:size(positions, 2)
                result[:, i] = $name(positions[:, i], cell, invcell; kwargs...)
            end
            result
        end
        function $name{T <: Number, N}(positions::PositionArray{T, N},
                                       cell::Matrix, invcell::Matrix;
                                       kwargs...)
            result = similar(positions)
            for i = 1:length(positions)
                result[i] = $name(positions[i], cell, invcell; kwargs...)
            end
            result
        end
        function $name{T <: Number, N}(positions::PositionDataArray{T, N},
                                       args...; kwargs...)
            DataArray($name(positions.data, args...), positions.na)
        end
    end
end

"""
    supercell(lattice::Crystal, supercell::Matrix;
              site_id::Bool=true, cell_id::Bool=true)

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
                   site_id::Bool=true, cell_id::Bool=true,
                   tolerance::Real=default_tolerance)
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


    const invcell = inv(supercell)
    all_atoms = Any[]
    for i = 1:quotient[1], j = 1:quotient[2], k = 1:quotient[3]
        for n in 1:nrow(atoms)
            atoms[n, :position] = into_cell(
                lattice[n, :position] + itransform * [i, j, k],
                supercell, invcell; tolerance=tolerance
            )
        end
        if cell_id
            atoms[:, :cell_id] = convert(Position, [i, j, k])
        end
        push!(all_atoms, deepcopy(atoms))
    end
    vcat(Crystal(eltype(lattice.cell), supercell), all_atoms...)
end

"""
    cell_parameters(a::Real, b::Real, c::Real,
                    α::Real=π/2, β::Real=π/2, γ::Real=π/2)

Computes cell from input cell parameters [a, b, c, α, β, γ].  α, β, γ are in
radian.
"""
function cell_parameters(a::Real, b::Real, c::Real,
                         α::Real=π/2, β::Real=π/2, γ::Real=π/2)
    cx = cos(β)
    cy = (cos(α) - cos(β)cos(γ))/sin(γ)
    [a b * cos(γ) c * cx;
     0 b * sin(γ) c * cy;
     0 0          c * √(1 - cx * cx - cy * cy)]
end

"""
    cell_parameters(cell::Matrix)

Parameters [a, b, c, α, β, γ] of the input cell. α, β, γ are in radian.
See also `cell_parameters°`.
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
"""
cell_parameters(crystal::Crystal) = cell_parameters(crystal.cell)

"""
    cell_parameters°(x::Matrix)

Cell parameters with angles in degrees. See also `cell_parameters`.
"""
function cell_parameters°(x::Matrix)
    p = cell_parameters(x)
    p[1], p[2], p[3], rad2deg(p[4]), rad2deg(p[5]), rad2deg(p[6])
end

"""
    cell_parameters°(a::Real, b::Real, c::Real,
                     α::Real=90, β::Real=90, γ::Real=90)

Cell from parameters with angles in degrees. See also `cell_parameters`.
"""
function cell_parameters°(a::Real, b::Real, c::Real,
                          α::Real=90, β::Real=90, γ::Real=90)
    cell_parameters(a, b, c, deg2rad(α), deg2rad(β), deg2rad(γ))
end

end
