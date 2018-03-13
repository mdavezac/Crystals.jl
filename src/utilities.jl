module Utilities
export hart_forcade, to_fractional, to_cartesian, is_periodic, into_cell, origin_centered
export into_voronoi, supercell, cell_parameters

using Crystals.Constants: default_tolerance
using Crystals.Structures: Crystal
using Crystals.SNF: smith_normal_form
using DataFrames: nrow
using NamedTuples: @NT
using Unitful
using Unitful: Dimensions, NoUnits
using MicroLogging
using DocStringExtensions
using Base.Iterators: cycle, take

"""
    eldimension(u::Any)

Dimension of the input, if it has one. Arrays have a dimension if their element
type has a uniform dimension.
"""
eldimension{T <: AbstractArray}(u::Type{T}) = eldimension(eltype(T))
eldimension{T, D, U}(::Type{Quantity{T, D, U}}) = D
eldimension{T <: Number}(::Type{T}) = Dimensions{()}
eldimension(u::DataType) = Dimensions{()}
eldimension(u::Any) = eldimension(typeof(u))

is_unitful(a::AbstractArray) = is_unitful(eltype(a))
is_unitful{T <: Number}(a::Type{T}) = Val{:unitless}()
is_unitful{T, D, U}(a::Type{Quantity{T, D, U}}) = Val{:unitful}()

""" Tuple holding Hart-Forcade transform """
const HartForcadeTransform = @NT(transform::Matrix, quotient::Vector)

"""
$(SIGNATURES)

Computes the cyclic group of a supercell with respect to a lattice. It makes it
possible to identify the class of periodically equivalent cell that a given
position within the supercell belongs to. The function returns a named tuple with the
transform and the quotient.

# Examples

```jldoctest
using Crystals, Unitful
fcc = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"
supercell = [0 2 2; 0 -4 2; -1 0 -2]
ht = hart_forcade(fcc, fcc * supercell)

println(ht)
println("Positions in supercell:")
for index in CartesianRange((ht.quotient...))
    position = inv(ht.transform) * [index[u] for u in eachindex(ht.quotient)]
    println("- ", ustrip.(position), " (", unit(eltype(position)), ")")
end

# output

Hart-Forcade transform
- transform (nm^-1): [-1.0 -1.0 1.0; -1.0 1.0 1.0; -1.0 1.0 3.0]
- quotient: [1,2,6]

Positions in supercell:
- [-1.0,0.0,0.0] (nm)
- [-2.0,0.5,-0.5] (nm)
- [-0.5,0.0,0.5] (nm)
- [-1.5,0.5,0.0] (nm)
- [0.0,0.0,1.0] (nm)
- [-1.0,0.5,0.5] (nm)
- [0.5,0.0,1.5] (nm)
- [-0.5,0.5,1.0] (nm)
- [1.0,0.0,2.0] (nm)
- [0.0,0.5,1.5] (nm)
- [1.5,0.0,2.5] (nm)
- [0.5,0.5,2.0] (nm)
```
"""
function hart_forcade(lattice::AbstractMatrix, supercell::AbstractMatrix; digits::Integer=8)
    fractional = convert(Matrix{Int64}, round.(inv(lattice) * supercell, digits))

    snf, left, right = smith_normal_form(fractional)

    HartForcadeTransform(left * inv(lattice), diag(snf))
end

function Base.show(io::IO, ht::HartForcadeTransform)
    if unit(eltype(ht.transform)) ≠ NoUnits
        trans_unit = " ($(unit(eltype(ht.transform))))"
    else
        trans_unit = ""
    end
    println(io, "Hart-Forcade transform")
    println(io, "- transform", trans_unit, ": ", ustrip.(ht.transform))
    println(io, "- quotient: ", ustrip.(ht.quotient))
end

"""
     to_fractional(positions::AbstractArray, cell::AbstractMatrix)

Converts positions to fractional units (not Unitful). If the positions have no units, they
are already fractional units, and this operation is a no-op.
"""
function to_fractional(pos::AbstractArray, cell::AbstractMatrix)
    to_fractional(is_unitful(pos), pos, cell)
end
to_fractional(::Val{:unitless}, positions::AbstractArray, ::AbstractMatrix) = positions
function to_fractional(v::Val{:unitful}, positions::AbstractArray, cell::AbstractMatrix)
    to_fractional(v, is_unitful(cell), positions, cell)
end
function to_fractional(::Val{:unitful}, ::Val{:unitless}, ::AbstractArray, ::AbstractMatrix)
    error("The cell is unitless. Cannot determine how to make positions fractional.")
end
function to_fractional(::Val{:unitful}, v::Val{:unitful},
                       positions::AbstractArray, cell::AbstractMatrix)
    inv(cell) * positions
end

"""
    to_cartesian(positions::AbstractArray, cell::AbstractMatrix)

Converts positions to Cartesian units (Unitful). If the positions have no units, they are
already Cartesian units, and this operation is a no-op.
"""
function to_cartesian(pos::AbstractArray, cell::AbstractMatrix)
    to_cartesian(is_unitful(pos), pos, cell)
end
to_cartesian(::Val{:unitful}, positions::AbstractArray, ::AbstractMatrix) = positions
function to_cartesian(v::Val{:unitless}, positions::AbstractArray, cell::AbstractMatrix)
    to_cartesian(v, is_unitful(cell), positions, cell)
end
function to_cartesian(::Val{:unitless}, ::Val{:unitless}, ::AbstractArray, ::AbstractMatrix)
    error("The cell has no physical units. " *
          "Cannot determine how to make positions Cartesian.")
end
function to_cartesian(::Val{:unitless}, ::Val{:unitful},
                      positions::AbstractArray, cell::AbstractMatrix)
    cell * positions
end


"""
    to_same_kind(positions::AbstractArray, same::AbstractArray, cell::AbstractMatrix)

Converts positions to same units, Cartesian or fractional, as `same`. If possible, it
returns a reference to `positions` without copy. The function helps maintain _kind_
stability in other utility functions.
"""
function to_same_kind(positions::AbstractArray, same::AbstractArray, cell::AbstractMatrix)
    to_same_kind(is_unitful(same), positions, cell)
end
function to_same_kind(::Val{:unitful}, positions::AbstractArray, cell::AbstractMatrix)
    to_cartesian(positions, cell)
end
function to_same_kind(::Val{:unitless}, positions::AbstractArray, cell::AbstractMatrix)
    to_fractional(positions, cell)
end

"""
    is_periodic(a::AbstractVector, b::AbstractVector, cell::AbstractMatrix;
                tolerance::Real=$(default_tolerance))

True if the positions are one-to-one periodic with respect to the input cell.
"""
function is_periodic(a::AbstractVector, b::AbstractVector, cell::AbstractMatrix;
                     tolerance::Real=default_tolerance)
    const result = abs.(
        mod.(to_fractional(a, cell) .- to_fractional(b, cell) + 0.5, 1) - 0.5
    ) .< tolerance
    all(result)
end

"""
    is_periodic(a::AbstractMatrix, b::AbstractVector, cell::AbstractMatrix;
                tolerance::Real=$(default_tolerance))

Array of boolean describing whether positions in `a` are periodic with positions in `b`.
"""
function is_periodic{T, D, U}(a::AbstractMatrix,
                              b::AbstractVector,
                              cell::AbstractMatrix{Quantity{T, D, U}};
                              tolerance::Real=default_tolerance)
    const result = abs.(
        mod.(to_fractional(a, cell) .- to_fractional(b, cell) + 0.5, 1) - 0.5
    ) .< tolerance
    all(result, 1)
end
function is_periodic{T, D, U}(a::AbstractVector,
                              b::AbstractMatrix,
                              cell::AbstractMatrix{Quantity{T, D, U}};
                              tolerance::Real=default_tolerance)
    is_periodic(b, a, cell, tolerance=tolerance)
end
function is_periodic{T, D, U}(a::AbstractMatrix,
                              b::AbstractMatrix,
                              cell::AbstractMatrix{Quantity{T, D, U}};
                              tolerance::Real=default_tolerance)
    const result = abs.(
        mod.(to_fractional(a, cell) .- to_fractional(b, cell) + 0.5, 1) - 0.5
    ) .< tolerance
    all(result, 1)
end

"""
    $(SIGNATURES)

Folds periodic positions into cell
"""
function into_cell(pos::AbstractArray, cell::AbstractMatrix;
                   tolerance::Real=default_tolerance)
    frac = mod.(to_fractional(pos, cell), 1)
    if tolerance > 0
        for i in 1:length(frac)
            if abs(frac[i] - 1e0) < tolerance
                frac[i] = 0e0
            end
        end
    end
    to_same_kind(frac, pos, cell)
end


"""
    origin_centered(positions::AbstractArrays, cell::AbstractMatrix)

Folds positions back to origin, such that each fractional component ``x_f`` is between
``-0.5\\leq x_f < 0.5``. If the input is in Cartesian (fractional) coordinates, then the
result is also in Cartesian (fractional) coordinates.
"""
function origin_centered(pos::AbstractArray, cell::AbstractMatrix)
    to_same_kind(mod.(to_fractional(pos, cell) .+ 0.5, -1) .+ 0.5, pos, cell)
end

"""
    into_voronoi(positions::AbstractArray, cell::AbstractMatrix; extent::Integer=1)

Folds positions into first Brillouin zone of the input cell. Makes a well-meaning effort at
returning the periodic image with the smallest possible norm. It recenter the atoms around
the origin and then looks for the smallest periodic images within `-extent:extent` cells. If
the cell is quite pathological, then the result will not be within the Voronoi cell.
"""
function into_voronoi(pos::AbstractArray, cell::AbstractMatrix; extent::Integer=1)
    zcentered = to_cartesian(origin_centered(pos, cell), cell)
    result = deepcopy(zcentered)
    norms = [norm(zcentered[:, z]) for z in 1:size(zcentered, 2)]
    for n in eachindex(norms)
        for i = -extent:extent, j = -extent:extent, k = -extent:extent
            translation = cell * [i, j, k]
            position = zcentered[:, n] + translation
            d = norm(position)
            if d < norms[n]
                result[:, n] = position
                norms[n] = d
            end
        end
    end
    to_same_kind(result, pos, cell)
end

"""
    $(SIGNATURES)

Creates a supercell from an input lattice.

# Parameters

* `lattice::Crystal`: the original lattice
* `supercell::AbstractMatrix`: the cell of the supercell in Cartesian coordinates
* `site_id::Bool`: Whether to add/modify an atomic property indicating the index
  of the site in the original lattice
* `cell_id::Bool`: Whether to add/modify an atomic property indicating the index
  of the cell the site belongs to
"""
function supercell(lattice::Crystal, supercell::AbstractMatrix;
                   site_id::Bool=true, tolerance::Real=default_tolerance)
    nrow(lattice) == 0 && error("Lattice is empty")

    newcell = to_cartesian(supercell, lattice.cell)
    transform, quotient = hart_forcade(lattice.cell, newcell)
    itransform = inv(transform)

    siteids = collect(take(cycle(1:nrow(lattice)), prod(quotient) * nrow(lattice)))
    result = lattice[siteids]
    result.cell = newcell

    positions = result[:cartesian]
    for (i, index) in zip(1:nrow(lattice):nrow(result), CartesianRange((quotient...)))
        for j in i:i+nrow(lattice)-1
            positions[:, j] += itransform * [index[u] for u in eachindex(quotient)]
        end
    end
    result[:position] = into_cell(positions, result.cell)
    if site_id
        result[:site_id] = siteids
    end
    result
end

""" Tuple holding cell parameters """
const CellParameters = @NT(a, b, c, α, β, γ)

"""
    cell_parameters(a::Quantity, b::Quantity, c::Quantity,
                    α::Quantity=(π/2)u"rad", β::Quantity=(π/2)u"rad",
                    γₒ::Quantity=(π/2)u"rad")

Computes the cell matrix from the cell parameters [a, b, c, α, β, γ].
"""
function cell_parameters(a::Quantity, b::Quantity, c::Quantity,
                         α::Quantity=(π/2)u"rad", β::Quantity=(π/2)u"rad",
                         γₒ::Quantity=(π/2)u"rad")
    cx = cos(β)
    cy = (cos(α) - cos(β)cos(γₒ))/sin(γₒ)
    units = unit(a)
    bb = uconvert(unit(a), b)
    cc = uconvert(unit(a), c)
    z₀ = typeof(a)(0)
    [a    bb * cos(γₒ) cc * cx;
     z₀   bb * sin(γₒ) cc * cy;
     z₀   z₀           cc * √(1 - cx * cx - cy * cy)]
end
cell_parameters(cp::CellParameters) = cell_parameters(cp...)

"""
    cell_parameters(cell::AbstractMatrix)
    cell_parameters(cell::Crystal)

Parameters (a, b, c, α, β, γ) of the input cell returned in a named tuple.
"""
function cell_parameters(cell::AbstractMatrix)
    G = transpose(cell) * cell
    a, b, c = sqrt.(diag(G))
    α = acos(0.5(G[2, 3] + G[3, 2])/(c * b))u"rad"
    β = acos(0.5(G[3, 1] + G[1, 3])/(a * c))u"rad"
    γₒ = acos(0.5(G[1, 2] + G[2, 1])/(a * b))u"rad"
    CellParameters(a, b, c, α, β, γₒ)
end

end
