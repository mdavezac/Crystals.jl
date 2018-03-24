module Structures
using MicroLogging
using DocStringExtensions
export AbstractCrystal, Crystal, is_fractional, volume, are_compatible_lattices, round!
using Unitful: Quantity, Dimensions, Units, unit, ustrip
using Missings: Missing, missing
using ArgCheck: @argcheck
using Base.Iterators: filter, drop

using DataFrames: DataFrame, nrow, missing, index, ncol, deleterows!, eltypes
using Crystals.Constants: default_tolerance
import Base
import Base: delete!
import Unitful
import DataFrames
using Unitful: dimension, unit

""" Reserved column names for special input and output """
const RESERVED_COLUMNS = [:position, :fractional, :cartesian, :x, :y, :z]

""" Top type node for Crystals """
abstract type AbstractCrystal end

RowIndices{T <: Integer} = Union{T, AbstractVector{T}, Range{T}, Colon}

""" Describe a crystalline structure """
type Crystal{T <: Number, D, U, P <: Number} <: AbstractCrystal
    """ Periodicity of the crystal structure """
    cell::Matrix{Quantity{T, D, U}}
    """ Atomic positions """
    positions::Matrix{P}
    """ Atomic properties """
    properties::DataFrame
end

function Crystal{C, D, U, P <: Number}(cell::Matrix{Quantity{C, D, U}},
                                       positions::Matrix{P},
                                       args...; kwargs...)
    @argcheck size(cell, 1) == size(cell, 2)
    @argcheck size(positions, 1) == size(cell, 1)

    properties = DataFrame(args...; kwargs...)

    @argcheck length(names(properties) ∩ RESERVED_COLUMNS) == 0
    @argcheck nrow(properties) == 0 || nrow(properties) == size(positions, 2)

    if P <: Quantity
        const Q = promote_type(Quantity{C, D, U}, P)
        p = convert(Matrix{Q}, positions)
        Crystal{Q.parameters..., Q}(convert(Matrix{Q}, cell),
                                    convert(Matrix{Q}, positions),
                                    properties)
    else
        const T = promote_type(C, P)
        const Q = Quantity{T, D, U}
        Crystal{T, D, U, T}(convert(Matrix{Q}, cell),
                            convert(Matrix{T}, positions),
                            properties)
    end
end


function Crystal{C, D, U}(cell::Matrix{Quantity{C, D, U}}; kwargs...)
    const dpositions = collect(filter(x -> x[1] ∈ (:position, :positions), kwargs))
    const tpositions = collect(filter(x -> x[1] ∈ (:tposition, :tpositions), kwargs))
    if length(dpositions) ≠ 0 && length(tpositions) ≠ 0
        const positions = hcat((x[2] for x in dpositions)...,
                               transpose(hcat([x[2] for x in tpositions]...)))
    elseif length(dpositions) ≠ 0
        const positions = hcat([x[2] for x in dpositions]...)
    elseif length(tpositions) ≠ 0
        const positions = transpose(hcat([x[2] for x in tpositions]...))
    else
        const positions = Matrix{Quantity{C, D, U}}(size(cell, 1), 0)
    end

    leftover = filter(x -> x[1] ∉ (:position, :positions, :tposition, :tpositions), kwargs)
    Crystal(cell, positions; leftover...)
end

""" Underlying physical units of the crystal """
Unitful.unit(crystal::Crystal) = typeof(crystal).parameters[3]()
Unitful.dimension(crystal::Crystal) = typeof(crystal).parameters[2]()
Base.length(crystal::Crystal) = size(crystal.positions, 2)
"""
    $(SIGNATURES)

True if the crystal structure is fractional.
"""
is_fractional(crystal::Crystal) = is_fractional(typeof(crystal))
is_fractional{T <: Crystal}(::Type{T}) = !(T.parameters[end] <: Quantity)
fractional_trait(crystal::Crystal) = fractional_trait(typeof(crystal))
function fractional_trait{T <: Crystal}(::Type{T})
    is_fractional(T) ? Val{:fractional}(): Val{:cartesian}()
end

volume(cell::Matrix) = abs(det(cell))
function volume{T, D, U}(cell::Matrix{Quantity{T, D, U}})
    abs(det(ustrip(cell))) * unit(eltype(cell))^3
end
"""
    $(SIGNATURES)

Returns the volume of a `Crystal` instance or of a cell. It comes down to computing
``|\det(A)|`` where ``A`` is the crystal cell.

# Examples

```jldoctest
julia> using Crystals

julia> volume(Crystal([0 2 2; 2 0 2; 2 2 0]u"nm"))
16.0 nm^3

julia> volume([0 2 2; 2 0 2; 2 2 0]u"nm")
16.0 nm^3

julia> volume([0 2 2; 2 0 2; 2 2 0])
16.0
```
"""
volume(crystal::Crystal) = volume(crystal.cell)

function are_compatible_lattices(cell_a::Matrix, cell_b::Matrix;
                                 digits=12, rtol=default_tolerance, atol=default_tolerance)
    all(isinteger.(round.(inv(cell_a) * cell_b, digits))) &&
        isapprox(ustrip(volume(cell_a)), ustrip(volume(cell_b)); rtol=rtol, atol=atol)
end
function are_compatible_lattices(crystal::Crystal, cell::Matrix; kwargs...)
    are_compatible_lattices(crystal.cell, cell; kwargs...)
end
function are_compatible_lattices(crystal_a::Crystal, crystal_b::Crystal; kwargs...)
    are_compatible_lattices(crystal_a.cell, crystal_b.cell; kwargs...)
end
function are_compatible_lattices(cell::Matrix, crystal::Crystal; kwargs...)
    are_compatible_lattices(cell, crystal.cell; kwargs...)
end
"""
    are_compatible_lattices(lattices...)

True if the lattices are mathematically equivalent. Two lattices are equivalent if they
represent the same periodicity. In practice, this means the two lattices have the same
volume, and their cell vectors are integer linear combinations of one another.

# Parameters

* `lattices::Vararg{Union{Matrix, Crystal}}`: Any number of lattices
* `digits::Integer`: when checking the cells are integer co-linear, the product ``A^{-1}B``
  is first rounded to this number of digits
* `rtol::Real`: relative tolerance when checking the volumes correspond (no units).
  Default to $(default_tolerance).
* `atol::Real`: absolute tolerance when checking the volumes correspond (no units).
  Default to $(default_tolerance).

# Examples

```jldoctest
using Crystals, Unitful
crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm")
cells = Matrix{Int64}[eye(3)]
while length(cells) < 5
    cell = rand(-5:5, (3, 3))
    volume(cell) == 1 && push!(cells, cell)
end
lattices = [crystal.cell * c for c in cells]
are_compatible_lattices(crystal, lattices...)

# output
true
```
"""
function are_compatible_lattices(lattices::Vararg{Union{Matrix, Crystal}})
    if length(lattices) < 2
        return true
    end
    all(are_compatible_lattices(lattices[1], u) for u in lattices[2:end])
end


"""
    position_for_crystal(crystal::Crystal, position::AbstractArray)
    position_for_crystal(crystal::Crystal, compatible_crystal::Crystal)

Converts position to same type as in the input crystal. If real positions are
given for a crystal with fractional positions, then the positions are converted.
And vice-versa.
"""
function position_for_crystal(crystal::Crystal, position::AbstractArray)
    position_for_crystal(fractional_trait(crystal), crystal.cell, position)
end
function position_for_crystal{T <: Quantity}(::Val{:fractional},
                                             cell::AbstractMatrix,
                                             position::AbstractArray{T})
    inv(cell) * position
end
position_for_crystal(::Val{:fractional}, ::AbstractMatrix, pos::AbstractArray) = pos
position_for_crystal(::Val{:cartesian}, cell::AbstractMatrix, p::AbstractArray) = cell * p
function position_for_crystal{T <: Quantity}(::Val{:cartesian},
                                             ::AbstractMatrix,
                                             positions::AbstractArray{T})
    positions
end
position_for_crystal(a::Val, c::Crystal) = position_for_crystal(a, c.cell, c.positions)
function position_for_crystal(crystal::Crystal, other::Crystal)
    @assert are_compatible_lattices(crystal, other)
    position_for_crystal(crystal, position_for_crystal(Val{:cartesian}(), other))
end

@inline _is_not_position(s::Symbol) = :position ≠ s
@inline _is_not_position(s::AbstractVector{Symbol}) = :position ∉ s

function Base.show(io::IO, crystal::Crystal)
    println(io, "cell(", unit(crystal), "):")
    for i in 1:size(crystal.cell, 1)
        print(io, "  ")
        join(io, ustrip.(crystal.cell[i, :]), ' ')
        println(io)
    end
    with_pos = copy(crystal.properties)
    const name = is_fractional(crystal) ? :fractional : :Cartesian
    const positions = [tuple(ustrip.(crystal.positions[:, i])...) for i in 1:length(crystal)]
    show(io, hcat(DataFrame(Any[positions], [name]), crystal.properties),
         false, :Atom, false)
end

"""
    $(SIGNATURES)

Appends an atom to a crystal structure. The position of the atom is a necessary
argument, whether in Cartesian or in fractional coordinates. If keyword arguments are
present, then they represent atomic properties for the atom being added. Properties that
are not explicitly given are set to `missing`. Similarly, new properties that were not
present in the crystal structure previously are `missing` except for the newly added atom.

# Examples

```jldoctest
using Crystals
crystal = Crystal(eye(2)u"km",
                  tpositions=[1 1; 2 3; 4 5]u"m",
                  species=["Al", "O", "O"],
                  label=[:+, :-, :-])
push!(crystal, [[10, 20]u"nm", "B", :a])
crystal

# output

cell(m):
  1000.0 0.0
  0.0 1000.0

│ Atom │ Cartesian        │ species │ label │
├──────┼──────────────────┼─────────┼───────┤
│ 1    │ (1.0, 1.0)       │ Al      │ +     │
│ 2    │ (2.0, 3.0)       │ O       │ -     │
│ 3    │ (4.0, 5.0)       │ O       │ -     │
│ 4    │ (1.0e-8, 2.0e-8) │ B       │ a     │
```
"""
function Base.push!(crystal::Crystal, associative::Associative{Symbol, <: Any})
    position = position_for_crystal(crystal, associative[:position])
    crystal.positions = hcat(crystal.positions, position)
    push!(crystal.properties, associative)
    crystal
end

function Base.push!(crystal::Crystal, associative::Associative)
    position_in = get(() -> associative["position"], associative, :position)
    position = position_for_crystal(crystal, position_in)
    push!(crystal.properties, associative)
    crystal.positions = hcat(crystal.positions, position)
    crystal
end

function Base.push!(crystal::Crystal, iterable::Any)
    position_in = first(iterable)
    position = position_for_crystal(crystal, position_in)
    push!(crystal.properties, drop(iterable, 1))
    crystal.positions = hcat(crystal.positions, position)
    crystal
end

function Base.getindex(crystal::Crystal, index::RowIndices)
    typeof(crystal)(crystal.cell, hcat(crystal.positions[:, index]),
                    crystal.properties[index, :])
end

function Base.getindex(crystal::Crystal, symbol::Symbol)
    if symbol == :position
        return crystal.positions
    elseif symbol == :cartesian
        return position_for_crystal(Val{:cartesian}(), crystal)
    elseif symbol == :fractional
        return position_for_crystal(Val{:fractional}(), crystal)
    elseif symbol == :x
        return crystal.positions[1, :]
    elseif symbol == :y
        return crystal.positions[2, :]
    elseif symbol == :z
        return crystal.positions[3, :]
    end
    crystal.properties[symbol]
end

Base.getindex(crystal::Crystal, ::Colon) = deepcopy(crystal)

function Base.getindex(crystal::Crystal, symbols::AbstractVector{Symbol})
    const specials = symbols ∩ RESERVED_COLUMNS
    if length(specials) == 0
        return crystal.properties[symbols]
    elseif length(specials) > 1 && length(setdiff(specials, (:x, :y, :z))) == 0
        args = collect(filter(x -> x ∉ specials, symbols))
        result = crystal.properties[args]
        for u in specials
            result[u] = crystal[u]
        end
        return result
    elseif length(specials) > 1
        allowed = setdiff(RESERVED_COLUMNS, (:x, :y, :z))
        error("Cannot use more than one of $allowed at a time")
    elseif specials[1] == :position
        args = collect(filter(x -> x ≠ :position, symbols))
        typeof(crystal)(crystal.cell, crystal.positions, crystal.properties[args])
    elseif specials[1] == :cartesian || specials[1] == :fractional
        const T, D, U = typeof(crystal).parameters[1:3]
        const positions = position_for_crystal(Val{specials[1]}(), crystal)
        const args = collect(filter(x -> x ≠ specials[1], symbols))
        const P = eltype(positions)
        Crystal{T, D, U, P}(crystal.cell, positions, crystal.properties[args])
    end
end

function Base.getindex(crystal::Crystal, index::RowIndices, symbol::Symbol)
    if symbol == :position
        return crystal.positions[:, index]
    elseif symbol == :x
        return crystal.positions[1, index]
    elseif symbol == :y
        return crystal.positions[2, index]
    elseif symbol == :z
        return crystal.positions[3, index]
    elseif symbol == :cartesian
        return position_for_crystal(Val{:cartesian}(),
                                    crystal.cell,
                                    crystal.positions[:, index])
    elseif symbol == :fractional
        return position_for_crystal(Val{:fractional}(),
                                    crystal.cell,
                                    crystal.positions[:, index])
    end
    crystal.properties[index, symbol]
end

function Base.getindex(crystal::Crystal, ::Colon, symbol::Symbol)
    if symbol == :position
        return crystal.positions
    elseif symbol == :x
        return crystal.positions[1, :]
    elseif symbol == :y
        return crystal.positions[2, :]
    elseif symbol == :z
        return crystal.positions[3, :]
    elseif symbol == :cartesian
        return position_for_crystal(Val{:cartesian}(), crystal)
    elseif symbol == :fractional
        return position_for_crystal(Val{:fractional}(), crystal)
    end
    crystal.properties[:, symbol]
end

Base.getindex(crystal::Crystal, ::Colon, ::Colon) = copy(crystal)
Base.getindex(crystal::Crystal, row::Any, ::Colon) = Base.getindex(crystal, row)

function Base.getindex(crystal::Crystal, index::RowIndices, symbols::AbstractVector{Symbol})
    const specials = symbols ∩ RESERVED_COLUMNS
    if length(specials) == 0
        return crystal.properties[index, symbols]
    elseif length(specials) > 1 && length(setdiff(specials, (:x, :y, :z))) == 0
        args = collect(filter(x -> x ∉ specials, symbols))
        result = crystal.properties[index, args]
        for u in specials
            result[u] = crystal[index, u]
        end
        return result
    elseif length(specials) > 1
        allowed = setdiff(RESERVED_COLUMNS, (:x, :y, :z))
        error("Cannot use more than one of $allowed at a time")
    elseif specials[1] == :position
        args = collect(filter(x -> x ≠ :position, symbols))
        typeof(crystal)(crystal.cell,
                        hcat(crystal.positions[:, index]),
                        crystal.properties[index, args])
    else
        const T, D, U = typeof(crystal).parameters[1:3]
        const positions = position_for_crystal(Val{specials[1]}(),
                                               crystal.cell,
                                               crystal.positions[:, index])
        const args = collect(filter(x -> x ≠ specials[1], symbols))
        const P = eltype(positions)
        Crystal{T, D, U, P}(crystal.cell, hcat(positions), crystal.properties[index, args])
    end
end

function Base.getindex(crystal::Crystal, symbol::Symbol, pos::Any)
    @assert symbol == :position
    crystal.positions[pos, :]
end

function Base.getindex(crystal::Crystal, atoms::Any, symbol::Symbol, pos::Any)
    @assert symbol == :position
    crystal.positions[pos, atoms]
end

DataFrames.eltypes(crystal::Crystal) = push!(eltypes(crystal.properties),
                                       Vector{eltype(crystal.positions)})
Base.names(crystal::Crystal) = push!(names(crystal.properties), :position)
Base.size(crystal::Crystal) = (r = size(crystal.properties); (r[1], r[2] + 1))
Base.size(crystal::Crystal, i::Integer) = size(crystal)[i]
Base.ndims(::Crystal) = 2
DataFrames.nrow(crystal::Crystal) = size(crystal.positions, 2)
DataFrames.ncol(crystal::Crystal) = ncol(crystal.properties) + 1
Base.endof(crystal::Crystal) = nrow(crystal)

function Base.copy(crystal::Crystal)
    typeof(crystal)(copy(crystal.cell), copy(crystal.positions), copy(crystal.properties))
end

function Base.deepcopy(crystal::Crystal)
    typeof(crystal)(deepcopy(crystal.cell),
                    deepcopy(crystal.positions),
                    deepcopy(crystal.properties))
end

function Base.setindex!(crystal::Crystal, v::DataFrame, cols::AbstractVector{Symbol})
    if :position ∈ names(v) || :position ∈ cols
        error("Positions cannot be set from a dataframe")
    end
    setindex!(crystal.properties, v, cols)
end

function Base.setindex!(crystal::Crystal,
                        v::DataFrame,
                        rows::Any,
                        cols::AbstractVector{Symbol})
    if :position ∈ names(v) || :position ∈ cols
        error("Positions cannot be set from a dataframe")
    end
    setindex!(crystal.properties, v, rows, cols)
end

function Base.setindex!(crystal::Crystal, v::Crystal, cols::AbstractVector{Symbol})
    if :position ∈ cols
        crystal.positions = position_for_crystal(crystal, v)
        setindex!(crystal.properties, v.properties,
                  collect(filter(x -> x ≠ :position, cols)))
    else
        setindex!(crystal.properties, v.properties, cols)
    end
end

function Base.setindex!(crystal::Crystal, v::Any, row::Any, col::Symbol)
    if col == :position
        setindex!(crystal.positions, position_for_crystal(crystal, v), :, row)
    elseif col == :x
        setindex!(crystal.positions, v, 1, row)
    elseif col == :y
        setindex!(crystal.positions, v, 2, row)
    elseif col == :z
        setindex!(crystal.positions, v, 3, row)
    else
        setindex!(crystal.properties, v, row, col)
    end
end

function Base.setindex!(crystal::Crystal, v::Any, col::Symbol)
    if col == :position
        if size(v) ≠ size(crystal.positions)
            error("Input has incorrect size")
        end
        crystal.positions = position_for_crystal(crystal, v)
    elseif col == :x
        if size(v) ≠ size(crystal.positions)
            error("Input has incorrect size")
        end
        crystal.positions[1, :] = position_for_crystal(crystal, v)
    elseif col == :y
        if size(v) ≠ size(crystal.positions)
            error("Input has incorrect size")
        end
        crystal.positions[2, :] = position_for_crystal(crystal, v)
    elseif col == :z
        if size(v) ≠ size(crystal.positions)
            error("Input has incorrect size")
        end
        crystal.positions[3, :] = position_for_crystal(crystal, v)
    else
        setindex!(crystal.properties, v, col)
    end
end

function Base.setindex!(crystal::Crystal, v::Crystal,
                        row::Any, cols::AbstractVector{Symbol}; check_periodicity=true)
    if :position ∈ cols
        if check_periodicity && !are_compatible_lattices(crystal, v)
            error("The two crystal structures do not have the same periodicity")
        end
        crystal.positions[:, row] = position_for_crystal(crystal, v)
        setindex!(crystal.properties, v.properties, row,
                  collect(filter(x -> x ≠ :position, cols)))
    else
        setindex!(crystal.properties, v.properties, row, cols)
    end
end

function Base.setindex!(crystal::Crystal, v::Any, col::Symbol, xyz::RowIndices)
    if col != :position
        error("crystal[:positions, integer] is only available for positions")
    end

    crystal.positions[xyz, :] = v
end

function Base.setindex!(crystal::Crystal, v::Any, rows::RowIndices,
                        col::Symbol, xyz::RowIndices)
    if col != :position
        error("crystal[:positions, integer] is only available for positions")
    end
    crystal.positions[xyz, rows] = v
end

Base.setindex!(crys::Crystal, v::Any, ::Colon, col::Symbol) = Base.setindex!(crys, v, col)

"""
    Base.delete!(crystal::Crystal, col::Symbol)
    Base.delete!(crystal::Crystal, col::AbstractVector{Symbol})

Deletes one or more atomic property. Positions cannot be deleted.
"""
function Base.delete!(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}})
    if !_is_not_position(col)
        error("Cannot delete position column from a structure")
    end
    delete!(crystal.properties, col)
    crystal
end

"""
    Base.delete!(crystal::Crystal, ::Colon)

Deletes all atomic properties except for positions.
"""
Base.delete!(crystal::Crystal, ::Colon) = (empty!(crystal.properties); crystal)

Base.setindex!(crystal::Crystal, ::Void, col::Any) = delete!(crystal, col)

"""
    $(SIGNATURES)

Deletes all atomic sites, both properties and positions.
"""
function Base.empty!(crystal::Crystal)
    empty!(crystal.properties)
    crystal.positions = crystal.positions[:, 1:0]
end

"""
    deleterows!(crystal::Crystal, rows::Integer)
    deleterows!(crystal::Crystal, rows::AbstractVector{Integer})
    deleterows!{T <: Integer}(crystal::Crystal, rows::Range{T})
    deleterows!(crystal::Crystal, rows::Colon)

Deletes one (single integer), a few (sequence or range), or all (colon) atoms in the
structure.
"""
function DataFrames.deleterows!(crystal::Crystal, row::Integer)
    deleterows!(crystal.properties, row)
    rows = collect(filter(i -> i ≠ row, 1:nrow(crystal)))
    crystal.positions = crystal.positions[:, rows]
end

function DataFrames.deleterows!(crystal::Crystal, ::Colon)
    deleterows!(crystal, 1:length(crystal))
end

function DataFrames.deleterows!(crystal::Crystal, rows::RowIndices)
    deleterows!(crystal.properties, rows)
    prows = collect(filter(i -> i ∉ rows, 1:nrow(crystal)))
    crystal.positions = crystal.positions[:, prows]
end

"""
    delete!(crystal::Crystal, rows::Integer)
    delete!(crystal::Crystal, rows::AbstractVector{Integer})
    delete!{T <: Integer}(crystal::Crystal, rows::Range{T})
    delete!(crystal::Crystal, rows::Colon)

Alias for `deleterows`[@ref].
"""
delete!(crystal::Crystal, rows::RowIndices) = deleterows!(crystal, rows)

"""
    $(SIGNATURES)

Concatenates crystals together. The lattices must be compatible.
"""
function Base.vcat(crystal::Crystal, other::Vararg{Crystal}; check_periodicity=true)
    if check_periodicity && !are_compatible_lattices(crystal, other...)
        error("Some crystal structures do not have the same periodicity")
    end
    typeof(crystal)(crystal.cell,
                    hcat(crystal.positions,
                         (position_for_crystal(crystal, u) for u in other)...),
                    vcat(crystal.properties, (u.properties for u in other)...))
end
Base.vcat(crystal::Crystal) = crystal

"""
    $(SIGNATURES)

Appends one or more crystal structure to the first structure. Unless `check_periodicity` is
`false`, the structures must have the exact same periodicity. An error will result
otherwise.
"""
function Base.append!(crystal::Crystal, other::Vararg{Crystal}; check_periodicity=true)
    if check_periodicity && !are_compatible_lattices(crystal, other...)
        error("Some crystal structures do not have the same periodicity")
    end
    crystal.positions = hcat(crystal.positions,
                             (position_for_crystal(crystal, u) for u in other)...)
    append!(crystal.properties, (u.properties for u in other)...)
end


"""
    $(SIGNATURES)

Rounds the cell and positions of a crystal. See `round` for possible parameters.
"""
function round!{T, D, U, TT, DD, UU}(crystal::Crystal{T, D, U, Quantity{TT, DD, UU}},
                                     args...)
    crystal.cell = round.(reinterpret(T, crystal.cell), args...) * unit(Quantity{T, D, U})
    const punit = unit(Quantity{TT, DD, UU})
    crystal[:position] = round.(reinterpret(TT, crystal[:position]), args...) * punit
    crystal
end

function round!{T, D, U, TT}(crystal::Crystal{T, D, U, TT}, args...)
    crystal.cell = round.(reinterpret(T, crystal.cell), args...) * unit(Quantity{T, D, U})
    crystal[:position] = round.(crystal[:position], args...)
    crystal
end

"""
    $(SIGNATURES)

Rounds the cell and positions of a crystal. See `Base.round` for possible parameters.

# Examples

```jldoctest
using Crystals
crystal = Crystal([0 0.501 0.501; 0.496 0.001 0.497; 0.497 0.497 0]u"nm",
                  position=[0.001, -0.001, -0.001]u"nm",
                  position=[0.25, 0.251, -0.247]u"nm")
round(crystal, 2)

# output

cell(nm):
  0.0 0.5 0.5
  0.5 0.0 0.5
  0.5 0.5 0.0

│ Atom │ Cartesian           │
├──────┼─────────────────────┤
│ 1    │ (0.0, -0.0, -0.0)   │
│ 2    │ (0.25, 0.25, -0.25) │
```
"""
Base.round(crystal::Crystal, args...) = round!(deepcopy(crystal), args...)

Base.eachindex(crystal::Crystal) = 1:nrow(crystal)

end
