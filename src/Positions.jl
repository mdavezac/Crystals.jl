module Positions
export Position, PositionArray, PositionDataArray

using Crystals: Log
using Unitful: Units, Quantity
import Unitful: ustrip
using FixedSizeArrays: FixedVector, NTuple
using DataFrames: isna, DataArray
using AffineTransforms: AffineTransform
import Base
import Base: *, /, .+, .-

" All acceptable types for positions "
immutable Position{T <:Number, N} <: FixedVector{N, T}
    """ Underlying data """
    _::NTuple{N, T}
end

" Alias to vectors of positions "
typealias PositionArray{T <: Number, N}  Vector{Position{T, N}}
" Alias to data array of positions "
typealias PositionDataArray{T <: Number, N} DataArray{Position{T, N}, 1}

# Add conversion rules from arrays
Base.convert{T <: Number}(::Type{PositionArray{T}}, x::Matrix) =
    convert(Vector{Position{T, size(x, 1)}}, x)
Base.convert{T <: Number, N}(::Type{Vector{Position{T, N}}}, x::Matrix) =
    Position{T, N}[Position{T, N}(x[:, u]) for u in 1:size(x, 2)]
Base.convert{T <: Number, N}(::Type{Array}, x::Vector{Position{T, N}}) =
    T[x[i][j] for j = 1:N, i = 1:length(x)]
function Base.convert(::Type{PositionArray}, x::Matrix)
    if eltype(x) <: Number
        const INNER = eltype(x)
    else
        const reducer = (x, y) -> promote_type(x, typeof(y))
        const INNER = reduce(reducer, typeof(x[1]), x[2:end])
    end
    convert(PositionArray{INNER}, x)
end

function Base.convert{T <: Number, N}(::Type{Array}, x::PositionDataArray{T, N})
    any(isna(x)) &&
        Log.error("Cannot convert DataArray with NA's to desired type")
    T[x[i][j] for j = 1:N, i = 1:length(x)]
end
function Base.convert{T <: Number, N}(::Type{Matrix}, x::PositionDataArray{T, N})
    convert(Array, x)
end
Base.convert(::Type{PositionDataArray}, x::Matrix) =
    DataArray(convert(PositionArray, x))
Base.convert{T <: Number}(::Type{PositionDataArray{T}}, x::Matrix) =
    DataArray(convert(PositionArray{T}, x))

Base.convert(::Type{PositionArray}, x::Vector) =
    convert(PositionArray, transpose(transpose(x)))
Base.convert{T <: Number}(::Type{PositionArray{T}}, x::Vector) =
    convert(PositionArray{T}, transpose(transpose(x)))
Base.convert(::Type{PositionDataArray}, x::Vector) =
    convert(PositionDataArray, transpose(transpose(x)))
Base.convert{T <: Number}(::Type{PositionDataArray{T}}, x::Vector) =
    convert(PositionDataArray{T}, transpose(transpose(x)))

function *{T <: Position}(matrix::Matrix, vectors::Array{T, 1})
    size(matrix, 1) == length(eltype(vectors)) ||
        Log.error("Inconsistent sizes")
    result = similar(vectors)
    for i in 1:length(vectors)
        result[i] = matrix * vectors[i]
    end
    result
end

function *{T <: Position}(matrix::Matrix, vectors::DataArray{T, 1})
    size(matrix, 1) == length(eltype(vectors)) ||
    Log.error("Inconsistent sizes")
    result = similar(vectors)
    for i in 1:length(vectors)
        if !isna(vectors, i)
            result[i] = matrix * vectors[i]
        end
    end
    result
end

*(A::AffineTransform, x::Position) = A * Vector(x)
/(A::AffineTransform, x::Position) = A / Vector(x)

for op in (:+, :-)
    @eval begin
        function Base.$op(p::Position, t::Number)
            const T = promote_type(eltype(p), typeof(t))
            const N = length(p)
            Position{T, N}([$op(v, t) for v in p])
        end
        function Base.$op(p::Position, t::Vector)
            length(p) ≠ length(t) && Log.error("Inconsistent input sizes")
            const T = promote_type(eltype(p), eltype(t))
            const N = length(p)
            Position{T, N}([$op(p[i], t[i]) for i in 1:length(p)])
        end
        function Base.$op{T <: Number, N}(p::PositionArray{T, N}, t::Matrix)
            const NN = length(eltype(p))
            size(t) == (NN, length(p)) ||
                Log.error("Inconsistent Position and matrix sizes")
            const TT = promote_type(eltype(eltype(p)), eltype(t))
            Position{TT, NN}[$op(p[i], t[:, i]) for i in 1:length(p)]
        end
        Base.$op{T <: Number, N}(p::PositionDataArray{T, N}, t::Matrix) =
            DataArray($op(p.data, t), p.na)
    end
end
for (dotop, op) in ((:.+, :+), (:.-, :-))
    @eval begin
        function $dotop{T <: Number, N}(p::PositionArray{T, N}, t::Position)
            length(eltype(p)) == length(t) || Log.error("Inconsistent sizes")
            const TT = promote_type(eltype(eltype(p)), eltype(t))
            const NN = length(t)
            Position{TT, NN}[$op(v, t) for v in p]
        end
        function $dotop{T <: Number, N}(p::PositionArray{T, N}, t::Vector)
            const TT = promote_type(eltype(eltype(p)), eltype(t))
            $dotop(p, Position{TT}(t))
        end
        Base.$dotop{T <: Number, N}(p::PositionDataArray{T, N}, t::Vector) =
            DataArray($dotop(p.data, t), p.na)
        Base.$op{T <: Number, N}(p::PositionDataArray{T, N}, t::Vector) =
            DataArray($dotop(p.data, t), p.na)
        Base.$dotop{T <: Number, N}(p::PositionDataArray{T, N}, t::Position) =
            DataArray($dotop(p.data, t), p.na)
        Base.$op{T <: Number, N}(p::PositionDataArray{T, N}, t::Position) =
            DataArray($dotop(p.data, t), p.na)
    end
end

# Makes it easier to create dimensionalize position
# see arraymath.jl/Unitful.jl for similar code
for op in (:.*, :*)
    @eval begin
        ($op)(A::Units, B::Position) = Position((A * u for u in B)...)
        ($op)(B::Position, A::Units) = Position((A * u for u in B)...)
    end
end

"""
is_fractional(pos::Position)

True if the positions is given in fractional coordinates.
A position is in fractional coordinates if it is dimensionless
"""
is_fractional(pos::Position) = !(eltype(position) <: Quantity)

ustrip(p::Position) = Position((ustrip(u) for u in p)...)

end
