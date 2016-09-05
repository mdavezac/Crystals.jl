module Positions
export Position, PositionArray, PositionDataArray

using Crystals: Log
using FixedSizeArrays: FixedVector, NTuple
using DataFrames: isna, DataArray
import Base
import Base: *

" All acceptable types for positions "
immutable Position{T <:Real, N} <: FixedVector{N, T}
    """ Underlying data """
    _::NTuple{N, T}
end

" Alias to vectors of positions "
typealias PositionArray{T <: Real, N}  Vector{Position{T, N}}
" Alias to data array of positions "
typealias PositionDataArray{T <: Real, N} DataArray{Position{T, N}, 1}

# Add conversion rules from arrays
Base.convert{T <: Real}(::Type{PositionArray{T}}, x::Matrix) =
    convert(Vector{Position{T, size(x, 1)}}, x)
Base.convert{T <: Real, N}(::Type{Vector{Position{T, N}}}, x::Matrix) =
    Position{T, N}[Position{T, N}(x[:, u]) for u in 1:size(x, 2)]
Base.convert{T <: Real, N}(::Type{Array}, x::Vector{Position{T, N}}) =
    T[x[i][j] for j = 1:N, i = 1:length(x)]
function Base.convert(::Type{PositionArray}, x::Matrix)
    if eltype(x) <: Real
        const INNER = eltype(x)
    else
        const reducer = (x, y) -> promote_type(x, typeof(y))
        const INNER = reduce(reducer, typeof(x[1]), x[2:end])
    end
    convert(PositionArray{INNER}, x)
end

function Base.convert{T <: Real, N}(::Type{Array}, x::PositionDataArray{T, N})
    any(isna(x)) &&
        Log.error("Cannot convert DataArray with NA's to desired type")
    T[x[i][j] for j = 1:N, i = 1:length(x)]
end
Base.convert(::Type{PositionDataArray}, x::Matrix) =
    DataArray(convert(PositionArray, x))
Base.convert{T <: Real}(::Type{PositionDataArray{T}}, x::Matrix) =
    DataArray(convert(PositionArray{T}, x))

Base.convert(::Type{PositionArray}, x::Vector) =
    convert(PositionArray, transpose(transpose(x)))
Base.convert{T <: Real}(::Type{PositionArray{T}}, x::Vector) =
    convert(PositionArray{T}, transpose(transpose(x)))
Base.convert(::Type{PositionDataArray}, x::Vector) =
    convert(PositionDataArray, transpose(transpose(x)))
Base.convert{T <: Real}(::Type{PositionDataArray{T}}, x::Vector) =
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
end
