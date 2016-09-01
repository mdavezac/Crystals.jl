module Positions
export Position, PositionArray, PositionDataArray

using FixedSizeArrays: FixedVectorNoTuple
using DataFrames: isna, DataArray
import Base

" All acceptable types for positions "
abstract Position{T <:Real, N} <: FixedVectorNoTuple{N, T}
const MAX_POSITION_SIZE = 6
const MIN_POSITION_SIZE = 2

for s in MIN_POSITION_SIZE:MAX_POSITION_SIZE
    local name = Symbol("Position$(s)D")
    local docstring = "Position for $s-dimensional crystals"
    local code = begin
        u = quote
            $docstring
            immutable $name{T <:Real} <: Position{T, $s}
            end
        end
        for (i, symb) in enumerate([:x, :y, :z, :u, :v, :w][1:s])
            push!(u.args[4].args[3].args, :($symb::T))
        end
        u
    end
    eval(code)
end

" All derived position types "
const PositionTypes = ([
    eval(parse("Position$(n)D"))
    for n in MIN_POSITION_SIZE:MAX_POSITION_SIZE
]...)

" Alias to vectors of positions "
typealias PositionArray{T <: Real, N}  Vector{Position{T, N}}
" Alias to data array of positions "
typealias PositionDataArray{T <: Real, N} DataArray{Position{T, N}, 1}

# Add conversion rules from arrays
Base.convert{T <: Position}(::Type{Vector{T}}, x::Matrix) =
    T[T(x[:, u]) for u in 1:size(x, 2)]
Base.convert{T <: Position}(::Type{Array}, x::Vector{T}) =
    eltype(eltype(x))[x[i][j] for j = 1:length(T), i = 1:length(x)]
function Base.convert{T <: Position}(::Type{Array}, x::DataArray{T, 1})
    any(isna(x)) && error("Cannot convert DataArray with NA's to desired type")
    eltype(eltype(x))[x[i][j] for j = 1:length(T), i = 1:length(x)]
end
function Base.convert{T <: Position}(::Type{Vector{T}}, x::Matrix)
    length(T) == size(x, 1) || error("Columns cannot be converted to one of $T")
    T[convert(T, x[:, u]) for u in 1:size(x, 2)]
end
function Base.convert{T <: Real}(::Type{Position}, x::Vector{T})
    if eltype(x) <: Real
        const INNER = eltype(x)
    else
        const reducer = (x, y) -> promote_type(x, typeof(y))
        const INNER = reduce(reducer, typeof(x[1]), x[2:end])
    end
    convert(Position{INNER}, x)
end
function Base.convert{T <: Real}(::Type{Position{T}}, x::Vector)
    k = findfirst(X -> length(X) == size(x, 1), PositionTypes)
    k ≠ 0 || error("Cannot create position of size $(size(x, 1))")
    convert(PositionTypes[k]{T}, x)
end
function Base.convert(::Type{PositionArray}, x::Matrix)
    if eltype(x) <: Real
        const INNER = eltype(x)
    else
        const reducer = (x, y) -> promote_type(x, typeof(y))
        const INNER = reduce(reducer, typeof(x[1]), x[2:end])
    end
    convert(PositionArray{INNER}, x)
end
function Base.convert{T <: Real}(::Type{PositionArray{T}}, x::Matrix)
    k = findfirst(X -> length(X) == size(x, 1), PositionTypes)
    k ≠ 0 || error("Cannot create position of size $(size(x, 1))")
    convert(Vector{PositionTypes[k]{T}}, x)
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
end
