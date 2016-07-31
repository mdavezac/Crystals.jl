for s in 2:6
  local name = Symbol("Position$(s)D")
  local docstring = "Position for $s-dimensional crystals"
  local code = begin
    u = quote
      $docstring
      immutable $name{T <:Real} <: FixedVectorNoTuple{$s, T}
      end
    end
    for (i, symb) in enumerate([:x, :y, :z, :u, :v, :w][1:s])
      push!(u.args[4].args[3].args, :($symb::T))
    end
    u
  end
  eval(code)
end

" All acceptable types for positions "
typealias PositionTypes Union{
  Position2D, Position3D,
  Position4D, Position5D, Position6D
}

const MAX_POSITION_SIZE = maximum([length(u) for u in PositionTypes.types])
const MIN_POSITION_SIZE = minimum([length(u) for u in PositionTypes.types])
function positions_type(x::Vector)
  MIN_POSITION_SIZE ≤ length(x) ≤ MAX_POSITION_SIZE ||
    error("Incorrect column vector size")
  position_index = findfirst(u -> length(u) == length(x), PositionTypes.types)
  PositionTypes.types[position_index]
end

# Add conversion rules from arrays
for position_type in PositionTypes.types
  @eval begin
    function Base.convert{T <: Real}(::Type{Vector{$(position_type){T}}}, x::Matrix)
      size(x, 1) == $(length(position_type)) || error("Incompatible vector size")
      $(position_type){T}[$(position_type)(x[:, u]) for u in 1:size(x, 2)]
    end
  end
end

" Converts inputs to something a DataArray of Positions will understand "
function Positions(x::Matrix)
  MIN_POSITION_SIZE ≤ size(x, 1) ≤ MAX_POSITION_SIZE ||
    error("Incorrect column vector size")
  return convert(Vector{positions_type(x[:, 1]){eltype(x)}}, x)
end
Positions(x::DataArray) = x
function Positions(x::Vector)
  T = positions_type(x){eltype(x)}
  T[T(x)]
end

function Position(x::Vector)
  T = positions_type(x){eltype(x)}
  convert(T, x)
end
function Position{T <: Real}(x::T...)
  position_index = findfirst(u -> length(u) == length(x), PositionTypes.types)
  PositionTypes.types[position_index](x...)
end

abstract AbstractCrystal

"""
Crystal structure

Holds all data necessary to define a crystal structure.
The cell and scale are held as fields.
The atomic properties, including the positions, are held in a DataFrame.
"""
type Crystal{T} <: AbstractCrystal
  """ Scale/unit of the positions and cell """
  scale::T
  """ Periodicity of the crystal structure """
  cell::Matrix{T}
  """ Atoms and atomic properties """
  atoms::DataFrame

  function Crystal(cell, scale, atoms::DataFrame)
    size(cell, 1) == size(cell, 2) || error("Cell matrix is not square")
    MIN_POSITION_SIZE ≤ size(cell, 1) ≤ MAX_POSITION_SIZE ||
      error("Incorrect column vector size")
    :position ∈ names(atoms) || nrow(atoms) == 0 ||
      error("Input Dataframe has atoms without positions")
    if nrow(atoms) == 0
      atoms[:position] = Vector{positions_type(cell[:, 1]){eltype(cell)}}()
    else
      eltype(atoms[:position]) <: PositionTypes ||
        error("Positions do not have acceptable type")
      size(cell, 1) == length(eltype(atoms[:position])) ||
        error("Cell and position dimensionality do not match")
    end
    new(scale, cell, atoms)
  end
end

function Crystal(cell::Matrix, scale=1::Real; kwargs...)
  nkwargs = Tuple{Symbol, Any}[]
  for i in 1:length(kwargs)
    if kwargs[i][1] ≠ :position
      push!(nkwargs, kwargs[i])
      continue
    end
    position = Positions(kwargs[i][2])
    size(cell, 1) == length(eltype(position)) ||
      error("Dimensionality of cell and positions do not match")
    push!(nkwargs, (:position, position))
  end

  atoms = DataFrame(; nkwargs...)
  Crystal{eltype(cell)}(cell, scale, atoms)
end

function Crystal(cell::Matrix, columns::Vector{Any}, names::Vector{Symbol},
                 scale=1::Real)
  for (i, (value, name)) in enumerate(zip(columns[:], names[:]))
    name == :position || continue

    columns[i] = Positions(value)
    size(cell, 1) == length(eltype(columns[i])) ||
      error("Dimensionality of cell and positions do not match")
  end

  atoms = DataFrame(columns, names)
  Crystal{eltype(cell)}(cell, scale, atoms)
end

# Forwards all indexing to the DataFrame
""" Equivalent to `getindex(crystal.atoms, ...)` """
Base.getindex(crystal::Crystal, col::Any) = crystal.atoms[col]
Base.getindex(crystal::Crystal, col::Any, row::Any) = crystal.atoms[col, row]
Base.names(crystal::Crystal) = names(crystal.atoms)
""" Equivalent to `endof(crystal.atoms)` """
Base.endof(crystal::Crystal) = endof(crystal.atoms)
""" Equivalent to `size(crystal.atoms)` """
Base.size(crystal::Crystal) = size(crystal.atoms)
Base.size(crystal::Crystal, i) = size(crystal.atoms, i)
""" Equivalent to `ndims(crystal.atoms)` """
Base.ndims(crystal::Crystal) = ndims(crystal.atoms)

# Single column
Base.setindex!(crystal::Crystal, v::Crystal, col::Any) =
  setindex!(crystal.atoms, v.atoms, col)
Base.setindex!(crystal::Crystal, v::Matrix, col::Any) =
  setindex!(crystal.atoms, Positions(v), col)
Base.setindex!(crystal::Crystal, v::Any, col::Any) =
  setindex!(crystal.atoms, v, col)
Base.setindex!(crystal::Crystal, v::Crystal, row::Any, col::Any) =
  setindex!(crystal.atoms, v.atoms, row, col)
Base.setindex!(crystal::Crystal, v::Matrix, row::Any, col::Any) =
  setindex!(crystal.atoms, Positions(v), row, col)
Base.setindex!(crystal::Crystal, v::Any, row::Any, col::Any) =
  setindex!(crystal.atoms, v, row, col)

# Special deletion assignment
Base.setindex!(crystal::Crystal, x::Void, col::Int) = delete!(crystal, col)
Base.empty!(crystal::Crystal) = (empty!(crystal.atoms); crystal)
Base.insert!(crystal::Crystal, col::Int, item::AbstractVector, name::Symbol) =
  (insert!(crystal.atoms, col, item, name); crystal)

Base.merge!(crystal::Crystal, others::DataFrame...) =
  (merge!(crystal.atoms, others...); crystal)
Base.copy!(crystal::Crystal) =
  Crystal{eltype(crystal.cell)}(
    copy(crystal.cell), crystal.scale, copy(crystal.atoms))
Base.deepcopy(crystal::Crystal) =
  Crystal{eltype(crystal.cell)}(
    deepcopy(crystal.cell), crystal.scale, deepcopy(crystal.atoms))

Base.delete!(crystal::Crystal, cols::Any) =
  (delete!(crystal.atoms, cols); crystal)
deleterows!(crystal::Crystal, cols::Any) =
  (deleterows!(crystal.atoms, cols); crystal)
hcat!(crystal::Crystal, x...) = (hcat!(crystal.atoms, x...); crystal)
Base.hcat(crystal::Crystal, x...) = hcat!(copy(crystal), x...)
nullable!(crystal::Crystal, x...) = (nullable!(crystal.atoms, x...); crystal)
pool!(crystal::Crystal, x::Any) = pool!(crystal.atoms, x::Any)
Base.append!(crystal::Crystal, atoms::DataFrame) =
  (append!(crystal.atoms, atoms); crystal)
Base.push!(crystal::Crystal, x) = push!(crystal.atoms, x)

function Base.push!(crystal::Crystal, position::PositionTypes; kwargs...)
  if length(position) ≠ size(crystal.cell, 1)
    error("Dimensionality of input position and crystal do not match")
  end
  row = Any[NA for u in 1:length(crystal.atoms)]
  row[index(crystal.atoms)[:position]] = position
  missing = Tuple{Symbol, Any}[]
  const colnames = names(crystal)
  for (name, value) in kwargs
    if name ∉ colnames
      push!(missing, (name, value))
    else
      row[index(crystal.atoms)[name]] = value
    end
  end
  push!(crystal, row)
  for (name, value) in missing
    # given twice, makes no sense
    @assert name ∉ names(crystal)
    crystal[name] = fill(value, size(crystal, 1))
    crystal[1:(size(crystal, 1) - 1), name] = NA
  end
end

function Base.push!{T <: Real}(crystal::Crystal, position::Vector{T}; kwargs...)
  pos = positions_type(position){eltype(crystal.cell)}(position)
  Base.push!(crystal, pos; kwargs...)
end

function Base.show(io::IO, crystal::Crystal, args...; kwargs...)
  println(io, typeof(crystal))
  println(io, "cell: ", crystal.cell)
  println(io, "scale: ", crystal.scale)
  println(io, crystal.atoms)
end


function Base.showcompact(io::IO, pos::PositionTypes)
  result = string(pos)
  print(io, result[findfirst(result, '('):end])
end

function ourshowcompact(io::IO, pos::PositionTypes)
  result = string(pos)
  print(io, result[findfirst(result, '(') + 1:end - 1])
end
